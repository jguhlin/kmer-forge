use crossbeam::channel::{Receiver, Sender, bounded, unbounded};
use dashmap::DashMap;
use needletail::{Sequence, parse_fastx_reader};
use pulp::Arch;
use xxhash_rust::xxh3::xxh3_64;
use rayon::prelude::*;

use std::borrow::Cow;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};
use std::thread::{self, JoinHandle};

use crate::rolling_encoder::RollingKmer3;

static DNA_CODES: [u8; 256] = {
    const X: u8 = 4;
    let mut table = [X; 256];
    table[b'A' as usize] = 0; table[b'a' as usize] = 0;
    table[b'C' as usize] = 1; table[b'c' as usize] = 1;
    table[b'G' as usize] = 2; table[b'g' as usize] = 2;
    table[b'T' as usize] = 3; table[b't' as usize] = 3;
    // everything else stays 4
    table
};

#[inline(always)]
fn encode_base(base: u8) -> u8 {
    DNA_CODES[base as usize]
}

#[inline(always)]
fn is_bad(mmer: &[u8]) -> bool {
    // skip if any base == 4 ('N')
    if mmer.contains(&4) {
        return true;
    }
    if mmer.len() >= 3 {
        // skip if prefix is AAA(0,0,0) or ACA(0,1,0)
        if &mmer[..3] == [0,0,0] || &mmer[..3] == [0,1,0] {
            return true;
        }
    }
    // skip if interior 'AA' -> [0,0] except at index 0
    // i.e. search for [0,0] starting at i>=1
    for i in 1..mmer.len() {
        if i < mmer.len() && mmer[i-1] == 0 && mmer[i] == 0 {
            // found [0,0] at index i-1..i
            // skip if i-1 > 0, i.e. i>1
            if i-1 > 0 {
                return true;
            }
        }
    }
    false
}

#[inline(always)]
fn canonical_mmer<'a>(fwd: &'a [u8], rev: &'a [u8]) -> &'a [u8] {
    if rev < fwd { rev } else { fwd }
}

fn compute_signature(kmer: &[u8], m: usize) -> Vec<u8> {
    let k_len = kmer.len();

    // 1) encode entire k-mer
    // For big k, do it with pulp in chunks:
    let arch = Arch::new();
    let mut encoded = vec![0u8; k_len];
    arch.dispatch(|| {
        for (i, &b) in kmer.iter().enumerate() {
            encoded[i] = encode_base(b);
        }
    });

    // 2) build reverse complement once
    let mut encoded_rc = vec![0u8; k_len];
    arch.dispatch(|| {
        let end = k_len - 1;
        for i in 0..k_len {
            let c = encoded[i];
            // complement
            let c_rc = match c {
                0 => 3,
                1 => 2,
                2 => 1,
                3 => 0,
                _ => 4, // N->N
            };
            encoded_rc[end - i] = c_rc;
        }
    });

    // 3) find lexicographically smallest “good” m-mer
    let mut best: Option<&[u8]> = None;

    for start in 0..=(k_len - m) {
        let fwd = &encoded[start..(start + m)];
        // corresponding rc slice
        //   forward window is start..start+m
        //   in rc => (k_len - (start+m))..(k_len - start)
        let rc_start = k_len - (start + m);
        let rev = &encoded_rc[rc_start..(rc_start + m)];

        let cand = canonical_mmer(fwd, rev);
        if is_bad(cand) {
            continue;
        }

        match best {
            Some(current) => {
                if cand < current {
                    best = Some(cand);
                }
            }
            None => {
                best = Some(cand);
            }
        }
    }

    best.unwrap().to_vec()
}

fn read_to_superkmers(read_seq: &[u8], k: usize, m: usize) -> Vec<(Vec<u8>, Vec<u8>)> {
    let mut result = Vec::new();
    if read_seq.len() < k {
        return result;
    }

    let first_sig = compute_signature(&read_seq[0..k], m);
    let mut current_sig = first_sig;
    let mut start_pos = 0;

    for i in 1..=(read_seq.len() - k) {
        let next_sig = compute_signature(&read_seq[i..i + k], m);
        if next_sig != current_sig {
            // close out
            let block = read_seq[start_pos..(i + k - 1)].to_vec();
            result.push((current_sig, block));

            // start new
            current_sig = next_sig;
            start_pos = i;
        }
    }
    // final
    let block = read_seq[start_pos..].to_vec();
    result.push((current_sig, block));
    result
}

fn expand_superkmer_to_kmers(superkm_seq: &[u8], k: usize) -> impl Iterator<Item = &[u8]> {
    // yields each k-mer slice
    superkm_seq.windows(k)
}

pub struct KmerCounter {
    k: u8,
    temp_path: String,
    bins: Arc<Vec<KmerBin>>,
    bin_count: u16,
    threads: usize,
    workers: Vec<JoinHandle<()>>,
    kmer_tx: Sender<Vec<(Vec<u8>, Vec<u8>)>>,
    kmer_rx: Receiver<Vec<(Vec<u8>, Vec<u8>)>>,
    compression_tx: Sender<(usize, Vec<Vec<u8>>)>,
    compression_rx: Receiver<(usize, Vec<Vec<u8>>)>,
    output_tx: Sender<(usize, Vec<u8>)>,
    output_rx: Receiver<(usize, Vec<u8>)>,
    shutdown_flag: Arc<AtomicBool>,
}

impl KmerCounter {
    pub fn new(
        k: u8,
        temp_path: String,
        threads: usize,
        buffer_flush_size: usize,
        bin_power: u8,
    ) -> Self {
        assert!(k < 32, "Kmer size must be less than 32");
        let bin_count: usize = 2_usize.pow(bin_power as u32);
        assert!(bin_count > 0, "Bin count must be greater than 0");
        assert!(
            bin_count < u16::MAX as usize,
            "Bin count must be less than u16::MAX"
        );

        let mut bins = Vec::with_capacity(bin_count as usize);

        // Confirm temp_path is writable, and create the directory if it doesn't exist
        let temp_path = std::path::Path::new(&temp_path);
        if !temp_path.exists() {
            std::fs::create_dir(temp_path).expect("Could not create temp directory");
        }

        // Confirm the temp_path is empty
        let temp_path = temp_path.to_str().unwrap();
        let temp_path = format!("{}/", temp_path);
        let temp_path = std::path::Path::new(&temp_path);
        if temp_path.read_dir().unwrap().count() > 0 {
            panic!("Temp directory is not empty");
        }

        // Create the bins
        for i in 0..bin_count {
            let bin_path = format!(
                "{}/bin_{}.kmer_forge_temp_bin",
                temp_path.to_str().unwrap(),
                i
            );
            let out_fh =
                BufWriter::new(File::create(bin_path.clone()).expect("Could not create bin file"));
            let out_fh = Mutex::new(out_fh);
            bins.push(KmerBin {
                number: i as u16,
                filename: bin_path,
                out_fh,
                buffer: Mutex::new(Vec::with_capacity(buffer_flush_size)),
                buffer_flush_size,
            });
        }

        let bins = Arc::new(bins);

        // Create the channels
        let (kmer_tx, kmer_rx)= bounded(64);
        let (compression_tx, compression_rx) = bounded(64);
        let (output_tx, output_rx) = unbounded();

        // Create the workers
        let shutdown_flag = Arc::new(AtomicBool::new(false));
        let mut workers = Vec::with_capacity(threads);
        for _ in 0..threads {
            let kmer_rx = kmer_rx.clone();
            let compression_tx = compression_tx.clone();
            let compression_rx = compression_rx.clone();
            let output_tx = output_tx.clone();
            let output_rx = output_rx.clone();
            let shutdown_flag = shutdown_flag.clone();
            let bins = bins.clone();
            let worker = thread::spawn(move || {
                kmer_worker(
                    kmer_rx,
                    compression_rx,
                    compression_tx,
                    output_rx,
                    output_tx,
                    shutdown_flag,
                    bins,
                    bin_power,
                );
            });
            workers.push(worker);
        }

        KmerCounter {
            k,
            temp_path: temp_path.to_str().unwrap().to_string(), // stupid....
            bins,
            bin_count: bin_count as u16,
            threads,
            shutdown_flag,
            workers,
            kmer_tx,
            kmer_rx,
            compression_tx,
            compression_rx,
            output_tx,
            output_rx,
        }
    }

    pub fn submit(&self, kmers: Vec<(Vec<u8>, Vec<u8>)>) {
        self.kmer_tx
            .send(kmers)
            .expect("Could not send kmers to worker");
    }

    pub fn try_submit(
        &self,
        kmers: Vec<(Vec<u8>, Vec<u8>)>,
    ) -> Result<(), crossbeam::channel::TrySendError<Vec<(Vec<u8>, Vec<u8>)>>> {
        self.kmer_tx.try_send(kmers)
    }

    pub fn stop_gathering(&mut self) {
        let backoff = crossbeam::utils::Backoff::new();
        while !self.kmer_rx.is_empty() {
            if backoff.is_completed() {
                // Sleep for 100ms
                std::thread::sleep(std::time::Duration::from_millis(100));
                backoff.reset();
            } else {
                backoff.snooze();
            }
        }

        // All kmers are processed, clear the bins
        for bin in self.bins.iter() {
            let mut bin_lock = bin.buffer.lock().unwrap();
            if !bin_lock.is_empty() {
                self.compression_tx
                    .send((bin.number as usize, bin_lock.clone()))
                    .expect("Could not send buffer to compressor");
            }
            bin_lock.clear();
            bin_lock.shrink_to_fit();
        }

        // Wait for all workers to finish
        while !self.compression_rx.is_empty() || !self.output_rx.is_empty() {
            if backoff.is_completed() {
                // Sleep for 100ms
                std::thread::sleep(std::time::Duration::from_millis(100));
                backoff.reset();
            } else {
                backoff.snooze();
            }
        }

        self.shutdown_flag.store(true, Ordering::Relaxed);
        for worker in self.workers.drain(..) {
            worker.join().expect("Could not join worker thread");
        }
    }

    pub fn merge_bins(self) {
        // Drain the bins, destruct, close the out_fh and open the file for reading instead

        let KmerCounter {
            temp_path,
            bins,
            threads,
            ..
        } = self;

        let mut bins = Arc::into_inner(bins).expect("Could not get bins");

        let mut counts: DashMap<u64, u32> = DashMap::new();

        let bins = bins.drain(..).into_iter().map(|x| (x.number, x.filename)).collect::<Vec<_>>();

        // Using rayon, merge the bins
        bins.into_par_iter().for_each(|(number, filename)| {
            let bin = File::open(filename).expect("Could not open bin file");
            let mut reader = BufReader::new(bin);

            let bincode_config = bincode::config::standard().with_fixed_int_encoding();

            let mut decompressor =
                zstd::bulk::Decompressor::new().expect("Could not create decompressor");

            loop {
                let kmers: Vec<u8> =
                    match bincode::decode_from_std_read(&mut reader, bincode_config) {
                        Ok(kmers) => kmers,
                        Err(_e) =>
                        // EOF most likely
                        {
                            break;
                        }
                    };

                if kmers.is_empty() {
                    break;
                }

                // Capacity set to 8Gb, should never be that high, ofc....
                let kmers = decompressor
                    .decompress(&kmers, 8 * 1024 * 1024 * 1024)
                    .expect("Could not decompress buffer");
                let kmers: Vec<u64> = bincode::decode_from_slice(&kmers, bincode_config)
                    .expect("Could not decode buffer")
                    .0;

                for kmer in kmers {
                    if counts.contains_key(&kmer) {
                        counts.alter(&kmer, |_, count| {
                            count.saturating_add(1)
                        });
                    } else {
                        counts.insert(kmer, 1);
                    }                    
                }
            }});


        println!("Counts: {:?}", counts.len());
        // Convert dashmap into a hashmap
        let counts: HashMap<u64, u32> = counts.into_iter().map(|(k, v)| (k, v)).collect();

        // Save to file
        let mut out_fh = BufWriter::new(
            File::create(format!("{}/counts.bin", temp_path)).expect("Could not create counts file"),
        );

        bincode::encode_into_std_write(
            counts,
            &mut out_fh,
            bincode::config::standard().with_fixed_int_encoding(),
        ).expect("Could not write to counts file");



    }
}

fn kmer_worker(
    kmer_rx: crossbeam::channel::Receiver<Vec<(Vec<u8>, Vec<u8>)>>,
    compression_rx: crossbeam::channel::Receiver<(usize, Vec<Vec<u8>>)>,
    compression_tx: crossbeam::channel::Sender<(usize, Vec<Vec<u8>>)>,
    output_rx: crossbeam::channel::Receiver<(usize, Vec<u8>)>,
    output_tx: crossbeam::channel::Sender<(usize, Vec<u8>)>,
    shutdown_flag: Arc<AtomicBool>,
    bins: Arc<Vec<KmerBin>>,
    bin_power: u8,
) {
    let bin_mask = (1 << bin_power) - 1;
    const FLUSH_THRESHOLD: usize = 16384;

    let mut compressor = zstd::bulk::Compressor::new(-1).expect("Could not create compressor");

    // did not seem to help
    // compressor.set_parameter(zstd::stream::raw::CParameter::Strategy(zstd::zstd_safe::zstd_sys::ZSTD_strategy::ZSTD_fast)).expect("Could not set compression level");

    // Create a thread-local buffer for each bin.
    let bin_count = bins.len();
    let mut local_buffers: Vec<Vec<Vec<u8>>> = (0..bin_count)
        .map(|_| Vec::with_capacity(FLUSH_THRESHOLD))
        .collect();

    let mut bins_to_submit = Vec::with_capacity(128);

    loop {
        if shutdown_flag.load(Ordering::Relaxed) {
            break;
        }

        while let Ok((bin, compressed)) = output_rx.try_recv() {
            let mut bin_lock = bins[bin].out_fh.lock().unwrap();
            bincode::encode_into_std_write(
                compressed,
                &mut *bin_lock,
                bincode::config::standard().with_fixed_int_encoding(),
            )
            .expect("Could not write to bin file");
            drop(bin_lock);
        }

        if !output_rx.is_full() {
            while let Ok((bin, mut kmers)) = compression_rx.try_recv() {
                kmers.sort_unstable();

                let encoded = bincode::encode_to_vec(
                    &kmers,
                    bincode::config::standard().with_fixed_int_encoding(),
                )
                .expect("Could not write to bin file");

                // zstd
                let compressed = compressor
                    .compress(&encoded)
                    .expect("Could not compress buffer");

                // lz4_flex
                // let compressed = compress(&encoded);

                // no compression
                // No real speed difference...
                // let compressed = encoded;

                output_tx
                    .send((bin, compressed))
                    .expect("Could not send compressed buffer to flusher");
            }
        }

        let kmers = match kmer_rx.try_recv() {
            Err(crossbeam::channel::TryRecvError::Empty) => continue,
            Err(crossbeam::channel::TryRecvError::Disconnected) => break,
            Ok(kmers) => kmers,
        };

        // For each kmer, calculate the bin and store it in the corresponding local buffer.
        for (minimzer, superkmer) in kmers {
            let bin_index = minimzer[0] as usize & bin_mask;

            local_buffers[bin_index].push(superkmer);

            // If the thread-local buffer has reached the threshold, flush it.
            if local_buffers[bin_index].len() >= FLUSH_THRESHOLD {
                // Acquire the lock for the global bin buffer.
                let mut global_buffer = bins[bin_index]
                    .buffer
                    .lock()
                    .expect("Could not acquire bin lock");
                // Move the local buffer's contents into the global buffer.
                global_buffer.append(&mut local_buffers[bin_index]);
                // Optionally, mark for flush if the global buffer exceeds its flush size.
                if global_buffer.len() > bins[bin_index].buffer_flush_size {
                    bins_to_submit.push(bin_index);
                }
            }
        }

        if !bins_to_submit.is_empty()
            && compression_rx.len() as f32 <= 0.8 * compression_rx.capacity().unwrap() as f32
        {
            for bin in bins_to_submit.drain(..) {
                let mut bin_lock = match bins[bin].buffer.try_lock() {
                    Ok(lock) => lock,
                    Err(_) => continue, // Assume another thread is flushing this buffer
                };

                // Confirm that another thread didn't flush it
                if bin_lock.len() < bins[bin].buffer_flush_size {
                    continue;
                }

                let mut bin_buffer = Vec::with_capacity(bin_lock.len());
                std::mem::swap(&mut *bin_lock, &mut bin_buffer);
                drop(bin_lock);

                compression_tx
                    .send((bin, bin_buffer))
                    .expect("Could not send buffer to compressor");
            }
        }

        bins_to_submit.clear();
    }

    // Final flush: after shutdown, flush any remaining items from the thread-local buffers.
    for (bin_index, local_buf) in local_buffers.iter_mut().enumerate() {
        if !local_buf.is_empty() {
            let mut global_buffer = bins[bin_index]
                .buffer
                .lock()
                .expect("Could not acquire bin lock");
            global_buffer.append(local_buf);
        }
    }
}

pub struct KmerBin {
    number: u16,
    out_fh: Mutex<BufWriter<std::fs::File>>,
    filename: String,
    buffer: Mutex<Vec<Vec<u8>>>,
    buffer_flush_size: usize,
}

// kmers up to 31 bases long
pub fn count_kmers_file(kmer_counter: &mut KmerCounter, file: &str, k: u8, min_quality: u8) {
    let file = File::open(file).expect("Could not open file");
    let reader = BufReader::with_capacity(8 * 1024 * 1024, file);
    let mut reader = parse_fastx_reader(reader).expect("Invalid file");

    // debugging
    let mut processed_reads = 0;

    let mut kmers_to_submit = Vec::with_capacity(8 * 1024);

    while let Some(record) = reader.next() {
        processed_reads += 1;
        if processed_reads % 100000 == 0 {
            println!("Processed {} reads", processed_reads);
        }
        let record = record.expect("Error reading record");
        let mut seq = record.seq();
        let qual = record.qual();

        let arch = Arch::new();

        if let Some(qual) = qual {
            // If average is less than min_quality, skip
            let total_qual: u64 = qual.iter().map(|q| *q as u64).sum();
            if (total_qual / qual.len() as u64) < min_quality as u64 {
                continue;
            }

            // Otherwise mask sequence when quality is less than min_quality

            let masked_seq: Vec<u8> = arch.dispatch(|| {
                seq.iter()
                    .zip(qual.iter())
                    .map(|(base, q)| if q < &min_quality { b'N' } else { *base })
                    .collect()
            });

            seq = Cow::Owned(masked_seq);
        } // If no quality, we don't worry about it

        // Is seq long enough to have a kmer?
        if seq.len() < k as usize {
            continue;
        }

        let seq = seq.strip_returns();
        let seq = seq.normalize(true);

        kmers_to_submit.extend(read_to_superkmers(&seq, k as usize, 7).drain(..));

        if kmers_to_submit.len() >= 8 * 1024 {
            match kmer_counter.try_submit(kmers_to_submit) {
                Ok(_) => {
                    kmers_to_submit = Vec::with_capacity(8 * 1024);
                }
                Err(crossbeam::channel::TrySendError::Full(kmers)) => {
                    kmers_to_submit = kmers;
                } // Process more reads then...
                Err(crossbeam::channel::TrySendError::Disconnected(_kmers)) => {
                    panic!("Kmer counter disconnected before shutdown");
                }
            }
        }
    }
}
