use bytes::{Buf, BufMut, Bytes};
use crossbeam::channel::{Receiver, Sender, bounded, unbounded};
use dashmap::DashMap;
use needletail::{kmer, parse_fastx_reader, Sequence};
use nthash::NtHashIterator;
use pulp::Arch;
use rayon::prelude::*;
use simd_minimizers::one_minimizer;
use simd_minimizers::packed_seq::AsciiSeq;
use xxhash_rust::xxh3::xxh3_64;

use std::borrow::Cow;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};
use std::thread::{self, JoinHandle};

use crate::SuperKmerStorage;

const LOCAL_FLUSH_THRESHOLD: usize = 8 * 1024;
const GLOBAL_FLUSH_THRESHOLD: usize = 32 * 1024;

// TODO: This is where I'm leaving off
// Need a way to store kmers as offsets rather than the actual sequence
// Maybe precompute bins and store the bin number and offset in the bin as well?
pub type SuperKmer = (u64, Vec<u8>, Vec<u8>);

pub struct KmerCounter {
    k: u8,
    min_quality: u8,
    temp_path: String,
    bins: Arc<Vec<KmerBin>>,
    bin_count: u16,
    threads: usize,
    workers: Vec<JoinHandle<()>>,
    reads_tx: Sender<Vec<(Bytes, Option<Bytes>)>>,
    reads_rx: Receiver<Vec<(Bytes, Option<Bytes>)>>,
    kmer_tx: Sender<Vec<(u64, Bytes)>>,
    kmer_rx: Receiver<Vec<(u64, Bytes)>>,
    compression_tx: Sender<(usize, SuperKmerStorage)>,
    compression_rx: Receiver<(usize, SuperKmerStorage)>,
    output_tx: Sender<(usize, Vec<u8>)>,
    output_rx: Receiver<(usize, Vec<u8>)>,
    shutdown_flag: Arc<AtomicBool>,
}

impl KmerCounter {
    pub fn new(k: u8, temp_path: String, threads: usize, bin_power: u8, min_quality: u8) -> Self {
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
                buffer: Mutex::new(SuperKmerStorage::with_capacity(
                    GLOBAL_FLUSH_THRESHOLD,
                    k as usize,
                )),
            });
        }

        let bins = Arc::new(bins);

        // Create the channels
        let (reads_tx, reads_rx) = bounded(threads + 4);
        let (kmer_tx, kmer_rx) = bounded(threads + 4);
        let (compression_tx, compression_rx) = bounded(threads + 4);
        let (output_tx, output_rx): (Sender<(usize, Vec<u8>)>, Receiver<(usize, Vec<u8>)>) =
            bounded(threads + 4);

        // Create the workers
        let shutdown_flag = Arc::new(AtomicBool::new(false));
        let mut workers = Vec::with_capacity(threads);
        for _ in 0..threads {
            let reads_rx = reads_rx.clone();
            let kmer_rx = kmer_rx.clone();
            let kmer_tx = kmer_tx.clone();
            let compression_tx = compression_tx.clone();
            let compression_rx = compression_rx.clone();
            let output_tx = output_tx.clone();
            let output_rx = output_rx.clone();
            let shutdown_flag = shutdown_flag.clone();
            let bins = bins.clone();
            let worker = thread::spawn(move || {
                kmer_worker(
                    reads_rx,
                    kmer_tx,
                    kmer_rx,
                    compression_rx,
                    compression_tx,
                    output_rx,
                    output_tx,
                    shutdown_flag,
                    bins,
                    bin_power,
                    k,
                    min_quality,
                    threads,
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
            reads_tx,
            reads_rx,
            min_quality,
        }
    }

    pub fn submit(&self, kmers: Vec<(Bytes, Option<Bytes>)>) {
        self.reads_tx
            .send(kmers)
            .expect("Could not send kmers to worker");
    }

    pub fn try_submit(
        &self,
        kmers: Vec<(u64, Bytes)>,
    ) -> Result<(), crossbeam::channel::TrySendError<Vec<(u64, Bytes)>>> {
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
                let mut new_storage = SuperKmerStorage::new();
                std::mem::swap(&mut *bin_lock, &mut new_storage);
                self.compression_tx
                    .send((bin.number as usize, new_storage))
                    .expect("Could not send buffer to compressor");
            }
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
            temp_path, bins, ..
        } = self;

        let mut bins = Arc::into_inner(bins).expect("Could not get bins");

        let mut counts: DashMap<Vec<u8>, u32> = DashMap::new();

        let bins = bins
            .drain(..)
            .into_iter()
            .map(|x| (x.number, x.filename))
            .collect::<Vec<_>>();

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
                let kmers: Vec<Vec<u8>> = bincode::decode_from_slice(&kmers, bincode_config)
                    .expect("Could not decode buffer")
                    .0;

                for kmer in kmers {
                    if counts.contains_key(&kmer) {
                        counts.alter(&kmer, |_, count| count.saturating_add(1));
                    } else {
                        counts.insert(kmer, 1);
                    }
                }
            }
        });

        println!("Counts: {:?}", counts.len());
        // Convert dashmap into a hashmap
        let counts: HashMap<Vec<u8>, u32> = counts.into_iter().map(|(k, v)| (k, v)).collect();

        // Save to file
        let mut out_fh = BufWriter::new(
            File::create(format!("{}/counts.bin", temp_path))
                .expect("Could not create counts file"),
        );

        bincode::encode_into_std_write(
            counts,
            &mut out_fh,
            bincode::config::standard().with_fixed_int_encoding(),
        )
        .expect("Could not write to counts file");
    }
}

fn kmer_worker(
    reads_rx: crossbeam::channel::Receiver<Vec<(Bytes, Option<Bytes>)>>,
    kmer_tx: crossbeam::channel::Sender<Vec<(u64, Bytes)>>,
    kmer_rx: crossbeam::channel::Receiver<Vec<(u64, Bytes)>>,
    compression_rx: crossbeam::channel::Receiver<(usize, SuperKmerStorage)>,
    compression_tx: crossbeam::channel::Sender<(usize, SuperKmerStorage)>,
    output_rx: crossbeam::channel::Receiver<(usize, Vec<u8>)>,
    output_tx: crossbeam::channel::Sender<(usize, Vec<u8>)>,
    shutdown_flag: Arc<AtomicBool>,
    bins: Arc<Vec<KmerBin>>,
    bin_power: u8,
    k: u8,
    min_quality: u8,
    threads: usize,
) {
    let bin_mask = (1 << bin_power) - 1;
    let mut compressor = zstd::bulk::Compressor::new(-3).expect("Could not create compressor");

    let mut kmer_messages = Vec::new();
    let mut compression_messages = Vec::new();
    let mut output_messages = Vec::new();

    // did not seem to help
    // compressor.set_parameter(zstd::stream::raw::CParameter::Strategy(zstd::zstd_safe::zstd_sys::ZSTD_strategy::ZSTD_fast)).expect("Could not set compression level");

    // Create a thread-local buffer for each bin.
    let bin_count = bins.len();
    let mut local_buffers: Vec<SuperKmerStorage> = (0..bin_count)
        .map(|_| SuperKmerStorage::with_capacity(LOCAL_FLUSH_THRESHOLD, k as usize))
        .collect();

    let mut superkmers = Vec::with_capacity(128 * 1024);

    let mut did_work;

    loop {
        did_work = false;

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
            did_work = true;
        }

        if output_messages.is_empty() {
            if let Ok((bin, kmers)) = compression_rx.try_recv() {
                did_work = true;

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

                // output_tx
                // .send((bin, compressed))
                // .expect("Could not send compressed buffer to flusher");

                output_messages.push((bin, compressed));
            }
        }

        if compression_messages.is_empty() {
            if let Ok(kmers) = kmer_rx.try_recv() {
                did_work = true;

                let mut buffers_to_check = std::collections::HashSet::new();

                // For each kmer, calculate the bin and store it in the corresponding local buffer.
                for (minimizer, kmer) in kmers.into_iter() {
                    let hash = xxh3_64(&minimizer.to_ne_bytes());
                    let bin = hash & bin_mask;
                    let bin_index = bin as usize;

                    local_buffers[bin_index].add_superkmer(&kmer);
                    buffers_to_check.insert(bin_index);
                }

                // Check local buffers
                for bin_index in buffers_to_check.drain() {
                    let local_buf = &mut local_buffers[bin_index];
                    if local_buf.len() >= LOCAL_FLUSH_THRESHOLD {
                        let mut bin_lock = bins[bin_index].buffer.lock().unwrap();
                        bin_lock.append(local_buf);

                        if bin_lock.len() > GLOBAL_FLUSH_THRESHOLD {
                            let mut bin_buffer =
                                SuperKmerStorage::with_capacity(GLOBAL_FLUSH_THRESHOLD, k as usize);
                            std::mem::swap(&mut *bin_lock, &mut bin_buffer);
                            drop(bin_lock);

                            compression_messages.push((bin_index, bin_buffer));
                        }
                    }
                }
            }
        }

        if kmer_messages.is_empty() {
            if let Ok(reads) = reads_rx.try_recv() {
                did_work = true;

                for (mut seq, qual) in reads {
                    if let Some(qual) = qual {
                        // If average is less than min_quality, skip
                        let total_qual: u64 = qual.iter().map(|q| *q as u64).sum();
                        if (total_qual / qual.len() as u64) < min_quality as u64 {
                            continue;
                        }

                        // Otherwise mask sequence when quality is less than min_quality

                        let masked_seq: Bytes = seq
                            .iter()
                            .zip(qual.iter())
                            .map(|(base, q)| if q < &min_quality { b'N' } else { *base })
                            .collect();

                        seq = masked_seq;
                    } // If no quality, we don't worry about it

                    if seq.len() < k as usize {
                        continue;
                    }

                    let seq = seq.strip_returns();
                    let seq = seq.normalize(true);
                    let rc = seq.reverse_complement();

                    let seq = Bytes::from(seq.into_owned());
                    let rc = Bytes::from(rc);

                    let mut kmers = seq.canonical_kmers(k, &rc);

                    let mut superkmers_count = 0;
                    let mut total_kmers = 0;

                    // let mut kmer_min = NtHashIterator::new(&kmers.next().unwrap().1, 7)
                    // .expect("Could not create NtHashIterator")
                    // .min()
                    // .expect("Could not get min");

                    let mut superkmer_start_kmer_start = 0;

                    // let kmer_min = one_minimizer(AsciiSeq(&kmers.next().unwrap().1), 7);
                    let mut kmer_min = one_minimizer_filtered(&kmers.next().unwrap().1, 7);
                    // let mut kmer_min = seq.slice(kmer_min..kmer_min + 7);

                    for (i, kmer, rc) in kmers {
                        total_kmers += 1;
                        // let iter = NtHashIterator::new(&kmer, 7).expect("Could not create NtHashIterator");
                        // let min = iter.min().expect("Could not get min");
                        // let min = one_minimizer(AsciiSeq(&kmer), 7);
                        // let min = seq.slice(min..min + 7);
                        let min = one_minimizer_filtered(&kmer, 7);

                        if min != kmer_min {
                            // superkmers.push((superkmer_start_pos, pos, kmer_min, rc));

                            // Get the actual sequence / superkmer
                            // let superkmer = &seq[superkmer_start_kmer_start..i + k as usize];
                            let superkmer = seq.slice(superkmer_start_kmer_start..i + k as usize);
                            superkmers.push((kmer_min, superkmer));

                            // Insert and all that
                            kmer_min = min;
                            superkmer_start_kmer_start = i;
                            superkmers_count += 1;
                        }
                    }

                    // Last one
                    let superkmer = seq.slice(superkmer_start_kmer_start..);
                    superkmers.push((kmer_min, superkmer));

                    if superkmers.len() > 64 * 1024 {
                        kmer_messages.push(superkmers);
                        superkmers = Vec::with_capacity(64 * 1024);
                    }
                }
            }
        }

        if output_tx.len() < threads && !output_messages.is_empty() {
            // Prepare a temporary vector to hold messages that couldn’t be sent.
            let mut unsent_messages = Vec::with_capacity(output_messages.len());
            // Drain moves each element out without cloning.
            for message in output_messages.drain(..) {
                if let Err(crossbeam::channel::TrySendError::Full(message)) =
                    output_tx.try_send(message)
                {
                    // If the send fails (e.g. channel is full), push the message back.
                    unsent_messages.push(message);
                }
            }
            // Replace the original queue with the unsent messages.
            output_messages = unsent_messages;
        }

        if compression_tx.len() < threads && !compression_messages.is_empty() {
            // Prepare a temporary vector to hold messages that couldn’t be sent.
            let mut unsent_messages = Vec::with_capacity(compression_messages.len());
            // Drain moves each element out without cloning.
            for message in compression_messages.drain(..) {
                if let Err(crossbeam::channel::TrySendError::Full(message)) =
                    compression_tx.try_send(message)
                {
                    // If the send fails (e.g. channel is full), push the message back.
                    unsent_messages.push(message);
                }
            }
            // Replace the original queue with the unsent messages.
            compression_messages = unsent_messages;
        }

        if kmer_tx.len() < threads && !kmer_messages.is_empty() {
            // Prepare a temporary vector to hold messages that couldn’t be sent.
            let mut unsent_messages = Vec::with_capacity(kmer_messages.len());
            // Drain moves each element out without cloning.
            for message in kmer_messages.drain(..) {
                if let Err(crossbeam::channel::TrySendError::Full(message)) =
                    kmer_tx.try_send(message)
                {
                    // If the send fails (e.g. channel is full), push the message back.
                    unsent_messages.push(message);
                }
            }
            // Replace the original queue with the unsent messages.
            kmer_messages = unsent_messages;
        }

        if !did_work {
            println!("Worker sleeping - {} {} {} {} - {} {} {}", reads_rx.len(), kmer_tx.len(), compression_tx.len(), output_tx.len(), kmer_messages.len(), compression_messages.len(), output_messages.len());
            std::thread::sleep(std::time::Duration::from_millis(100));
        }
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
    buffer: Mutex<SuperKmerStorage>,
}

// kmers up to 31 bases long
pub fn count_kmers_file(kmer_counter: &mut KmerCounter, file: &str, k: u8, min_quality: u8) {
    let file = File::open(file).expect("Could not open file");
    let reader = BufReader::with_capacity(8 * 1024 * 1024, file);
    let mut reader = parse_fastx_reader(reader).expect("Invalid file");

    const FLUSH_BUFFER: usize = 16 * 1024;

    let mut reads = Vec::with_capacity(FLUSH_BUFFER);

    // debugging
    let mut processed_reads = 0;

    // let mut kmers_to_submit = Vec::with_capacity(8 * 1024);

    while let Some(record) = reader.next() {
        processed_reads += 1;
        if processed_reads % 100000 == 0 {
            println!("Processed {} reads", processed_reads);
        }
        let record = record.expect("Error reading record");
        let seq = record.seq();
        let qual = record.qual();

        let qual = if let Some(qual) = qual {
            Some(Bytes::from(qual.to_vec()))
        } else {
            None
        };

        let seq = Bytes::from(seq.into_owned());

        reads.push((seq, qual));

        if reads.len() >= FLUSH_BUFFER {
            kmer_counter.submit(reads);
            reads = Vec::with_capacity(FLUSH_BUFFER);
        }
    }

    if !reads.is_empty() {
        kmer_counter.submit(reads);
    }
}

#[inline(always)]
pub fn one_minimizer_filtered(seq: &[u8], m: usize) -> u64 {
    assert!(seq.len() >= m, "Sequence length must be at least m");

    // Precompute 256^(m-1) to remove the contribution of the dropped byte.
    let multiplier = 256u64.pow((m - 1) as u32);

    // Compute the candidate for the first window.
    let mut candidate = seq.iter().take(m).fold(0u64, |acc, &b| (acc << 8) | b as u64);
    let mut best_overall = candidate;
    let mut best_filtered: Option<u64> = if passes_rules(candidate, m) {
        Some(candidate)
    } else {
        None
    };

    // Slide the window over the sequence, updating the candidate in O(1) time.
    for i in 1..=seq.len() - m {
        candidate = (candidate - (seq[i - 1] as u64) * multiplier) * 256 + seq[i + m - 1] as u64;
        best_overall = best_overall.min(candidate);
        if passes_rules(candidate, m) {
            best_filtered = Some(match best_filtered {
                Some(current) => current.min(candidate),
                None => candidate,
            });
        }
    }

    // If a candidate passed the rules, return it. Otherwise, return the overall minimizer.
    best_filtered.unwrap_or(best_overall)
}

/// Checks whether the candidate m‑mer passes filtering rules:
/// - It must not start with "AA"
/// - It must not be a homopolymer (all bases identical)
#[inline(always)]
fn passes_rules(candidate: u64, m: usize) -> bool {
    // Extract the first two bytes (most-significant bytes, as the candidate is big‑endian).
    let first_two = candidate >> ((m - 2) * 8);
    if first_two == ((b'A' as u64) << 8 | (b'A' as u64)) {
        return false;
    }

    // Check for a homopolymer: ensure at least one base differs.
    let first_byte = (candidate >> ((m - 1) * 8)) & 0xFF;
    for i in 1..m {
        let byte = (candidate >> ((m - 1 - i) * 8)) & 0xFF;
        if byte != first_byte {
            return true; // At least one base is different.
        }
    }
    false // All bases are the same → too simple.
}
