use bumpalo::Bump;
use bumpalo::collections::Vec as BumpVec;
use crossbeam::channel::{Receiver, Sender, bounded, unbounded};
// use growable_bloom_filter::GrowableBloom;
use needletail::{Sequence, parse_fastx_file, parse_fastx_reader};
use pulp::Arch;
use xxhash_rust::xxh3::xxh3_64;
use lz4_flex::block::*;

use std::borrow::Cow;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::sync::atomic::{AtomicBool, AtomicU32, Ordering};
use std::sync::{Arc, Mutex};
use std::thread::{self, JoinHandle};

mod utils;

use utils::*;

pub struct KmerCounter {
    k: u8,
    temp_path: String,
    bins: Arc<Vec<KmerBin>>,
    bin_count: u16,
    threads: usize,
    workers: Vec<JoinHandle<()>>,
    kmer_tx: Sender<Vec<u64>>,
    kmer_rx: Receiver<Vec<u64>>,
    compression_tx: Sender<(usize, Vec<u64>)>,
    compression_rx: Receiver<(usize, Vec<u64>)>,
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
            let out_fh = BufWriter::new(File::create(bin_path).expect("Could not create bin file"));
            let out_fh = Mutex::new(out_fh);
            bins.push(KmerBin {
                number: i as u16,
                out_fh,
                buffer: Mutex::new(Vec::with_capacity(buffer_flush_size)),
                buffer_flush_size,
            });
        }

        let bins = Arc::new(bins);

        // Create the channels
        let (kmer_tx, kmer_rx): (Sender<Vec<u64>>, Receiver<Vec<u64>>) = bounded(64);
        let (compression_tx, compression_rx): (Sender<(usize, Vec<u64>)>, Receiver<(usize, Vec<u64>)>) = bounded(64);
        let (output_tx, output_rx): (Sender<(usize, Vec<u8>)>, Receiver<(usize, Vec<u8>)>) = unbounded();

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
                kmer_worker(kmer_rx, compression_rx, compression_tx, output_rx, output_tx, shutdown_flag, bins, bin_power);
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

    pub fn submit(&self, kmers: Vec<u64>) {
        self.kmer_tx
            .send(kmers)
            .expect("Could not send kmers to worker");
    }

    pub fn try_submit(&self, kmers: Vec<u64>) -> Result<(), crossbeam::channel::TrySendError<Vec<u64>>> {
        self.kmer_tx.try_send(kmers)
    }

    pub fn shutdown(&mut self) {

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
                self.compression_tx.send((bin.number as usize, bin_lock.clone())).expect("Could not send buffer to compressor");
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
}

fn kmer_worker(
    kmer_rx: crossbeam::channel::Receiver<Vec<u64>>,
    compression_rx: crossbeam::channel::Receiver<(usize, Vec<u64>)>,
    compression_tx: crossbeam::channel::Sender<(usize, Vec<u64>)>,
    output_rx: crossbeam::channel::Receiver<(usize, Vec<u8>)>,
    output_tx: crossbeam::channel::Sender<(usize, Vec<u8>)>,
    shutdown_flag: Arc<AtomicBool>,
    bins: Arc<Vec<KmerBin>>,
    bin_power: u8,
) {

    let bin_mask = (1 << bin_power) - 1;
    const FLUSH_THRESHOLD: usize = 8192;

    let mut compressor = zstd::bulk::Compressor::new(-1).expect("Could not create compressor");

    // did not seem to help
    // compressor.set_parameter(zstd::stream::raw::CParameter::Strategy(zstd::zstd_safe::zstd_sys::ZSTD_strategy::ZSTD_fast)).expect("Could not set compression level");
    
    // Create a thread-local buffer for each bin.
    let bin_count = bins.len();
    let mut local_buffers: Vec<Vec<u64>> = (0..bin_count)
        .map(|_| Vec::with_capacity(FLUSH_THRESHOLD))
        .collect();

    let mut bins_to_submit = Vec::with_capacity(128);

    loop {
        if shutdown_flag.load(Ordering::Relaxed) {
            break;
        }

        match compression_rx.try_recv() {
            Ok((bin, mut kmers)) => {
                kmers.sort_unstable();

                let encoded = bincode::encode_to_vec(&kmers, bincode::config::standard().with_fixed_int_encoding()).expect("Could not write to bin file");

                // zstd
                // let compressed = compressor.compress(&encoded).expect("Could not compress buffer");

                // lz4_flex
                // let compressed = compress(&encoded);

                // no compression
                let compressed = encoded;

                output_tx.send((bin, compressed)).expect("Could not send compressed buffer to flusher");
            },
            Err(crossbeam::channel::TryRecvError::Disconnected) => break, // Or panic? Something has gone wrong
            Err(crossbeam::channel::TryRecvError::Empty) => {}, // Do nothing, continue on
        }

        match output_rx.try_recv() {
            Ok((bin, compressed)) => {
                let mut bin_lock = bins[bin].out_fh.lock().unwrap();
                bincode::encode_into_std_write(
                    compressed,
                    &mut *bin_lock,
                    bincode::config::standard().with_fixed_int_encoding(),
                )
                .expect("Could not write to bin file");
                drop(bin_lock);
            },
            Err(crossbeam::channel::TryRecvError::Disconnected) => break, // Or panic? Something has gone wrong
            Err(crossbeam::channel::TryRecvError::Empty) => {}, // Do nothing, continue on
        }

        let kmers = match kmer_rx.try_recv() {
            Err(crossbeam::channel::TryRecvError::Empty) => continue,
            Err(crossbeam::channel::TryRecvError::Disconnected) => break,
            Ok(kmers) => kmers,
        };

        // For each kmer, calculate the bin and store it in the corresponding local buffer.
        for kmer in kmers {
            let hash = xxh3_64(&kmer.to_ne_bytes());
            let bin = hash & bin_mask;
            let bin_index = bin as usize;

            local_buffers[bin_index].push(kmer);

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

        if !bins_to_submit.is_empty() && compression_rx.len() as f32 <= 0.8 * compression_rx.capacity().unwrap() as f32 {
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

                compression_tx.send((bin, bin_buffer)).expect("Could not send buffer to compressor");
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
    buffer: Mutex<Vec<u64>>,
    buffer_flush_size: usize,
}

// Returns
// 0: Number of reads
// 1: Hashmap Kmer, Count

// kmers up to 255
pub fn parse_file(file: &str, k: u8, min_quality: u8) {
    // -> (u64, KmerCounter) {
    let mut count = 0;
    let file = File::open(file).expect("Could not open file");
    let reader = BufReader::with_capacity(8 * 1024 * 1024, file);
    let mut reader = parse_fastx_reader(reader).expect("Invalid file");
    // let mut reader = parse_fastx_file(file).expect("Invalid file");
    let mut kmer_counter = KmerCounter::new(k, "temp".to_string(), 32, 512 * 1024, 8);

    println!("Kmer counter created");

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

        // We are here. So process
        count += 1;

        let seq = seq.strip_returns();
        let seq = seq.normalize(true);

        let rc = seq.reverse_complement();

        let mut kmers = seq.kmers(k);
        let mut rolling_encoder =
            RollingKmer3::new(&kmers.next().unwrap(), k as usize).expect("Sequence too short");
        // Add the first code to the vec
        kmers_to_submit.push(rolling_encoder.code);
        // Roll the rest
        for kmer in kmers {
            rolling_encoder.roll(kmer[k as usize - 1]);
            kmers_to_submit.push(rolling_encoder.code);
        }

        let mut rc_kmers = rc.kmers(k);
        let mut rolling_encoder =
            RollingKmer3::new(&rc_kmers.next().unwrap(), k as usize).expect("Sequence too short");
        // Add the first code to the vec
        kmers_to_submit.push(rolling_encoder.code);
        // Roll the rest
        for kmer in rc_kmers {
            rolling_encoder.roll(kmer[k as usize - 1]);
            kmers_to_submit.push(rolling_encoder.code);
        }

        if kmers_to_submit.len() >= 8 * 1024 {
            match kmer_counter.try_submit(kmers_to_submit) {
                Ok(_) => {
                    kmers_to_submit = Vec::with_capacity(8 * 1024);
                }
                Err(crossbeam::channel::TrySendError::Full(kmers)) => {
                    kmers_to_submit = kmers;
                }, // Process more reads then...
                Err(crossbeam::channel::TrySendError::Disconnected(_kmers)) => {
                    panic!("Kmer counter disconnected before shutdown");
                },
            }
            
        }
    }

    kmer_counter.shutdown();

    // (count, kmer_counter)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rolling_kmer() {
        let k: usize = 21;
        // Example sequence: must be at least k nucleotides.
        let seq = b"AACGTNACGTNACGTNACGTN"; // exactly 21 nucleotides
        let mut rk = RollingKmer3::new(seq, k).expect("Sequence too short");
        println!("Initial encoding: {:#018x}", rk.code);

        // Roll the k-mer: remove the leftmost nucleotide ('A') and append 'T'.
        rk.roll(b'T');
        println!("After rolling:   {:#018x}", rk.code);

        // Roll the k-mer: remove the leftmost nucleotide ('C') and append 'G'.
        rk.roll(b'G');

        // Start a new one to confirm that the previous one was rolled correctly.
        let mut rk2 = RollingKmer3::new(b"CGTNACGTNACGTNACGTNTG", k).expect("Sequence too short");
        println!("Initial encoding: {:#018x}", rk2.code);

        assert!(rk.code == rk2.code);

        // Now for k=9
        let k: usize = 9;
        let seq = b"ACGTNACGT";
        let mut rk = RollingKmer3::new(seq, k).expect("Sequence too short");
        rk.roll(b'G');
        rk.roll(b'C');

        let mut rk2 = RollingKmer3::new(b"GTNACGTGC", k).expect("Sequence too short");
        assert!(rk.code == rk2.code);
    }
}
