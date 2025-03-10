use bumpalo::Bump;
use bumpalo::collections::Vec as BumpVec;
use needletail::{Sequence, parse_fastx_file};
use std::{borrow::Cow, collections::HashMap};

enum KmerStorage {
    Singleton,
    LowFreq,
    HighFreq,
}

pub struct KmerCounter {
    kmers: Vec<(u64, Vec<u8>, usize, KmerStorage)>, // hash, kmer, index, storage loc
    buffer: Vec<Vec<u8>>,
    // Singletons don't need to be stored!
    low_freq: Vec<(Vec<u8>, u8)>,

    high_freq: Vec<(Vec<u8>, u64)>,

    // clean-up count
    clean_up: u64,

    // todo: reuse low freq indices when they move to high freq?
    // But after 640k reads only 4.6k high freq kmers
    // so probably not worth it
}

impl KmerCounter {
    pub fn new() -> Self {
        KmerCounter {
            kmers: Vec::with_capacity(64 * 1024 * 1024),
            buffer: Vec::with_capacity(1024 * 1024),
            low_freq: Vec::new(),
            high_freq: Vec::new(),
            clean_up: 32 * 1024 * 1024, // By default, clean up every 1024 reads
        }
    }

    pub fn len(&self) -> usize {
        self.kmers.len()
    }

    pub fn buffer_len(&self) -> usize {
        self.buffer.len()
    }

    pub fn record_from_iter(&mut self, kmers: impl Iterator<Item = Vec<u8>>) {
        self.buffer.extend(kmers);

        if self.buffer.len() > self.clean_up as usize {
            self.process_buffer();
        }
    }

    pub fn record_kmer(&mut self, kmer: &[u8]) {
        let mut kmer_b = Vec::default();
        kmer_b.extend_from_slice(kmer);
        self.buffer.push(kmer_b);

        if self.buffer.len() > self.clean_up as usize {
            self.process_buffer();
        }
    }

    pub fn process_buffer(&mut self) {

        println!("Processing Buffer");
        self.buffer.sort_unstable();

        println!("Counting");
        // Count kmers in buffer
        // Buffer is sorted, so we can count them in one pass
        let mut counts: HashMap<Vec<u8>, u64> = HashMap::new();
        let mut last_kmer: Option<&[u8]> = None;
        for kmer in &self.buffer {
            if Some(kmer.as_slice()) == last_kmer {
                *counts.get_mut(kmer).unwrap() += 1;
            } else {
                counts.insert(kmer.clone(), 1);
            }
            last_kmer = Some(kmer.as_slice());
        }
        
        println!("Counted {} kmers", counts.len());

        // todo generate entries in a separate vec, then add and sort all at once later

        let mut kmer_additions: Vec<(u64, Vec<u8>, usize, KmerStorage)> = Vec::new();

        // Process counts
        for (kmer, count) in counts {
            let x_hash = xxhash_rust::xxh3::xxh3_64(&kmer);

            // Find kmer in kmer list
            match self
                .kmers
                .binary_search_by(|(hash, _, _, _)| hash.cmp(&x_hash))
            {
                Ok(i) => {
                    // Kmer already exists
                    match self.kmers[i].3 {
                        KmerStorage::Singleton => {
                            // If singleton, move to low freq
                            if count <= u8::MAX as u64 {
                                self.kmers[i].3 = KmerStorage::LowFreq;
                                let index = self.low_freq.len();
                                self.low_freq
                                    .push((self.kmers[i].1.clone(), 1 + count as u8));
                                self.kmers[i].2 = index;
                            } else {
                                let index = self.high_freq.len();
                                self.kmers[i].3 = KmerStorage::HighFreq;
                                self.kmers[i].2 = index;
                                self.high_freq.push((self.kmers[i].1.clone(), 1 + count));
                            }
                        }
                        KmerStorage::LowFreq => {
                            // If low freq, increment count

                            // Confirm it will remain low freq
                            // todo: now it's duplicated in both low and high freq
                            let index = self.kmers[i].2;
                            if self.low_freq[index].1 as u64 + count > u8::MAX as u64 {
                                let index = self.high_freq.len();
                                self.kmers[i].2 = index;
                                self.kmers[i].3 = KmerStorage::HighFreq;
                                self.high_freq.push((
                                    self.kmers[i].1.clone(),
                                    self.low_freq[index].1 as u64 + count,
                                ));
                                // Do not remove as it changes the index
                            } else {
                                self.low_freq[index].1 += count as u8;
                            }
                        }
                        KmerStorage::HighFreq => {
                            let index = self.kmers[i].2;
                            // If high freq, increment count
                            self.high_freq[index].1 += count;
                        }
                    }
                }
                Err(i) => {
                    // Kmer does not exist
                    if count == 1 {
                        // If singleton, don't store
                        kmer_additions.push((x_hash, kmer.to_vec(), 0, KmerStorage::Singleton));
                    } else {
                        // If not singleton, store
                        if count <= u8::MAX as u64 {
                            let index = self.low_freq.len();
                            self.low_freq.push((kmer.to_vec(), 1 + count as u8));
                            kmer_additions.push((x_hash, kmer.to_vec(), index, KmerStorage::LowFreq));
                            
                        } else {
                            let index = self.high_freq.len();
                            kmer_additions.push((x_hash, kmer.to_vec(), index, KmerStorage::HighFreq));
                            self.high_freq.push((kmer.to_vec(), 1 + count));
                        }
                    }
                }
            }
        }

        println!("Adding {} kmers", kmer_additions.len());
        self.kmers.extend(kmer_additions);


        println!("Sorting kmers");
        // Sort kmers
        self.kmers.sort_unstable_by(|(hash1, _, _, _), (hash2, _, _, _)| hash1.cmp(hash2));

        println!("Sorted kmers");

        // Clear buffer
        self.buffer.clear();

        println!("Processed {} kmers", self.kmers.len());
        println!("Singletons: {}", self.kmers.len() - self.low_freq.len() - self.high_freq.len());
        println!("Low Freq: {}", self.low_freq.len());
        println!("High Freq: {}", self.high_freq.len());
        
    }
}

// Returns
// 0: Number of reads
// 1: Hashmap Kmer, Count

// kmers up to 255
pub fn parse_file(file: &str, k: u8, min_quality: u8) -> (u64, KmerCounter) {
    let mut count = 0;
    let mut kmer_counter = KmerCounter::new();
    let mut reader = parse_fastx_file(file).expect("Invalid file");

    // debugging
    let mut processed_reads = 0;

    while let Some(record) = reader.next() {
        processed_reads += 1;
        if processed_reads % 10000 == 0 {
            println!("Processed {} reads", processed_reads);
        }
        let record = record.expect("Error reading record");
        let mut seq = record.seq();
        let qual = record.qual();

        if let Some(qual) = qual {
            // If average is less than min_quality, skip
            let total_qual: u64 = qual.iter().map(|q| *q as u64).sum();
            if (total_qual / qual.len() as u64) < min_quality as u64 {
                continue;
            }

            // Otherwise mask sequence when quality is less than min_quality
            let masked_seq: Vec<u8> = seq
                .iter()
                .zip(qual.iter())
                .map(|(base, q)| if q < &min_quality { b'N' } else { *base })
                .collect();

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

        kmer_counter.record_from_iter(seq.kmers(k).map(|kmer| kmer.to_vec()));
        kmer_counter.record_from_iter(rc.kmers(k).map(|kmer| kmer.to_vec()));
    }

    (count, kmer_counter)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kmer_counts() {
        
    }
}
