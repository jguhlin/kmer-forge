use needletail::{parse_fastx_file, Sequence};
use std::{borrow::Cow, collections::HashMap};

// Returns
// 0: Number of reads
// 1: Hashmap Kmer, Count

// kmers up to 255
pub fn parse_file(file: &str, k: u8, min_quality: u8) -> (u64, HashMap<Vec<u8>, u64>) {
    let mut count = 0;
    let mut kmer_count_old = HashMap::new();
    let mut reader = parse_fastx_file(file).expect("Invalid file");

    while let Some(record) = reader.next() {
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

        // Kmers for both seq and rc
        for kmer in seq.kmers(k) {
            *kmer_count.entry(kmer.to_vec()).or_insert(0) += 1;
        }

        for kmer in rc.kmers(k) {
            *kmer_count.entry(kmer.to_vec()).or_insert(0) += 1;
        }
    }

    (count, kmer_count_old)
    

}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kmer_counts() {
        let filename = "E6_22T2Y5LT3_TCTAGGCGCG-CGAAGGTTAA_L001_R1.fastq.gz";
        let (count, kmer_count) = parse_file(filename, 21, 30);

        println!("Total reads: {}", count);
        println!("Unique kmers: {}", kmer_count.len());
        panic!();
        
    }
}
