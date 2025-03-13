use kmer_forge::*;

fn main() {
    let filename = "E6_22T2Y5LT3_TCTAGGCGCG-CGAAGGTTAA_L001_R1.fastq.gz";
    let k = 21;

    let mut kmer_counter = KmerCounter::new(k, "temp".to_string(), 32, 512 * 1024, 8);

    count_kmers_file(&mut kmer_counter, filename, k, 20);

    kmer_counter.stop_gathering();
    kmer_counter.merge_bins();



}
