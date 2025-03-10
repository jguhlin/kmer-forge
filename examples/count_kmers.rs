use kmer_forge::*;

fn main() {
    let filename = "E6_22T2Y5LT3_TCTAGGCGCG-CGAAGGTTAA_L001_R1.fastq.gz";
    let (count, kmer_count) = parse_file(filename, 21, 30);

    println!("Total reads: {}", count);
    println!("Unique kmers: {}", kmer_count.len());
    panic!();
}