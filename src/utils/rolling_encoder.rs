#[inline(always)]
fn nucleotide_to_3bit(nuc: u8) -> u64 {
    match nuc {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'T' | b't' => 2,
        b'G' | b'g' => 3,
        b'N' | b'n' => 4,
        _ => 4, // default unknown bases to N
    }
}

/// A rolling k-mer encoder using 3 bits per nucleotide.
/// With 3 bits per base, a u64 can hold up to 21 nucleotides (63 bits).
#[derive(Debug, Clone, Copy)]
pub struct RollingKmer3 {
    pub k: usize,  // k-mer length (must be <= 21)
    pub code: u64, // current encoded value
    mask: u64,     // precomputed mask for rolling (lower 3*(k-1) bits)
}

impl RollingKmer3 {
    /// Creates a new RollingKmer3 by encoding the first `k` nucleotides of `seq`.
    /// Returns None if the sequence is too short, k is 0, or k > 21.
    #[inline(always)]
    pub fn new(seq: &[u8], k: usize) -> Option<Self> {
        if seq.len() < k || k == 0 || k > 21 {
            return None;
        }
        let mut code: u64 = 0;
        // Use unsafe indexing since we've already checked the length.
        unsafe {
            for i in 0..k {
                code = (code << 3) | nucleotide_to_3bit(*seq.get_unchecked(i));
            }
        }
        let bits_per_base = 3;
        // Precompute the mask to keep the lower 3*(k-1) bits.
        let mask = (1u64 << (bits_per_base * (k - 1))) - 1;
        Some(Self { k, code, mask })
    }

    /// Rolls the k-mer by dropping the leftmost nucleotide and appending the new nucleotide.
    /// This is done by masking out the top 3 bits (for a k-mer of length k, keeping 3*(k-1) bits),
    /// shifting left by 3, then OR-ing in the new nucleotide's 3-bit code.
    #[inline(always)]
    pub fn roll(&mut self, new: u8) {
        let bits_per_base = 3;
        self.code = ((self.code & self.mask) << bits_per_base) | nucleotide_to_3bit(new);
    }
}
