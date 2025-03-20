use bincode::{Decode, Encode};

#[derive(Encode, Decode)]
pub struct SuperKmerStorage {
    pub bytes: Vec<u8>, // Todo, convert to Bytes from tokio?
    pub lengths: Vec<u8>, // TODO: Make changeable
}

impl SuperKmerStorage {
    pub fn new() -> Self {
        Self {
            bytes: Vec::new(),
            lengths: Vec::new(),
        }
    }

    pub fn len(&self) -> usize {
        self.lengths.len()
    }
    
    pub fn is_empty(&self) -> bool {
        self.lengths.is_empty()
    }

    pub fn add_superkmer(&mut self, superkmer: &[u8]) {
        self.bytes.extend_from_slice(superkmer);
        self.lengths.push(superkmer.len() as u8);
    }

    pub fn append(&mut self, other: &mut Self) {
        self.bytes.append(&mut other.bytes);
        self.lengths.append(&mut other.lengths);
    }

    pub fn clear(&mut self) {
        self.bytes.clear();
        self.lengths.clear();
    }
}

impl Iterator for SuperKmerStorage {
    type Item = Vec<u8>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.bytes.is_empty() {
            return None;
        }

        let length = self.lengths.remove(0) as usize;
        let superkmer = self.bytes.drain(..length).collect();

        Some(superkmer)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let len = self.lengths.len();
        (len, Some(len))
    }

    fn count(self) -> usize {
        self.lengths.len()
    }
}