use crate::util::{ACGT_TO_0123, DNA_ALPHABET};
use crate::{util, SbwtIndex, SubsetSeq};
use byteorder::{LittleEndian, ReadBytesExt};
use std::fmt::Debug;
use std::io::{Read, Write};

/// A table storing the SBWT intervals of all 4^p possible p-mers.
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct PrefixLookupTable {
    /// ranges\[i\] is the interval of the p-mer with colexicographic rank
    /// i in the sorted list of all possible p-mers.
    /// If the p-mer does not exist in the SBWT, the range is [0..0).
    pub ranges: Vec<std::ops::Range<usize>>,

    /// Prefix length p.
    pub prefix_length: usize,
}
impl PrefixLookupTable {
    // Clippy false positive. It's actually intended like this
    #[allow(clippy::single_range_in_vec_init)]
     pub fn new_empty(n_sets_in_sbwt: usize) -> PrefixLookupTable {
        Self{ranges: vec![0..n_sets_in_sbwt], prefix_length: 0}
    }


    pub fn new<SS: SubsetSeq>(sbwt: &SbwtIndex<SS>, prefix_length: usize) -> PrefixLookupTable {
        let mut pmer = vec![0u8; prefix_length];
        let mut ranges = vec![0..0; num::pow(4_usize, prefix_length)];
        for x in 0..num::pow(4, prefix_length) as u64{
            pmer.clear();

            // Construct the p-mer string
            for i in 0..prefix_length {
                let char_idx = (x >> (2*(prefix_length - 1 - i))) & 0x3;
                let c = DNA_ALPHABET[char_idx as usize];
                pmer.push(c);
            }

            if let Some(range) = sbwt.search(&pmer) {
                log::trace!("range [{}, {})", range.start, range.end);
                ranges[x as usize] = range;
             
            } // Else left as 0..0
        }
        PrefixLookupTable {ranges, prefix_length}
    }
    pub fn lookup(&self, prefix: &[u8]) -> std::ops::Range<usize> {
        assert!(prefix.len() == self.prefix_length);
        let mut table_idx = 0_usize;
        for (i, c) in prefix.iter().rev().enumerate() {
            let char_idx = ACGT_TO_0123[*c as usize];
            if char_idx == 255 {
                return 0..0; // Not a DNA character
            }
            table_idx |= ((char_idx as u64) << (2*i)) as usize;
        }
        self.ranges[table_idx].clone()
    }

    /// Loads a a prefix lookup table that was previosly serialized with
    /// [PrefixLookupTable::serialize].
    pub fn load<R: Read>(input: &mut R) -> std::io::Result<Self> {

        let prefix_length = byteorder::ReadBytesExt::read_u64::<LittleEndian>(input).unwrap()
            as usize;
        let ranges_len = byteorder::ReadBytesExt::read_u64::<LittleEndian>(input).unwrap()
            as usize;
        let mut ranges = vec![0..0; ranges_len];
        for range in ranges.iter_mut(){
            let start = byteorder::ReadBytesExt::read_u64::<LittleEndian>(input).unwrap()
                as usize;
            let end = byteorder::ReadBytesExt::read_u64::<LittleEndian>(input).unwrap()
                as usize;
            *range = start..end;
        }

        Ok(Self{ranges, prefix_length})
    }

    /// Write the lookup table to the given writer.
    /// The lookup table can be then later loaded with [PrefixLookupTable::load].
    /// Returns number of bytes written.
    pub fn serialize<W: Write>(&self, out: &mut W) -> std::io::Result<usize> {
        let mut n_written = 0_usize;
        n_written += util::write_bytes(out, &(self.prefix_length as u64).to_le_bytes())?;
        n_written += util::write_bytes(out, &(self.ranges.len() as u64).to_le_bytes())?;
        for range in self.ranges.iter(){
            n_written += util::write_bytes(out, &(range.start as u64).to_le_bytes())?;
            n_written += util::write_bytes(out, &(range.end as u64).to_le_bytes())?;
        }
        Ok(n_written)
    }
}

