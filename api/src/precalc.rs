use crate::util::{ACGT_TO_0123, DNA_ALPHABET};
use crate::{util, SbwtIndex, SubsetSeq};
use byteorder::{LittleEndian, ReadBytesExt};
use std::fmt::Debug;
use std::io::{Read, Write};

pub enum SerializeMode {
    Raw,
    VByte,
}

const MODE: SerializeMode = SerializeMode::VByte;

/// A table storing the SBWT intervals of all 4^p possible p-mers.
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct PrefixLookupTable {
    /// ranges[i] is the interval of the p-mer with colexicographic rank
    /// i in the sorted list of all possible p-mers.
    /// If the p-mer does not exist in the SBWT, the range is [0..0).
    pub ranges: Vec<std::ops::Range<usize>>,

    /// Prefix length p.
    pub prefix_length: usize,

    /// Size of the serialized table in bytes when loaded/serialized.
    pub size_bytes: usize,
}

impl PrefixLookupTable {
    // Clippy false positive. It's actually intended like this
    #[allow(clippy::single_range_in_vec_init)]
    pub fn new_empty(n_sets_in_sbwt: usize) -> PrefixLookupTable {
        Self {
            ranges: vec![0..n_sets_in_sbwt],
            prefix_length: 0,
            size_bytes: 0,
        }
    }

    pub fn new<SS: SubsetSeq>(sbwt: &SbwtIndex<SS>, prefix_length: usize) -> PrefixLookupTable {
        let mut pmer = vec![0u8; prefix_length];
        let mut ranges = vec![0..0; num::pow(4_usize, prefix_length)];
        for x in 0..num::pow(4, prefix_length) as u64 {
            pmer.clear();

            // Construct the p-mer string
            for i in 0..prefix_length {
                let char_idx = (x >> (2 * (prefix_length - 1 - i))) & 0x3;
                let c = DNA_ALPHABET[char_idx as usize];
                pmer.push(c);
            }

            if let Some(range) = sbwt.search(&pmer) {
                ranges[x as usize] = range;
            } // Else left as 0..0
        }
        PrefixLookupTable {
            ranges,
            prefix_length,
            size_bytes: 0,
        }
    }

    pub fn lookup(&self, prefix: &[u8]) -> std::ops::Range<usize> {
        assert!(prefix.len() == self.prefix_length);
        let mut table_idx = 0_usize;
        for (i, c) in prefix.iter().rev().enumerate() {
            let char_idx = ACGT_TO_0123[*c as usize];
            if char_idx == 255 {
                return 0..0; // Not a DNA character
            }
            table_idx |= ((char_idx as u64) << (2 * i)) as usize;
        }
        self.ranges[table_idx].clone()
    }

    /// Loads a prefix lookup table that was previously serialized with
    /// [PrefixLookupTable::serialize].
    pub fn load<R: Read>(input: &mut R) -> std::io::Result<Self> {
        let mut bytes_read = 0_usize;

        let (prefix_length, ranges_len) = match MODE {
            SerializeMode::Raw => {
                let prefix_length = input.read_u64::<LittleEndian>()? as usize;
                let ranges_len = input.read_u64::<LittleEndian>()? as usize;
                bytes_read += 16; // 8 bytes for prefix_length + 8 bytes for ranges_len
                (prefix_length, ranges_len)
            }
            SerializeMode::VByte => {
                let (prefix_length, n1) = Self::vbyte_read(input)?;
                let (ranges_len, n2) = Self::vbyte_read(input)?;
                bytes_read += n1 + n2;
                (prefix_length as usize, ranges_len as usize)
            }
        };

        let mut ranges = vec![0..0; ranges_len];

        match MODE {
            SerializeMode::Raw => {
                for range in ranges.iter_mut() {
                    let start = input.read_u64::<LittleEndian>()? as usize;
                    let end = input.read_u64::<LittleEndian>()? as usize;
                    *range = start..end;
                }
                bytes_read += ranges_len * 16; // 16 bytes per range (2 x u64)
            }
            SerializeMode::VByte => {
                for range in ranges.iter_mut() {
                    let (start, n1) = Self::vbyte_read(input)?;
                    let (len, n2) = Self::vbyte_read(input)?;
                    bytes_read += n1 + n2;
                    let end = start + len; // Reconstruct end using stored delta
                    *range = (start as usize)..(end as usize);
                }
            }
        }

        Ok(Self {
            ranges,
            prefix_length,
            size_bytes: bytes_read,
        })
    }

    pub fn vbyte_write<W: Write>(out: &mut W, mut val: u64) -> std::io::Result<usize> {
        let mut n_written = 0;
        loop {
            let mut byte = (val & 0x7F) as u8;
            val >>= 7;

            if val == 0 {
                byte |= 0x80;
            }

            n_written += util::write_bytes(out, &[byte])?;

            if val == 0 {
                break;
            }
        }
        Ok(n_written)
    }

    pub fn vbyte_read<R: Read>(input: &mut R) -> std::io::Result<(u64, usize)> {
        let mut val: u64 = 0;
        let mut shift = 0;
        let mut n_read = 0;
        let mut buf = [0u8; 1];

        loop {
            if shift >= 64 {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "VByte integer overflow for u64: sequence exceeds 64 bits",
                ));
            }

            input.read_exact(&mut buf)?;
            n_read += 1;
            let byte = buf[0];

            if shift == 63 && (byte & 0x7E) != 0 {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "VByte integer overflow for u64: value exceeds u64::MAX",
                ));
            }

            let chunk = (byte & 0x7F) as u64;
            val |= chunk << shift;

            if (byte & 0x80) != 0 {
                break;
            }

            shift += 7;
        }

        Ok((val, n_read))
    }

    /// Write the lookup table to the given writer.
    /// Returns number of bytes written.
    pub fn serialize<W: Write>(&self, out: &mut W) -> std::io::Result<usize> {
        let mut n_written = 0_usize;

        match MODE {
            SerializeMode::Raw => {
                n_written += util::write_bytes(out, &(self.prefix_length as u64).to_le_bytes())?;
                n_written += util::write_bytes(out, &(self.ranges.len() as u64).to_le_bytes())?;

                for range in self.ranges.iter() {
                    n_written += util::write_bytes(out, &(range.start as u64).to_le_bytes())?;
                    n_written += util::write_bytes(out, &(range.end as u64).to_le_bytes())?;
                }
            }
            SerializeMode::VByte => {
                n_written += Self::vbyte_write(out, self.prefix_length as u64)?;
                n_written += Self::vbyte_write(out, self.ranges.len() as u64)?;

                for range in self.ranges.iter() {
                    let start = range.start as u64;
                    let len = (range.end - range.start) as u64; // Delta compression

                    n_written += Self::vbyte_write(out, start)?;
                    n_written += Self::vbyte_write(out, len)?;
                }
            }
        }

        Ok(n_written)
    }
}