//! Prefix lookup tables: precomputed SBWT intervals of all 4^p possible p-mers, used to
//! skip the first p steps of a search.
//!
//! The table implementation is chosen at run time with [PrefixTableKind]. An index stores a
//! [PrefixTableVariant], which serializes a tag identifying its kind, so a loaded index knows
//! which kind of table it was built with.

use crate::util::{ACGT_TO_0123, DNA_ALPHABET};
use crate::{SbwtIndex, SubsetSeq};
use byteorder::{LittleEndian, ReadBytesExt};
use std::io::{Read, Write};
use std::ops::Range;

/// Identifies a prefix lookup table implementation, without holding a table.
/// Use this to choose which kind of table to build, e.g. from a command line flag.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum PrefixTableKind {
    /// Ranges stored as a flat vector, serialized as raw 64-bit integers. See [VectorLookupTable].
    #[default]
    Vector,
    /// Ranges stored as a flat vector, serialized with VByte + delta encoding. See [VByteLookupTable].
    VByte,
}

impl PrefixTableKind {
    /// Every kind.
    pub const ALL: [PrefixTableKind; 2] = [Self::Vector, Self::VByte];

    /// The names accepted by [str::parse], in the same order as [PrefixTableKind::ALL].
    pub const NAMES: [&'static str; 2] = [Self::Vector.name(), Self::VByte.name()];

    /// The name of this kind, as accepted by [str::parse] and printed by [std::fmt::Display].
    pub const fn name(self) -> &'static str {
        match self {
            Self::Vector => "vector",
            Self::VByte => "vbyte",
        }
    }

    /// The byte that identifies this kind in a serialized [PrefixTableVariant].
    const fn tag(self) -> u8 {
        match self {
            Self::Vector => 0,
            Self::VByte => 1,
        }
    }

    fn from_tag(tag: u8) -> std::io::Result<Self> {
        Self::ALL.into_iter().find(|kind| kind.tag() == tag).ok_or_else(|| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("Unknown prefix lookup table type tag {tag}"),
            )
        })
    }
}

impl std::fmt::Display for PrefixTableKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.name())
    }
}

impl std::str::FromStr for PrefixTableKind {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, String> {
        Self::ALL.into_iter().find(|kind| kind.name() == s).ok_or_else(|| {
            format!(
                "unknown prefix table type '{s}' (expected one of: {})",
                Self::NAMES.join(", ")
            )
        })
    }
}

/// A prefix lookup table of any of the implementations in this module.
/// This is what [SbwtIndex] stores.
#[derive(Clone, Eq, PartialEq, Debug)]
pub enum PrefixTableVariant {
    Vector(VectorLookupTable),
    VByte(VByteLookupTable),
}

/// The name this type had before there were several table implementations.
pub type PrefixLookupTable = PrefixTableVariant;

/// Calls `t.$method$args` on the table inside whichever variant `$self` is.
macro_rules! forward {
    ($self:expr, $method:ident $args:tt) => {
        match $self {
            PrefixTableVariant::Vector(t) => t.$method $args,
            PrefixTableVariant::VByte(t) => t.$method $args,
        }
    };
}

impl PrefixTableVariant {
    /// Builds a table of the given kind for all 4^p p-mers of `sbwt`.
    pub fn new_with_kind<SS: SubsetSeq>(kind: PrefixTableKind, sbwt: &SbwtIndex<SS>, prefix_length: usize) -> Self {
        match kind {
            PrefixTableKind::Vector => Self::Vector(VectorLookupTable::new(sbwt, prefix_length)),
            PrefixTableKind::VByte => Self::VByte(VByteLookupTable::new(sbwt, prefix_length)),
        }
    }

    /// Builds a table of the default kind ([PrefixTableKind::default]) for all 4^p p-mers of `sbwt`.
    pub fn new<SS: SubsetSeq>(sbwt: &SbwtIndex<SS>, prefix_length: usize) -> Self {
        Self::new_with_kind(PrefixTableKind::default(), sbwt, prefix_length)
    }

    /// The trivial table of the given kind with p = 0, holding only the range `0..n_sets` of the empty string.
    pub fn new_empty_with_kind(kind: PrefixTableKind, n_sets: usize) -> Self {
        match kind {
            PrefixTableKind::Vector => Self::Vector(VectorLookupTable::new_empty(n_sets)),
            PrefixTableKind::VByte => Self::VByte(VByteLookupTable::new_empty(n_sets)),
        }
    }

    /// The trivial table of the default kind with p = 0. See [PrefixTableVariant::new_empty_with_kind].
    pub fn new_empty(n_sets: usize) -> Self {
        Self::new_empty_with_kind(PrefixTableKind::default(), n_sets)
    }

    /// Which implementation this table is.
    pub fn kind(&self) -> PrefixTableKind {
        match self {
            Self::Vector(_) => PrefixTableKind::Vector,
            Self::VByte(_) => PrefixTableKind::VByte,
        }
    }

    /// The SBWT interval of `prefix`, which must have length [PrefixTableVariant::prefix_length].
    /// Returns `0..0` if the prefix does not occur or contains a non-ACGT character.
    pub fn lookup(&self, prefix: &[u8]) -> Range<usize> {
        forward!(self, lookup(prefix))
    }

    /// The prefix length p.
    pub fn prefix_length(&self) -> usize {
        forward!(self, prefix_length())
    }

    /// The number of bytes the table took on disk when it was loaded, or 0 if it was built in memory.
    pub fn size_bytes(&self) -> usize {
        forward!(self, size_bytes())
    }

    /// Writes a tag byte identifying the kind, followed by the table itself.
    /// Returns the number of bytes written.
    pub fn serialize<W: Write>(&self, out: &mut W) -> std::io::Result<usize> {
        let n_written = write_counted(out, &[self.kind().tag()])?;
        Ok(n_written + forward!(self, serialize(out))?)
    }

    /// Loads a table written by [PrefixTableVariant::serialize].
    pub fn load<R: Read>(input: &mut R) -> std::io::Result<Self> {
        let mut tag = [0u8; 1];
        input.read_exact(&mut tag)?;
        Ok(match PrefixTableKind::from_tag(tag[0])? {
            PrefixTableKind::Vector => Self::Vector(VectorLookupTable::load(input)?),
            PrefixTableKind::VByte => Self::VByte(VByteLookupTable::load(input)?),
        })
    }
}

/// Writes all of `bytes` and returns how many were written.
/// Used instead of `util::write_bytes`, which reports 8 bytes more than it writes.
fn write_counted<W: Write>(out: &mut W, bytes: &[u8]) -> std::io::Result<usize> {
    out.write_all(bytes)?;
    Ok(bytes.len())
}

/// The index of a p-mer in a table of all 4^p p-mers, with the first character
/// most significant. Returns None if the p-mer contains a non-ACGT character.
fn pmer_index(pmer: &[u8]) -> Option<usize> {
    let mut idx = 0_usize;
    for &c in pmer {
        let char_idx = ACGT_TO_0123[c as usize];
        if char_idx == 255 {
            return None;
        }
        idx = (idx << 2) | char_idx as usize;
    }
    Some(idx)
}

/// The SBWT intervals of all 4^p p-mers, in the order given by [pmer_index].
/// Absent p-mers get `0..0`.
fn all_pmer_ranges<SS: SubsetSeq>(sbwt: &SbwtIndex<SS>, prefix_length: usize) -> Vec<Range<usize>> {
    if prefix_length == 0 {
        return vec![sbwt.interval_of_empty_string()];
    }

    let total_pmers = 4_usize.pow(prefix_length as u32);
    let mut ranges = vec![0..0; total_pmers];
    let mut pmer = Vec::with_capacity(prefix_length);

    for (x, range) in ranges.iter_mut().enumerate() {
        pmer.clear();
        for i in 0..prefix_length {
            let char_idx = (x >> (2 * (prefix_length - 1 - i))) & 0x3;
            pmer.push(DNA_ALPHABET[char_idx]);
        }

        if let Some(r) = sbwt.search(&pmer) {
            *range = r;
        }
    }

    ranges
}

/// A table storing the ranges of all 4^p p-mers in a flat `Vec`, serialized as raw
/// little-endian 64-bit integers.
#[derive(Clone, Debug)]
pub struct VectorLookupTable {
    /// `ranges[i]` is the interval of the p-mer with index i, first character most significant.
    /// If the p-mer does not exist in the SBWT, the range is `0..0`.
    pub ranges: Vec<Range<usize>>,

    /// Prefix length p.
    pub prefix_length: usize,

    /// Number of bytes the table took on disk when loaded, 0 if built in memory.
    pub size_bytes: usize,
}

/// Equality ignores `size_bytes`, which differs between a built table and the same table loaded from disk.
impl PartialEq for VectorLookupTable {
    fn eq(&self, other: &Self) -> bool {
        self.prefix_length == other.prefix_length && self.ranges == other.ranges
    }
}

impl Eq for VectorLookupTable {}

impl VectorLookupTable {
    // Clippy false positive. It's actually intended like this
    #[allow(clippy::single_range_in_vec_init)]
    pub fn new_empty(n_sets_in_sbwt: usize) -> Self {
        Self {
            ranges: vec![0..n_sets_in_sbwt],
            prefix_length: 0,
            size_bytes: 0,
        }
    }

    pub fn new<SS: SubsetSeq>(sbwt: &SbwtIndex<SS>, prefix_length: usize) -> Self {
        Self {
            ranges: all_pmer_ranges(sbwt, prefix_length),
            prefix_length,
            size_bytes: 0,
        }
    }

    pub fn lookup(&self, prefix: &[u8]) -> Range<usize> {
        assert_eq!(prefix.len(), self.prefix_length);
        match pmer_index(prefix) {
            Some(idx) => self.ranges[idx].clone(),
            None => 0..0, // Not a DNA character
        }
    }

    pub fn prefix_length(&self) -> usize {
        self.prefix_length
    }

    pub fn size_bytes(&self) -> usize {
        self.size_bytes
    }

    /// Writes the table and returns the number of bytes written.
    pub fn serialize<W: Write>(&self, out: &mut W) -> std::io::Result<usize> {
        let mut n_written = 0_usize;
        n_written += write_counted(out, &(self.prefix_length as u64).to_le_bytes())?;
        n_written += write_counted(out, &(self.ranges.len() as u64).to_le_bytes())?;

        for range in self.ranges.iter() {
            n_written += write_counted(out, &(range.start as u64).to_le_bytes())?;
            n_written += write_counted(out, &(range.end as u64).to_le_bytes())?;
        }

        Ok(n_written)
    }

    /// Loads a table written by [VectorLookupTable::serialize].
    pub fn load<R: Read>(input: &mut R) -> std::io::Result<Self> {
        let prefix_length = input.read_u64::<LittleEndian>()? as usize;
        let ranges_len = input.read_u64::<LittleEndian>()? as usize;
        let mut bytes_read = 16; // prefix_length + ranges_len

        let mut ranges = vec![0..0; ranges_len];
        for range in ranges.iter_mut() {
            let start = input.read_u64::<LittleEndian>()? as usize;
            let end = input.read_u64::<LittleEndian>()? as usize;
            *range = start..end;
        }
        bytes_read += ranges_len * 16; // 2 x u64 per range

        Ok(Self {
            ranges,
            prefix_length,
            size_bytes: bytes_read,
        })
    }
}

/// A table storing the ranges of all 4^p p-mers in a flat `Vec`, serialized with
/// VByte + delta encoding (each range as its start and its length).
#[derive(Clone, Debug)]
pub struct VByteLookupTable {
    /// `ranges[i]` is the interval of the p-mer with index i, first character most significant.
    /// If the p-mer does not exist in the SBWT, the range is `0..0`.
    pub ranges: Vec<Range<usize>>,

    /// Prefix length p.
    pub prefix_length: usize,

    /// Number of bytes the table took on disk when loaded, 0 if built in memory.
    pub size_bytes: usize,
}

/// Equality ignores `size_bytes`, which differs between a built table and the same table loaded from disk.
impl PartialEq for VByteLookupTable {
    fn eq(&self, other: &Self) -> bool {
        self.prefix_length == other.prefix_length && self.ranges == other.ranges
    }
}

impl Eq for VByteLookupTable {}

impl VByteLookupTable {
    // Clippy false positive. It's actually intended like this
    #[allow(clippy::single_range_in_vec_init)]
    pub fn new_empty(n_sets_in_sbwt: usize) -> Self {
        Self {
            ranges: vec![0..n_sets_in_sbwt],
            prefix_length: 0,
            size_bytes: 0,
        }
    }

    pub fn new<SS: SubsetSeq>(sbwt: &SbwtIndex<SS>, prefix_length: usize) -> Self {
        Self {
            ranges: all_pmer_ranges(sbwt, prefix_length),
            prefix_length,
            size_bytes: 0,
        }
    }

    pub fn lookup(&self, prefix: &[u8]) -> Range<usize> {
        assert_eq!(prefix.len(), self.prefix_length);
        match pmer_index(prefix) {
            Some(idx) => self.ranges[idx].clone(),
            None => 0..0, // Not a DNA character
        }
    }

    pub fn prefix_length(&self) -> usize {
        self.prefix_length
    }

    pub fn size_bytes(&self) -> usize {
        self.size_bytes
    }

    /// Writes the table and returns the number of bytes written.
    pub fn serialize<W: Write>(&self, out: &mut W) -> std::io::Result<usize> {
        let mut n_written = 0_usize;
        n_written += Self::vbyte_write(out, self.prefix_length as u64)?;
        n_written += Self::vbyte_write(out, self.ranges.len() as u64)?;

        for range in self.ranges.iter() {
            let start = range.start as u64;
            let len = (range.end - range.start) as u64; // Delta compression

            n_written += Self::vbyte_write(out, start)?;
            n_written += Self::vbyte_write(out, len)?;
        }

        Ok(n_written)
    }

    /// Loads a table written by [VByteLookupTable::serialize].
    pub fn load<R: Read>(input: &mut R) -> std::io::Result<Self> {
        let (prefix_length, n1) = Self::vbyte_read(input)?;
        let (ranges_len, n2) = Self::vbyte_read(input)?;
        let mut bytes_read = n1 + n2;

        let mut ranges = vec![0..0; ranges_len as usize];
        for range in ranges.iter_mut() {
            let (start, r1) = Self::vbyte_read(input)?;
            let (len, r2) = Self::vbyte_read(input)?;
            bytes_read += r1 + r2;
            let end = start + len; // Reconstruct end using stored delta
            *range = (start as usize)..(end as usize);
        }

        Ok(Self {
            ranges,
            prefix_length: prefix_length as usize,
            size_bytes: bytes_read,
        })
    }

    fn vbyte_write<W: Write>(out: &mut W, mut val: u64) -> std::io::Result<usize> {
        let mut n_written = 0;
        loop {
            let mut byte = (val & 0x7F) as u8;
            val >>= 7;

            if val == 0 {
                byte |= 0x80;
            }

            n_written += write_counted(out, &[byte])?;

            if val == 0 {
                break;
            }
        }
        Ok(n_written)
    }

    fn vbyte_read<R: Read>(input: &mut R) -> std::io::Result<(u64, usize)> {
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
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kind_names_round_trip() {
        for kind in PrefixTableKind::ALL {
            assert_eq!(kind.name().parse::<PrefixTableKind>(), Ok(kind));
            assert_eq!(kind.to_string(), kind.name());
            assert_eq!(PrefixTableKind::from_tag(kind.tag()).unwrap(), kind);
        }
        assert!("nonsense".parse::<PrefixTableKind>().is_err());
        assert!(PrefixTableKind::from_tag(255).is_err());
    }

    #[test]
    fn pmer_index_order() {
        assert_eq!(pmer_index(b""), Some(0));
        assert_eq!(pmer_index(b"AA"), Some(0));
        assert_eq!(pmer_index(b"AC"), Some(1));
        assert_eq!(pmer_index(b"CA"), Some(4));
        assert_eq!(pmer_index(b"TT"), Some(15));
        assert_eq!(pmer_index(b"AN"), None);
    }

    #[test]
    fn serialize_and_load_each_kind() {
        let ranges = vec![0..0, 3..7, 7..7, 1000..123456];
        let tables = [
            PrefixTableVariant::Vector(VectorLookupTable { ranges: ranges.clone(), prefix_length: 1, size_bytes: 0 }),
            PrefixTableVariant::VByte(VByteLookupTable { ranges, prefix_length: 1, size_bytes: 0 }),
        ];
        for table in tables {
            let mut buf = Vec::<u8>::new();
            let n_written = table.serialize(&mut buf).unwrap();
            assert_eq!(n_written, buf.len());

            let loaded = PrefixTableVariant::load(&mut buf.as_slice()).unwrap();
            assert_eq!(loaded.kind(), table.kind());
            assert_eq!(loaded, table); // size_bytes is ignored by equality
            assert_eq!(loaded.size_bytes() + 1, n_written); // +1 for the tag byte
        }
    }
}