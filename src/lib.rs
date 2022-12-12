use serde::{Deserialize, Serialize};

/// The non-negative cost of an alignment.
pub type Cost = i32;

/// A size, in bytes.
pub type Bytes = u64;

/// A single base
// NOTE: This is also part of rust-bio-types.
pub type Base = u8;

/// A sequence
// NOTE: This is also part of rust-bio-types.
pub type Sequence = Vec<Base>;

/// A non-owning sequence
pub type Seq<'a> = &'a [Base];

/// Takes an input string and returns the corresponding number of bytes. See the
/// documentation of the parse-size crates for details.
///
/// Use `K/M/G` for base `1000`, and `Ki/Mi/Gi` for base `2^10`. Case is
/// ignored. The trailing `B` is optional.
///
/// "1234"  => Bytes(1234)
/// "1 KB"  => Bytes(1000)
/// "1 kb"  => Bytes(1000)
/// "1 k"   => Bytes(1000)
/// "1 KiB" => Bytes(1024)
/// "1 kib" => Bytes(1024)
/// "1 ki"  => Bytes(1024)
/// ...
pub fn parse_bytes(input: &str) -> Result<Bytes, parse_size::Error> {
    parse_size::parse_size(input)
}

// TODO(ragnar): Define which direction is insertion and which is deletion.
#[derive(Serialize, Deserialize, Debug)]
pub enum CigarOp {
    Match,
    // TODO(ragnar): Choose between substitution and mismatch and use consistently.
    Sub,
    Del,
    Ins,
}

/// Types representation of a Cigar string.
// This is similar to https://docs.rs/bio/1.0.0/bio/alignment/struct.Alignment.html,
// but more specific for our use case.
#[derive(Serialize, Deserialize, Debug)]
pub struct Cigar {
    pub operations: Vec<(CigarOp, u32)>,
}

/// Different cost models.
/// All values must be non-negative.
// TODO(ragnar): Find a suitable name.
// TODO(ragnar): I am not sure of the best representation. This enum is
// conceptually nice, but possibly annoying in practice. Another option is to
// always have an in-memory representation in the most general way, and make
// additional constructors that fill fields for simpler variants. In this case
// we also need a way to go back from the general case to more specific cases
// via e.g. `is_unit()` and `is_linear()`.
#[derive(Serialize, Deserialize, Debug)]
pub enum CostModel {
    /// Levenshtein distance / Edit distance
    Unit,
    /// Different cost for substitutions and indels.
    Linear {
        /// >= 0
        r#match: Cost,
        /// > 0
        sub: Cost,
        /// > 0
        indel: Cost,
    },
    Affine {
        /// >= 0
        r#match: Cost,
        /// > 0
        sub: Cost,
        /// >= 0
        /// When 0, equivalent to Linear.
        open: Cost,
        /// > 0
        extend: Cost,
    },
}
