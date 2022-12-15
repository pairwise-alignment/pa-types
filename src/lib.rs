use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::fmt::Write;

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
#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, Clone, Copy)]
pub enum CigarOp {
    Match,
    Sub,
    Del,
    Ins,
}

impl CigarOp {
    pub fn to_char(&self) -> char {
        match self {
            CigarOp::Match => 'M',
            CigarOp::Sub => 'X',
            CigarOp::Ins => 'I',
            CigarOp::Del => 'D',
        }
    }
}

/// Types representation of a Cigar string.
// This is similar to https://docs.rs/bio/1.0.0/bio/alignment/struct.Alignment.html,
// but more specific for our use case.
#[derive(Serialize, Deserialize, Debug)]
pub struct Cigar {
    pub operations: Vec<(CigarOp, u32)>,
}

impl ToString for Cigar {
    fn to_string(&self) -> String {
        let mut s = String::new();
        for (op, cnt) in &self.operations {
            write!(&mut s, "{}{}", cnt, op.to_char()).unwrap();
        }
        s
    }
}

impl Cigar {
    pub fn from_ops(ops: impl Iterator<Item = CigarOp>) -> Self {
        Cigar {
            operations: ops
                .group_by(|&op| op)
                .into_iter()
                .map(|(op, group)| (op, group.count() as _))
                .collect(),
        }
    }

    pub fn push(&mut self, op: CigarOp) {
        if let Some(s) = self.operations.last_mut() {
            if s.0 == op {
                s.1 += 1;
                return;
            }
        }
        self.operations.push((op, 1));
    }

    pub fn push_matches(&mut self, cnt: usize) {
        if let Some(s) = self.operations.last_mut() {
            if s.0 == CigarOp::Match {
                s.1 += cnt as u32;
                return;
            }
        }
        self.operations.push((CigarOp::Match, cnt as _));
    }

    pub fn verify(&self, cm: &CostModel, a: Seq, b: Seq) -> Cost {
        let mut pos: (usize, usize) = (0, 0);
        let mut cost: Cost = 0;

        for &(op, cnt) in &self.operations {
            match op {
                CigarOp::Match => {
                    for _ in 0..cnt {
                        assert_eq!(a.get(pos.0), b.get(pos.1));
                        pos.0 += 1;
                        pos.1 += 1;
                    }
                }
                CigarOp::Sub => {
                    for _ in 0..cnt {
                        assert_ne!(a.get(pos.0), b.get(pos.1));
                        pos.0 += 1;
                        pos.1 += 1;
                        cost += cm.sub;
                    }
                }
                CigarOp::Ins => {
                    pos.1 += cnt as usize;
                    cost += cm.open + cnt as Cost * cm.extend;
                }
                CigarOp::Del => {
                    pos.0 += cnt as usize;
                    cost += cm.open + cnt as Cost * cm.extend;
                }
            }
        }
        assert!(pos == (a.len(), b.len()));

        cost
    }
}

/// Different cost models.
/// All values must be non-negative.
#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, Clone)]
pub struct CostModel {
    /// >= 0
    pub r#match: Cost,
    /// > 0
    pub sub: Cost,
    /// >= 0
    /// When 0, equivalent to Linear.
    pub open: Cost,
    /// > 0
    pub extend: Cost,
}

impl CostModel {
    pub fn unit() -> Self {
        Self {
            r#match: 0,
            sub: 1,
            open: 0,
            extend: 1,
        }
    }
    pub fn is_unit(&self) -> bool {
        self == &Self::unit()
    }
    pub fn linear(sub: Cost, indel: Cost) -> Self {
        Self {
            r#match: 0,
            sub,
            open: 0,
            extend: indel,
        }
    }
    pub fn is_linear(&self) -> bool {
        self.r#match == 0 && self.open == 0
    }
    pub fn affine(sub: Cost, open: Cost, extend: Cost) -> Self {
        Self {
            r#match: 0,
            sub,
            open,
            extend,
        }
    }
    pub fn is_affine(&self) -> bool {
        self.r#match == 0
    }
}
