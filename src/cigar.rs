use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::fmt::Write;

use crate::*;

// TODO(ragnar): Define which direction is insertion and which is deletion.
#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, Clone, Copy)]
pub enum CigarOp {
    Match,
    Sub,
    Del,
    Ins,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct CigarElem {
    pub op: CigarOp,
    pub cnt: I,
}

impl CigarElem {
    pub fn new(op: CigarOp, cnt: I) -> Self {
        Self { op, cnt }
    }
}

/// Types representation of a Cigar string.
// This is similar to https://docs.rs/bio/1.0.0/bio/alignment/struct.Alignment.html,
// but more specific for our use case.
#[derive(Serialize, Deserialize, Debug)]
pub struct Cigar {
    pub ops: Vec<CigarElem>,
}

impl CigarOp {
    pub fn to_char(&self) -> char {
        match self {
            CigarOp::Match => '=',
            CigarOp::Sub => 'X',
            CigarOp::Ins => 'I',
            CigarOp::Del => 'D',
        }
    }
}
impl From<u8> for CigarOp {
    fn from(op: u8) -> Self {
        match op {
            b'=' | b'M' => CigarOp::Match,
            b'X' => CigarOp::Sub,
            b'I' => CigarOp::Ins,
            b'D' => CigarOp::Del,
            _ => panic!("Invalid CigarOp"),
        }
    }
}

impl ToString for Cigar {
    fn to_string(&self) -> String {
        let mut s = String::new();
        for elem in &self.ops {
            if elem.cnt == 1 {
                write!(&mut s, "{}", elem.op.to_char()).unwrap();
            } else {
                write!(&mut s, "{}{}", elem.cnt, elem.op.to_char()).unwrap();
            }
        }
        s
    }
}

impl Cigar {
    pub fn from_ops(ops: impl Iterator<Item = CigarOp>) -> Self {
        Cigar {
            ops: ops
                .group_by(|&op| op)
                .into_iter()
                .map(|(op, group)| CigarElem::new(op, group.count() as _))
                .collect(),
        }
    }

    pub fn from_path(a: Seq, b: Seq, path: &Path) -> Cigar {
        if path[0] != Pos(0, 0) {
            panic!("Path must start at (0,0)!");
        }
        Self::resolve_matches(
            path.iter().tuple_windows().map(|(&a, &b)| match b - a {
                Pos(0, 1) => CigarElem {
                    op: CigarOp::Ins,
                    cnt: 1,
                },
                Pos(1, 0) => CigarElem {
                    op: CigarOp::Del,
                    cnt: 1,
                },
                Pos(1, 1) => CigarElem {
                    op: CigarOp::Match,
                    cnt: 1,
                },
                _ => panic!("Path elements are not consecutive."),
            }),
            a,
            b,
        )
    }

    pub fn to_path(&self) -> Path {
        let mut pos = Pos(0, 0);
        let mut path = vec![pos];
        for el in &self.ops {
            for _ in 0..el.cnt {
                pos += match el.op {
                    CigarOp::Match => Pos(1, 1),
                    CigarOp::Sub => Pos(1, 1),
                    CigarOp::Del => Pos(1, 0),
                    CigarOp::Ins => Pos(0, 1),
                };
                path.push(pos);
            }
        }
        path
    }

    pub fn to_path_with_costs(&self, cm: CostModel) -> Vec<(Pos, Cost)> {
        let mut pos = Pos(0, 0);
        let mut cost = 0;
        let mut path = vec![(pos, cost)];

        for el in &self.ops {
            match el.op {
                CigarOp::Match => {
                    for _ in 0..(el.cnt as Cost) {
                        pos.0 += 1;
                        pos.1 += 1;
                        path.push((pos, cost));
                    }
                }
                CigarOp::Sub => {
                    for _ in 0..(el.cnt as Cost) {
                        pos.0 += 1;
                        pos.1 += 1;
                        cost += cm.sub;
                        path.push((pos, cost));
                    }
                }
                CigarOp::Ins => {
                    for len in 1..=(el.cnt as Cost) {
                        pos.1 += 1;
                        path.push((pos, cost + cm.ins(len)));
                    }
                    cost += cm.ins(el.cnt);
                }
                CigarOp::Del => {
                    for len in 1..=(el.cnt as Cost) {
                        pos.0 += 1;
                        path.push((pos, cost + cm.del(len)));
                    }
                    cost += cm.del(el.cnt);
                }
            }
        }
        path
    }

    pub fn push(&mut self, op: CigarOp) {
        if let Some(s) = self.ops.last_mut() {
            if s.op == op {
                s.cnt += 1;
                return;
            }
        }
        self.ops.push(CigarElem { op, cnt: 1 });
    }

    pub fn push_elem(&mut self, e: CigarElem) {
        if let Some(s) = self.ops.last_mut() {
            if s.op == e.op {
                s.cnt += e.cnt;
                return;
            }
        }
        self.ops.push(e);
    }

    pub fn push_matches(&mut self, cnt: I) {
        if let Some(s) = self.ops.last_mut() {
            if s.op == CigarOp::Match {
                s.cnt += cnt;
                return;
            }
        }
        self.ops.push(CigarElem {
            op: CigarOp::Match,
            cnt: cnt as _,
        });
    }

    pub fn verify(&self, cm: &CostModel, a: Seq, b: Seq) -> Cost {
        let mut pos: (usize, usize) = (0, 0);
        let mut cost: Cost = 0;

        for &CigarElem { op, cnt } in &self.ops {
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

    /// Splits all 'M'/Matches into matches and substitutions, and joins consecutive equal elements.
    pub fn resolve_matches(ops: impl Iterator<Item = CigarElem>, a: Seq, b: Seq) -> Self {
        let Pos(mut i, mut j) = Pos(0, 0);
        let mut c = Cigar { ops: vec![] };
        for CigarElem { op, cnt } in ops {
            match op {
                CigarOp::Match => {
                    for _ in 0..cnt {
                        c.push(if a[i as usize] == b[j as usize] {
                            CigarOp::Match
                        } else {
                            CigarOp::Sub
                        });
                        i += 1;
                        j += 1;
                    }
                    continue;
                }
                CigarOp::Sub => {
                    i += cnt;
                    j += cnt;
                }
                CigarOp::Ins => {
                    j += cnt;
                }
                CigarOp::Del => {
                    i += cnt;
                }
            };
            c.push_elem(CigarElem { op, cnt });
        }
        c
    }

    /// A simpler parsing function that only parses strings of characters `M=XID`, without preceding counts.
    /// Consecutive characters are grouped, and `M` and `=` chars are resolved into `=` and `X`.
    pub fn parse_without_counts(s: &str, a: Seq, b: Seq) -> Self {
        Self::resolve_matches(
            s.as_bytes().iter().map(|&op| CigarElem {
                op: op.into(),
                cnt: 1,
            }),
            a,
            b,
        )
    }

    /// A simpler parsing function that only parses strings of characters `MXID`, without preceding counts.
    /// Consecutive characters are grouped. `M` chars are *not* resolved and assumed to mean `=`.
    pub fn parse_without_resolving(s: &str) -> Self {
        let mut c = Cigar { ops: vec![] };
        for &op in s.as_bytes() {
            c.push(op.into())
        }
        c
    }

    /// A more generic (and slower) parsing function that also allows optional counts, e.g. `5M2X3M`.
    /// Consecutive characters are grouped, and `M` and `=` chars are resolved into `=` and `X`.
    pub fn parse(s: &str, a: Seq, b: Seq) -> Self {
        Self::resolve_matches(
            s.as_bytes()
                .split_inclusive(|b| b.is_ascii_alphabetic())
                .map(|slice| {
                    let (&op, cnt) = slice.split_last().unwrap();
                    let cnt = if cnt.is_empty() {
                        1
                    } else {
                        unsafe { std::str::from_utf8_unchecked(cnt) }
                            .parse()
                            .unwrap()
                    };
                    CigarElem { op: op.into(), cnt }
                }),
            a,
            b,
        )
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn to_string() {
        let c = Cigar {
            ops: vec![
                CigarElem {
                    op: CigarOp::Ins,
                    cnt: 1,
                },
                CigarElem {
                    op: CigarOp::Match,
                    cnt: 2,
                },
            ],
        };
        assert_eq!(c.to_string(), "I2=");
    }

    #[test]
    fn from_path() {
        let c = Cigar::from_path(
            b"aaa",
            b"aabc",
            &vec![Pos(0, 0), Pos(1, 1), Pos(2, 2), Pos(3, 3), Pos(3, 4)],
        );
        assert_eq!(c.to_string(), "2=XI");
    }
}
