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
    pub fn new(op: CigarOp, cnt: u32) -> Self {
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

impl ToString for Cigar {
    fn to_string(&self) -> String {
        let mut s = String::new();
        for elem in &self.ops {
            write!(&mut s, "{}{}", elem.cnt, elem.op.to_char()).unwrap();
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
                Pos(0, 1) => (CigarOp::Ins, 1),
                Pos(1, 0) => (CigarOp::Del, 1),
                Pos(1, 1) => (CigarOp::Match, 1),
                _ => panic!("Path elements are not consecutive."),
            }),
            a,
            b,
        )
    }

    pub fn to_path(&self) -> Path {
        let mut position = Pos(0, 0);
        let mut path = vec![position];
        for el in &self.ops {
            for _ in 0..el.cnt {
                position += match el.op {
                    CigarOp::Match => Pos(1, 1),
                    CigarOp::Sub => Pos(1, 1),
                    CigarOp::Del => Pos(1, 0),
                    CigarOp::Ins => Pos(0, 1),
                };
                path.push(position);
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
                    for _ in 0..(el.cnt as Cost) {
                        pos.1 += 1;
                        cost += cm.ins_cost(el.cnt);
                        path.push((pos, cost));
                    }
                }
                CigarOp::Del => {
                    for _ in 0..(el.cnt as Cost) {
                        pos.0 += 1;
                        cost += cm.del_cost(el.cnt);
                        path.push((pos, cost));
                    }
                }
            }
        }
        path
    }

    pub fn push(&mut self, op: CigarOp) {
        if let Some(s) = self.ops.last_mut() {
            if s.0 == op {
                s.1 += 1;
                return;
            }
        }
        self.ops.push((op, 1));
    }

    pub fn push_matches(&mut self, cnt: usize) {
        if let Some(s) = self.ops.last_mut() {
            if s.0 == CigarOp::Match {
                s.1 += cnt as u32;
                return;
            }
        }
        self.ops.push((CigarOp::Match, cnt as _));
    }

    pub fn verify(&self, cm: &CostModel, a: Seq, b: Seq) -> Cost {
        let mut pos: (usize, usize) = (0, 0);
        let mut cost: Cost = 0;

        for &(op, cnt) in &self.ops {
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

    /// Splits all 'M'/Matches into matches and substitutions.
    pub fn resolve_matches(ops: impl Iterator<Item = (CigarOp, u32)>, a: Seq, b: Seq) -> Self {
        let (mut i, mut j) = (0, 0);
        let mut operations = vec![];
        for (op, cnt) in ops {
            let cnt = cnt as usize;
            match op {
                CigarOp::Match => {
                    std::iter::zip(i..i + cnt, j..j + cnt)
                        .map(|(i, j)| a[i] == b[j])
                        .group_by(|&eq| eq)
                        .into_iter()
                        .for_each(|(eq, group)| {
                            operations.push((
                                if eq { CigarOp::Match } else { CigarOp::Sub },
                                group.count() as _,
                            ));
                        });

                    i += cnt;
                    j += cnt;
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
            operations.push((op, cnt as _));
        }
        Cigar { ops: operations }
    }

    pub fn parse(s: &str, a: Seq, b: Seq) -> Self {
        Self::resolve_matches(
            s.as_bytes()
                .split_inclusive(|b| b.is_ascii_alphabetic())
                .map(|slice| {
                    let (op, cnt) = slice.split_last().unwrap();
                    let cnt = if cnt.is_empty() {
                        1
                    } else {
                        unsafe { std::str::from_utf8_unchecked(cnt) }
                            .parse()
                            .unwrap()
                    };
                    let op = match *op {
                        b'M' | b'=' => CigarOp::Match,
                        b'X' => CigarOp::Sub,
                        b'I' => CigarOp::Ins,
                        b'D' => CigarOp::Del,
                        _ => panic!(),
                    };
                    (op, cnt as _)
                }),
            a,
            b,
        )
    }
}
