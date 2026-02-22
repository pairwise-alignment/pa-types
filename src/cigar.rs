/*!
+----------------+----------------+----------------------------------------------+---------------+---------------+
| Symbol         | Name           | Brief Description                            | Consumes Query| Consumes Ref  |
+----------------+----------------+----------------------------------------------+---------------+---------------+
| M              | Match          | No insertion or deletions, bases may differ | ✓             | ✓             |
| I              | Insertion      | Additional base in query (not in reference) | ✓             | ✗             |
| D              | Deletion       | Query is missing base from reference        | ✗             | ✓             |
| =              | Equal          | No insertions or deletions, bases agree     | ✓             | ✓             |
| X              | Not Equal      | No insertions or deletions, bases differ    | ✓             | ✓             |
| N              | None           | No query bases to align (spliced read)      | ✗             | ✓             |
| S              | Soft-Clipped   | Bases on end of read not aligned but stored | ✓             | ✗             |
| H              | Hard-Clipped   | Bases on end of read not aligned, not stored| ✗             | ✗             |
| P              | Padding        | Neither read nor reference has a base       | ✗             | ✗             |
+----------------+----------------+----------------------------------------------+---------------+---------------+

Uses SAM cigar conventions. See https://timd.one/blog/genomics/cigar.php
*/

use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::fmt::Write;

use crate::*;

#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, Clone, Copy)]
pub enum CigarOp {
    Match,
    Sub,
    Del,
    Ins,
}

#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, Clone, Copy)]
pub enum CigarOpChars {
    Match(u8),
    /// (from, to)
    Sub(u8, u8),
    Del(u8),
    Ins(u8),
}

#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, Clone, Copy)]
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
#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, Clone, Default)]
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

    /// Gives index change for operation (text_pos, query_pos)
    #[inline(always)]
    pub fn delta(&self) -> Pos {
        match self {
            CigarOp::Match | CigarOp::Sub => Pos(1, 1),
            CigarOp::Del => Pos(1, 0), // deletion: consumes ref/text only
            CigarOp::Ins => Pos(0, 1), // insertion: consumes pattern/query only
        }
    }

    /// Given index change (text_pos, query_pos) -> CigarOp
    /// Note ignores sub, (1,1) is always Match
    #[inline(always)]
    pub fn from_delta(delta: Pos) -> Self {
        match delta {
            Pos(0, 1) => CigarOp::Ins,
            Pos(1, 0) => CigarOp::Del,
            Pos(1, 1) => CigarOp::Match,
            _ => panic!("Invalid delta: {:?}", delta),
        }
    }
}

// Mostly to make verify/resolve matches clearer while still using the delta mapping
impl std::ops::Mul<I> for Pos {
    type Output = Pos;
    fn mul(self, rhs: I) -> Pos {
        Pos(self.0 * rhs, self.1 * rhs)
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
            write!(&mut s, "{}{}", elem.cnt, elem.op.to_char()).unwrap();
        }
        s
    }
}

impl Cigar {
    pub fn from_ops(ops: impl Iterator<Item = CigarOp>) -> Self {
        Cigar {
            ops: ops
                .chunk_by(|&op| op)
                .into_iter()
                .map(|(op, group)| CigarElem::new(op, group.count() as _))
                .collect(),
        }
    }

    /// Create Cigar from path and corresponding sequences.
    ///
    /// Assumes that path has (text_pos, pattern_pos) pairs.
    /// To distinguish between match and sub it uses simple
    /// equality (i.e. c1==c2) via `resolve_matches`
    pub fn from_path(text: Seq, pattern: Seq, path: &Path) -> Cigar {
        if path[0] != Pos(0, 0) {
            panic!("Path must start at (0,0)!");
        }
        Self::resolve_matches(
            path.iter()
                .tuple_windows()
                .map(|(&text_pos, &pattern_pos)| {
                    CigarElem::new(CigarOp::from_delta(pattern_pos - text_pos), 1)
                }),
            text,
            pattern,
        )
    }

    /// Return the diff from pattern to text.
    /// pos is always Pos(text_pos, query_pos).
    pub fn to_char_pairs<'s>(&'s self, text: &'s [u8], pattern: &'s [u8]) -> Vec<CigarOpChars> {
        let mut pos = Pos(0, 0);
        let fix_case = !(b'A' ^ b'a');
        let mut out = vec![];
        for el in &self.ops {
            for _ in 0..el.cnt {
                let c = match el.op {
                    CigarOp::Match => {
                        // NOTE: IUPAC characters can be matching even when they're not equal.
                        // assert_eq!(
                        //     (pattern[pos.0 as usize] & fix_case) as char,
                        //     (text[pos.1 as usize] & fix_case) as char,
                        //     "mismatch for {pos:?}"
                        // );
                        let c = CigarOpChars::Match(text[pos.0 as usize]);
                        pos += el.op.delta();
                        c
                    }
                    CigarOp::Sub => {
                        assert_ne!(
                            (text[pos.0 as usize] & fix_case) as char,
                            (pattern[pos.1 as usize] & fix_case) as char,
                            "cigar {:?}\npattern {:?}\ntext    {:?}\nmismatch for {pos:?}",
                            self.to_string(),
                            String::from_utf8_lossy(pattern),
                            String::from_utf8_lossy(text)
                        );
                        let c = CigarOpChars::Sub(text[pos.0 as usize], pattern[pos.1 as usize]);
                        pos += el.op.delta();
                        c
                    }
                    CigarOp::Del => {
                        // Note deletion consumes text hence text slice
                        let c = CigarOpChars::Del(text[pos.0 as usize]);
                        pos += el.op.delta();
                        c
                    }
                    CigarOp::Ins => {
                        // Note insertion consumes pattern hence pattern slice
                        let c = CigarOpChars::Ins(pattern[pos.1 as usize]);
                        pos += el.op.delta();
                        c
                    }
                };
                out.push(c);
            }
        }
        out
    }

    pub fn to_path(&self) -> Path {
        let mut pos = Pos(0, 0);
        let mut path = vec![pos];
        for el in &self.ops {
            for _ in 0..el.cnt {
                pos += el.op.delta();
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
                    for _ in 0..el.cnt {
                        pos += el.op.delta();
                        path.push((pos, cost));
                    }
                }
                CigarOp::Sub => {
                    for _ in 0..el.cnt {
                        pos += el.op.delta();
                        cost += cm.sub;
                        path.push((pos, cost));
                    }
                }
                CigarOp::Ins => {
                    for len in 1..=(el.cnt as Cost) {
                        pos += el.op.delta();
                        path.push((pos, cost + cm.ins(len)));
                    }
                    cost += cm.ins(el.cnt);
                }
                CigarOp::Del => {
                    for len in 1..=(el.cnt as Cost) {
                        pos += el.op.delta();
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

    pub fn verify(&self, cm: &CostModel, text: Seq, pattern: Seq) -> Result<Cost, &str> {
        let mut pos = Pos(0, 0);
        let mut cost: Cost = 0;

        for &CigarElem { op, cnt } in &self.ops {
            pos += op.delta() * cnt;
            match op {
                CigarOp::Match => {
                    for _ in 0..cnt {
                        if text.get(pos.0 as usize) != pattern.get(pos.1 as usize) {
                            return Err("Expected match but found substitution.");
                        }
                    }
                }
                CigarOp::Sub => {
                    for _ in 0..cnt {
                        if text.get(pos.0 as usize) == pattern.get(pos.1 as usize) {
                            return Err("Expected substitution but found match.");
                        }
                        cost += cm.sub;
                    }
                }
                CigarOp::Ins => {
                    cost += cm.open + cnt as Cost * cm.extend;
                }
                CigarOp::Del => {
                    cost += cm.open + cnt as Cost * cm.extend;
                }
            }
        }
        if pos != Pos(text.len() as I, pattern.len() as I) {
            return Err("Wrong alignment length.");
        }

        Ok(cost)
    }

    /// Splits all 'M'/Matches into matches and substitutions, and joins consecutive equal elements.
    pub fn resolve_matches(ops: impl Iterator<Item = CigarElem>, text: Seq, pattern: Seq) -> Self {
        let mut pos = Pos(0, 0);
        let mut c = Cigar { ops: vec![] };
        for CigarElem { op, cnt } in ops {
            match op {
                CigarOp::Match => {
                    for _ in 0..cnt {
                        c.push(if text[pos.0 as usize] == pattern[pos.1 as usize] {
                            CigarOp::Match
                        } else {
                            CigarOp::Sub
                        });
                        pos += op.delta();
                    }
                    continue;
                }
                _ => {
                    pos += op.delta() * cnt;
                }
            };
            c.push_elem(CigarElem { op, cnt });
        }
        c
    }

    /// A simpler parsing function that only parses strings of characters `M=XID`, without preceding counts.
    /// Consecutive characters are grouped, and `M` and `=` chars are resolved into `=` and `X`.
    pub fn parse_without_counts(s: &str, text: Seq, pattern: Seq) -> Self {
        Self::resolve_matches(
            s.as_bytes().iter().map(|&op| CigarElem {
                op: op.into(),
                cnt: 1,
            }),
            text,
            pattern,
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

    /// Parse a Cigar string with optional counts
    pub fn from_string(s: &str) -> Self {
        let mut c = Cigar { ops: vec![] };
        for slice in s.as_bytes().split_inclusive(|b| !b.is_ascii_digit()) {
            let (&op, cnt_bytes) = slice.split_last().expect("Cigar string cannot be empty");
            let cnt = if cnt_bytes.is_empty() {
                1
            } else {
                unsafe { std::str::from_utf8_unchecked(cnt_bytes) }
                    .parse()
                    .expect("Invalid Cigar count")
            };
            c.push_elem(CigarElem { op: op.into(), cnt });
        }
        c
    }

    /// A more generic (and slower) parsing function that also allows optional counts, e.g. `5M2X3M`.
    /// Consecutive characters are grouped, and `M` and `=` chars are resolved into `=` and `X`.
    pub fn parse(s: &str, text: Seq, pattern: Seq) -> Self {
        Self::resolve_matches(
            s.as_bytes()
                .split_inclusive(|pattern| pattern.is_ascii_alphabetic())
                .map(|pattern_slice| {
                    let (&op, cnt) = pattern_slice.split_last().unwrap();
                    let cnt = if cnt.is_empty() {
                        1
                    } else {
                        unsafe { std::str::from_utf8_unchecked(cnt) }
                            .parse()
                            .unwrap()
                    };
                    CigarElem { op: op.into(), cnt }
                }),
            text,
            pattern,
        )
    }

    /// Clear all cigar operations.
    pub fn clear(&mut self) {
        self.ops.clear();
    }

    pub fn reverse(&mut self) {
        self.ops.reverse();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_delta() {
        for op in [CigarOp::Match, CigarOp::Del, CigarOp::Ins] {
            assert_eq!(CigarOp::from_delta(op.delta()), op);
        }
        assert_eq!(CigarOp::from_delta(Pos(1, 1)), CigarOp::Match);
        //  assert_eq!(CigarOp::from_delta(Pos(1, 1)), CigarOp::Sub);
        assert_eq!(CigarOp::from_delta(Pos(1, 0)), CigarOp::Del); // Consume text
        assert_eq!(CigarOp::from_delta(Pos(0, 1)), CigarOp::Ins); // Consume pattern
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
        assert_eq!(c.to_string(), "1I2=");
    }

    #[test]
    fn from_path() {
        let c = Cigar::from_path(
            b"aaa",
            b"aabc",
            &vec![Pos(0, 0), Pos(1, 1), Pos(2, 2), Pos(3, 3), Pos(3, 4)],
        );
        assert_eq!(c.to_string(), "2=1X1I");
    }

    #[test]
    fn from_string_with_count() {
        let c = Cigar::from_string("24=");
        assert_eq!(c.ops.len(), 1);
        assert_eq!(c.ops[0], CigarElem::new(CigarOp::Match, 24));
    }

    #[test]
    fn from_string_mixed_ops() {
        let c = Cigar::from_string("2=3I1X");
        assert_eq!(
            c.ops,
            vec![
                CigarElem::new(CigarOp::Match, 2),
                CigarElem::new(CigarOp::Ins, 3),
                CigarElem::new(CigarOp::Sub, 1)
            ]
        );
    }

    #[test]
    fn from_string_no_counts() {
        let c = Cigar::from_string("=XIDDD");
        assert_eq!(
            c.ops,
            vec![
                CigarElem::new(CigarOp::Match, 1),
                CigarElem::new(CigarOp::Sub, 1),
                CigarElem::new(CigarOp::Ins, 1),
                CigarElem::new(CigarOp::Del, 3),
            ]
        );
    }

    #[test]
    #[rustfmt::skip]
    fn push_to_path() {
        let mut c = Cigar::default();
                                // 0 0
        c.push(CigarOp::Match); // 1  1
        c.push(CigarOp::Del);   // 2  1  (text +1)
        c.push(CigarOp::Ins);   // 2  2  (pattern +1)
        c.push(CigarOp::Sub);   // 3  3

        assert_eq!(
            c.to_path(),
            [
                Pos(0, 0),
                Pos(1, 1),
                Pos(2, 1),
                Pos(2, 2),
                Pos(3, 3),
            ]
        );
    }
}
