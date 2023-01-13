use serde::{Deserialize, Serialize};

use crate::I;

/// The non-negative cost of an alignment.
pub type Cost = i32;

/// The score of an alignment. Can be positive or negative.
pub type Score = i32;

/// Different cost models.
/// All values must be non-negative.
#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, Clone, Copy)]
pub struct CostModel {
    /// > 0
    ///
    /// TODO: Support -infinity to disallow subtitutions entirely.
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
            sub,
            open: 0,
            extend: indel,
        }
    }
    pub fn is_linear(&self) -> bool {
        self.open == 0
    }
    pub fn affine(sub: Cost, open: Cost, extend: Cost) -> Self {
        Self { sub, open, extend }
    }
    pub fn is_affine(&self) -> bool {
        self.open > 0
    }

    /// The cost of a substitution.
    pub fn sub(&self) -> Cost {
        self.sub
    }
    /// The cost of a substitution, or None if not allowed.
    pub fn maybe_sub(&self) -> Option<Cost> {
        Some(self.sub)
    }
    /// The cost of an insertion of given length.
    pub fn ins(&self, len: I) -> Cost {
        self.open + len * self.extend
    }
    /// The cost of a deletion of given length.
    pub fn del(&self, len: I) -> Cost {
        self.open + len * self.extend
    }
}

#[derive(Debug, Eq, PartialEq, Clone, Copy)]
pub struct ScoreModel {
    /// > 0
    pub r#match: Score,
    /// < 0
    pub sub: Score,
    /// <= 0
    pub open: Score,
    /// < 0
    pub extend: Score,

    pub factor: i32,
}

impl ScoreModel {
    const OFFSET: i32 = 1;

    pub fn from_costs(cm: CostModel) -> Self {
        let factor;
        if cm.sub > 2 && cm.extend > 1 {
            factor = 1;
        } else if cm.sub == 1 {
            factor = 3;
        } else {
            factor = 2;
        }

        Self {
            r#match: Self::OFFSET * 2,
            // < 0
            sub: -cm.sub * factor + Self::OFFSET * 2,
            // <= 0
            open: -cm.open * factor,
            // < 0
            extend: -cm.extend * factor + Self::OFFSET,
            factor,
        }
    }

    pub fn global_cost(&self, score: Score, a_len: usize, b_len: usize) -> Cost {
        let path_len = (a_len + b_len) as i32;
        let s = -score + path_len * Self::OFFSET;
        assert!(s % self.factor == 0);
        s / self.factor
    }
}
