use std::{path::PathBuf, time::Duration};

use serde::{Deserialize, Serialize};

/// The non-negative cost of an alignment.
pub type Cost = i32;

/// A size, in bytes.
pub type Bytes = u64;

pub fn parse_bytes(input: &str) -> Result<Bytes, parse_size::Error> {
    parse_size::parse_size(input)
}

// TODO: Define which is insertion and which is deletion.
#[derive(Serialize, Deserialize, Debug)]
pub enum CigarOp {
    Match,
    Sub,
    Del,
    Ins,
}

// This is similar to https://docs.rs/bio/1.0.0/bio/alignment/struct.Alignment.html,
// but more specific for our use case.
#[derive(Serialize, Deserialize, Debug)]
pub struct Cigar {
    pub operations: Vec<(CigarOp, u32)>,
}

/// Different cost models.
/// All values must be non-negative.
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
        /// > 0
        sub: Cost,
        /// > 0
        indel: Cost,
    },
    Affine {
        /// > 0
        sub: Cost,
        /// >= 0
        /// When 0, equivalent to Linear.
        open: Cost,
        /// > 0
        extend: Cost,
    },
}

/// An alignment job: a single task for the runner to execute and benchmark.
#[derive(Serialize, Deserialize, Debug)]
pub struct Job {
    /// Path to a `.seq` file.
    pub dataset: PathBuf,
    /// The cost model to use.
    pub costs: CostModel,
    /// Return the full alignment/cigar?
    pub traceback: bool,

    /// The algorithm/parameters to use.
    pub algo: Algorithm,
}

/// The output of an alignment job.
#[derive(Serialize, Deserialize, Debug)]
pub struct JobOutput {
    pub runtime: Duration,
    pub memory: Bytes,
    pub costs: Vec<Cost>,
    pub cigars: Vec<Cigar>,
}

/// The result of an alignment job, containing the input and output.
#[derive(Serialize, Deserialize, Debug)]
pub struct JobResult {
    pub job: Job,
    pub output: Option<JobOutput>,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct Costs {
    /// Match cost >= 0.
    pub match_cost: i32,
    /// Mismatch cost < 0.
    pub mismatch_cost: i32,
    /// Gap open cost <= 0.
    pub gap_open: i32,
    /// Gap extend cost <= 0.
    pub gap_extend: i32,
}

#[derive(Serialize, Deserialize, Debug)]
pub enum Algorithm {
    BlockAligner {
        // TODO(ragnar): I think the costs are part of the input rather than the algorithm parameters.
        costs: Costs,
        min_size: usize,
        max_size: usize,
    },
    // Add more algorithms here!
}
