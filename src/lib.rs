use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Debug)]
pub struct Parameters {
    algo: Algorithm,
    /// Path to a `.seq` file.
    dataset: String,
    traceback: bool,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct ParametersList {
    list: Vec<Parameters>,
}

#[derive(Serialize, Debug)]
pub struct Results {
    params: Parameters,
    runtime_secs: Option<f64>,
    memory_bytes: Option<usize>,
    scores: Option<Vec<i32>>,
    cigar: Option<Vec<String>>,
}

#[derive(Serialize, Debug)]
pub struct ResultsList {
    list: Vec<Results>,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct Costs {
    /// Match cost >= 0.
    match_cost: i32,
    /// Mismatch cost < 0.
    mismatch_cost: i32,
    /// Gap open cost <= 0.
    gap_open: i32,
    /// Gap extend cost <= 0.
    gap_extend: i32,
}

#[derive(Serialize, Deserialize, Debug)]
pub enum Algorithm {
    BlockAligner {
        costs: Costs,
        min_size: usize,
        max_size: usize,
    },
    // Add more algorithms here!
}
