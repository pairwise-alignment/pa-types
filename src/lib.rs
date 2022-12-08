use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Debug)]
pub struct Parameters {
    pub algo: Algorithm,
    /// Path to a `.seq` file.
    pub dataset: String,
    pub traceback: bool,
}

#[derive(Serialize, Debug)]
pub struct Results {
    pub params: Parameters,
    pub runtime_secs: Option<f64>,
    pub memory_bytes: Option<usize>,
    pub scores: Option<Vec<i32>>,
    pub cigars: Option<Vec<String>>,
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
        costs: Costs,
        min_size: usize,
        max_size: usize,
    },
    // Add more algorithms here!
}
