pub mod index;
pub mod parser;
pub mod config;
pub mod blob;

use std::error::Error;
use std::fmt;

pub type BoxedResult<T> = std::result::Result<T, Box<dyn Error>>;

#[derive(Hash, Eq, PartialEq, Debug, Clone)]
pub enum Parameter {
    Alpha,
    Beta,
    Gamma,
    non_sub,
}

// Implement the `fmt::Display` trait for the `Parameter` enum variants
impl fmt::Display for Parameter {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Parameter::Alpha => {write!(f, "Alpha")},
            Parameter::Beta => {write!(f, "Beta")},
            Parameter::Gamma => {write!(f, "Gamma")},
            Parameter::non_sub => {write!(f, "non_sub")},
        }
    }
}

