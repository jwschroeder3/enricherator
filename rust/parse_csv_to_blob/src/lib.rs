pub mod index;
pub mod parser;
pub mod config;
pub mod blob;

use std::error::Error;

pub type BoxedResult<T> = std::result::Result<T, Box<dyn Error>>;

#[derive(Hash, Eq, PartialEq, Debug, Clone)]
pub enum Parameter {
    Alpha,
    Beta,
    Gamma,
    non_sub,
}


