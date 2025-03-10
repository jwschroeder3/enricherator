use clap::Parser;
use std::path::Path;
use std::fs::File;
use std::io::Write;
use std::error::Error;
use std::collections::HashMap;
use std::fmt;
use serde_json;
use serde::{Serialize, Deserialize};
use crate::{BoxedResult, Parameter};
use crate::config::ExptKey;

#[allow(dead_code)]
#[cfg(test)]
mod tests {
    use super::*;
}

/// Contains information for serializing/deserializing an index to
/// quickly access the blob.
/// Structure of blob (each "_") is an f32 sample (4 bytes) from approx posterior
///
/// p: param
/// g: genotype
/// q: strand
/// l: position
/// s: sample number = 16 in this example
///| p1, g1, q1, l1 | p1, g1, q1, l2 | p2, g1, q1, l1 | p2, g1, q1, l2 |
///|________________|________________|________________|________________|
///
/// So we need to write p * g * q * l * s * 4 bytes, plus 3 * s * 4 bytes for cubic relationship
///
/// param offset is G * Q * L * S * 4
/// geno offset is Q * L * S * 4
/// strand offset is L * S * 4
/// position offset is S * 4
///
/// Above example:
///  P = 2
///  G = 1
///  Q = 1
///  L = 2
///  S = 16
///
#[derive(Serialize, Deserialize, Debug)]
pub struct Index {
    pub offsets: HashMap<String, Offset>,
    pub total_bytes: usize,
    pub buffer_size: usize,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct Offset {
    pub param_id: ParamID,
    pub offset: u64,
}

#[derive(Serialize, Deserialize, Hash, Ord, PartialOrd, Eq, PartialEq, Clone)]
pub struct ParamID {
    pub parameter: String,
    pub geno_idx: usize,
    pub strand_idx: Option<usize>,
}

// Implement the `fmt::Debug` trait for the `ParamID` struct
impl fmt::Debug for ParamID {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // You can use the `write!` macro to format and write output
        if let Some(strand_idx) = self.strand_idx {
            write!(f, "ParamID {{ parameter: {}, genotype: {}, strand: {} }}", self.parameter, self.geno_idx, strand_idx)
        } else {
            write!(f, "ParamID {{ parameter: {}, genotype: {} }}", self.parameter, self.geno_idx)
        }
    }
}

// Implement the `fmt::Display` trait for the `ParamID` struct
impl fmt::Display for ParamID {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if let Some(strand_idx) = self.strand_idx {
            write!(f, "{}_{}_{}", self.parameter, self.geno_idx, strand_idx)
        } else {
            write!(f, "{}_{}", self.parameter, self.geno_idx)
        }
    }
}

impl ParamID {
    pub fn to_string(&self) -> String {
        if let Some(strand_idx) = self.strand_idx {
           format!("{}.{}.{}", self.parameter, self.geno_idx, strand_idx)
        } else {
            format!("{}.{}", self.parameter, self.geno_idx)
        }
    }
}

impl Offset{
    fn new(param_id: ParamID, offset: u64) -> Self {
        Self{
            param_id,
            offset,
        }
    }
}

impl Index {
    pub fn new(expt_key: &ExptKey) -> BoxedResult<Self> {

        //let mut offset_vec: Vec<Offset> = Vec::new();
        let mut offset_map: HashMap<String, Offset> = HashMap::new();
        let var_map = &expt_key.var_map;
        let geno_map = &expt_key.geno_map;
        let strand_names = &expt_key.strand_map;
        let pos_num: u64 = expt_key.pos_num[0].try_into()?;
        let samp_num: u64 = expt_key.sample_num[0].try_into()?;
        let buffer_size: usize = (samp_num * 4).try_into()?;

        let mut offset_bytes: u64 = 0;
        for (param_idx, (param_name,param_genos)) in var_map.iter().zip(geno_map).enumerate() {
            for (geno_idx, geno_name) in param_genos.iter().enumerate() {
                for (strand_idx, strand_name) in strand_names.iter().enumerate() {
                    let param_id = ParamID{
                        parameter: param_name.to_string(),
                        geno_idx: geno_idx,
                        strand_idx: Some(strand_idx),
                    };
                    let offset = Offset::new(param_id.clone(), offset_bytes.clone());
                    //offset_vec.push(offset);
                    offset_map.insert(param_id.to_string(), offset);
                    // increment offset by
                    // number of genome positions
                    // times number of samples
                    // times 4 bytes for f32
                    offset_bytes += pos_num * samp_num * 4;
                }
            }
        }

        // each Gamma genotype will have an a_cubic, b_cubic, and c_cubic associated
        // with it, so set the offsets for those here. Strand is None, because
        // we don't 
        for (param_idx, (param_name,param_genos)) in var_map.iter().zip(geno_map).enumerate() {
            if param_name == "Gamma" {
                for (geno_idx,geno_name) in param_genos.iter().enumerate() {
                    for var in ["a_cubic", "b_cubic", "c_cubic"] {
                        let param_id = ParamID{
                            parameter: var.to_string(),
                            geno_idx: geno_idx,
                            strand_idx: None,
                        };
                        let offset = Offset::new(param_id.clone(), offset_bytes.clone());
                        //offset_vec.push(offset);
                        offset_map.insert(param_id.to_string(), offset);
                        // increment offset by number of samples times 4 bytes for f32
                        // no pos_num on RHS here
                        offset_bytes += samp_num * 4;
                    }
                }
            }
        }

        let total_bytes = offset_bytes.try_into()?;

        Ok(Self{offsets: offset_map, total_bytes: total_bytes, buffer_size: buffer_size})
    }

    pub fn to_file(&self, outdir: &str) -> BoxedResult<()> {
        let fname = Path::new(outdir).join("index.json");
        let mut file = File::create(fname)?;
        let index_json = serde_json::to_string_pretty(self)?;
        file.write_all(index_json.as_bytes())?;
        Ok(())
    }
}

