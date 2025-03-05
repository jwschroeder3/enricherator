use clap::Parser;
use serde::{Deserialize};
use serde_json;
use toml;
use std::fmt;
use std::fs::File;
use std::collections::{HashMap, HashSet};
use std::io::{self, Write, BufRead, BufReader, Read};
use std::error::Error;
use crate::BoxedResult;

#[allow(dead_code)]
#[cfg(test)]
mod tests {
    use super::*;
    //use assert_approx_eq::*;
    use std::ops::Deref;

    fn set_up_tests() -> HashMap<String, &'static str> {
        let into_colname_text = "#header text \nlp__,log_p__,log_g__,sig_noise.1,sub_Alpha.1.1.1,sub_Beta.1.1.1,sub_Gamma.1.1.1,Alpha.1.1.1,Alpha.2.1.1,Alpha.3.1.1,Alpha.1.1.2,Alpha.2.1.2,Alpha.3.1.2,Gamma.3.2.2\n";
        let into_sample_text = "#header text \nlp__,log_p__,log_g__,sig_noise.1,sub_Alpha.1.1.1,sub_Beta.1.1.1,sub_Gamma.1.1.1,Alpha.1.1.1,Alpha.2.1.1,Alpha.3.1.1,Alpha.1.1.2,Alpha.2.1.2,Alpha.3.1.2,Gamma.3.1.2\n#next header \n0.0,0.00,0.0,0.99999,1.0,1.2,1.2,1.0,2.0,3.0,4.0,5.0,6.0,7.0\n";
        let first_batch = "#header text \nlp__,log_p__,log_g__,sig_noise.1,sub_Alpha.1.1.1,sub_Beta.1";
        let leftover = "sub_Beta.1";
        let the_rest = ".1.1,sub_Gamma.1.1.1,Alpha.1.1.1,Alpha.2.1.1,Alpha.3.1.1,Alpha.1.1.2,Alpha.2.1.2,Alpha.3.1.2\n";
        let colname_text = "lp__,log_p__,log_g__,sig_noise.1,sub_Alpha.1.1.1,sub_Beta.1.1.1,sub_Gamma.1.1.1,Alpha.1.1.1,Alpha.2.1.1,Alpha.3.1.1,Alpha.1.1.2,Alpha.2.1.2,Alpha.3.1.2\n";
        let fname = "test_data/test_parse.csv";
        let expt_key = "test_data/expt_key.json";
        let parse_key = "test_data/parse_key.toml";
        let cfg_file = "test_data/config.toml";

        let contrast_text = "#header text \nlp__,log_p__,log_g__,sig_noise.1,sub_Alpha.1.1.1,sub_Beta.1.1.1,Alpha.1.1.1,Alpha.2.1.1,Alpha.3.1.1,Alpha.1.1.2,Alpha.2.1.2,Alpha.3.1.2,Beta.1.1.1,Beta.2.1.1,Beta.3.1.1,Beta.1.1.2,Beta.2.1.2,Beta.3.1.2,Gamma.1.1.1,Gamma.2.1.1,Gamma.3.1.1,Gamma.1.1.2,Gamma.2.1.2,Gamma.3.1.2\n#next header \n0.0,0.00,0.0,0.99999,1.0,1.2,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0\n0.0,0.00,0.0,0.99999,1.0,1.2,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5\n";

        let mut hash = HashMap::new();
        hash.insert(String::from("config"), cfg_file);
        hash.insert(String::from("into_colname"), into_colname_text);
        hash.insert(String::from("into_sample"), into_sample_text);
        hash.insert(String::from("first_batch"), first_batch);
        hash.insert(String::from("leftover"), leftover);
        hash.insert(String::from("the_rest"), the_rest);
        hash.insert(String::from("colname_text"), colname_text);
        hash.insert(String::from("fname"), fname);
        hash.insert(String::from("expt_key"), expt_key);
        hash.insert(String::from("parse_key"), parse_key);
        hash
    }

    #[test]
    fn test_read_config_file() {
        let hash = set_up_tests();
        let cfg = ConfigFile::read(
            hash.get("config").unwrap(),
        ).unwrap();
        assert_eq!(cfg.input.expt_key_file, String::from("expt_key.json"));
        assert_eq!(cfg.input.draws_file, String::from("draws.csv"));
        assert_eq!(cfg.input.pos_file_name, String::from("position_map.tsv"));
        assert_eq!(cfg.output.out_direc, String::from("out_direc"));
    }

    #[test]
    fn test_read_expt_key() {
        let hash = set_up_tests();
        let expt_key = ExptKey::read(
            hash.get("expt_key").unwrap()
        ).unwrap();
        assert_eq!(expt_key.var_map.get("Alpha").unwrap()[0], String::from("wt"));
        assert_eq!(expt_key.var_map.get("Alpha").unwrap()[1], String::from("vpsR"));
        assert_eq!(expt_key.var_map.get("Alpha").unwrap()[2], String::from("vpsT"));
        assert_eq!(expt_key.strand_map[0], String::from("plus"));
        assert_eq!(expt_key.strand_map[1], String::from("minus"));
    }
}

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct Args {
    /// toml file containing information about variables to retain and contrasts to perform, made
    /// by user
    #[arg(long)]
    pub cfg_file: String,
    /// json file from prep_data.R
    #[arg(long)]
    pub expt_key: String,
}

#[derive(Deserialize, Debug, Clone)]
struct ParserConfig {
    vars: Vec<String>,
    threshold: f64,
}

#[derive(Hash, Eq, PartialEq, Debug, Clone)]
enum Parameter {
    Alpha,
    Beta,
    Gamma,
    non_sub,
}

#[derive(Deserialize, Debug, Clone)]
pub struct ConfigFile {
    pub input: InputConfig,
    pub output: OutputConfig,
}

impl ConfigFile {
    pub fn read(fname: &str) -> BoxedResult<ConfigFile> {
        let mut file = File::open(fname)?;
        let mut data = String::new();
        file.read_to_string(&mut data)?;
        let cfg_file: ConfigFile = toml::from_str(&data)?;
        Ok(cfg_file)
    }
}

#[derive(Deserialize, Debug, Clone)]
pub struct InputConfig {
    pub expt_key_file: String,
    pub draws_file: String,
    pub pos_file_name: String,
}

#[derive(Deserialize, Debug, Clone)]
pub struct OutputConfig {
    pub out_direc: String,
}

#[derive(Deserialize, Debug)]
pub struct ExptKey {
    pub var_map: HashMap<String, Vec<String>>,
    pub strand_map: Vec<String>,
    pub resolution: Vec<u8>,
    pub pos_num: Vec<usize>,
    pub strand_num: Vec<usize>,
    pub sample_num: Vec<usize>,
}

impl ExptKey {
    pub fn read(fname: &str) -> BoxedResult<ExptKey> {
        let mut file = File::open(fname)?;
        let mut data = String::new();
        file.read_to_string(&mut data)?;
        let expt_key: ExptKey = serde_json::from_str(&data)?;
        Ok(expt_key)
    }
}

