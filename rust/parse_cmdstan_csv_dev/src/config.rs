use std::str::FromStr;
use strum_macros::EnumString;
use clap::Parser;
use serde::{Deserialize};
use serde_json;
use toml;
use std::fmt;
use std::fs::File;
use std::collections::{HashMap, HashSet};
use std::io::{self, Write, BufRead, BufReader, Read};
use std::error::Error;
use ndarray::{Array1, Array2};
use ndarray::prelude::*;

type BoxedResult<T> = std::result::Result<T, Box<dyn Error>>;

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
    fn test_strand_contrast_parser() {
        let hash = set_up_tests();
        let expt_key = ExptKey::read(hash.get("expt_key").unwrap()).unwrap();
        let mut contrast = StrandContrast::new("plus - minus", Parameter::Alpha, 50, 500);

        assert_eq!(contrast.parameter, Parameter::Alpha);
        assert_eq!(contrast.text, "plus - minus");

        contrast.classify_substrings("Alpha", &expt_key).unwrap();
        //println!("plus: {:?}", &contrast.numerators);
        //println!("minus: {:?}", &contrast.denominators);

        let plus_idx_target = vec![0];
        let minus_idx_target = vec![1];

        assert_eq!(contrast.parameter, Parameter::Alpha);
        assert_eq!(contrast.text, "plus - minus");
        assert_eq!(contrast.contrast_str, "plus-minus".to_string());
        assert_eq!(plus_idx_target[0], contrast.numerators.clone()[0]);
        assert_eq!(plus_idx_target.len(), contrast.numerators.clone().len());
        assert_eq!(minus_idx_target[0], contrast.denominators.clone()[0]);
        assert_eq!(minus_idx_target.len(), contrast.denominators.clone().len());

        let mut contrast = StrandContrast::new("minus - plus", Parameter::Alpha, 50, 500);

        assert_eq!(contrast.parameter, Parameter::Alpha);
        assert_eq!(contrast.text, "minus - plus");

        contrast.classify_substrings("Alpha", &expt_key).unwrap();
        //println!("plus: {:?}", &contrast.numerators);
        //println!("minus: {:?}", &contrast.denominators);

        let plus_idx_target = vec![1];
        let minus_idx_target = vec![0];

        assert_eq!(contrast.parameter, Parameter::Alpha);
        assert_eq!(contrast.text, "minus - plus");
        assert_eq!(contrast.contrast_str, "minus-plus".to_string());
        assert_eq!(plus_idx_target[0], contrast.numerators.clone()[0]);
        assert_eq!(plus_idx_target.len(), contrast.numerators.clone().len());
        assert_eq!(minus_idx_target[0], contrast.denominators.clone()[0]);
        assert_eq!(minus_idx_target.len(), contrast.denominators.clone().len());
 
    }

    #[test]
    fn test_geno_contrast_parser() {
        let hash = set_up_tests();
        let expt_key = ExptKey::read(hash.get("expt_key").unwrap()).unwrap();
        let pos_num = expt_key.pos_num[0];
        let mut contrast = GenoContrast::new("-wt - vpsR+vpsT - wt", Parameter::Gamma, pos_num, 500);
        contrast.classify_substrings("Gamma", &expt_key).unwrap();
        println!("contrast: {:?}", &contrast);

        let plus_idx_target = vec![2];
        let minus_idx_target = vec![0,1,0];

        assert_eq!(contrast.text, "-wt - vpsR+vpsT - wt");
        assert_eq!(contrast.contrast_str, "-wt-vpsR+vpsT-wt".to_string());
        assert_eq!(contrast.parameter, Parameter::Gamma);
        assert_eq!(plus_idx_target[0], contrast.numerators.clone()[0]);
        assert_eq!(plus_idx_target.len(), contrast.numerators.clone().len());
        assert_eq!(minus_idx_target[0], contrast.denominators.clone()[0]);
        assert_eq!(minus_idx_target[1], contrast.denominators.clone()[1]);
        assert_eq!(minus_idx_target[2], contrast.denominators.clone()[2]);
        assert_eq!(minus_idx_target.len(), contrast.denominators.clone().len());
        let mut contrast = GenoContrast::new("-wt - vpsR+vpsT - wt + sdfsdfd", Parameter::Gamma, pos_num, 500);
    }

    #[test]
    fn test_prep_config() {
        let hash = set_up_tests();
        let expt_key = ExptKey::read(
            hash.get("expt_key").unwrap()
        ).unwrap();
        let cfg = Config::from_file(
            hash.get("config").unwrap(),
            &expt_key,
        ).unwrap();
        assert_eq!(cfg.input.expt_key_file, String::from("expt_key.json"));
        assert_eq!(cfg.input.draws_file, String::from("draws.csv"));
        assert_eq!(cfg.input.pos_file_name, String::from("position_map.tsv"));
        assert_eq!(cfg.output.out_direc, String::from("out_direc"));
        assert_eq!(cfg.parse.vars[0], String::from("Alpha"));
        assert_eq!(cfg.parse.vars[3], String::from("non_sub"));
        assert_eq!(cfg.parse.threshold, 1.0);

        let geno_contrasts = cfg.geno_contrasts.clone().unwrap();
        let strand_contrasts = cfg.strand_contrasts.clone().unwrap();
        assert_eq!(geno_contrasts.threshold, 1.0);
        assert_eq!(strand_contrasts.threshold, 0.0);
        assert_eq!(
            geno_contrasts
                .var_key.get(&Parameter::Beta).unwrap()[0]
                .numerators.clone()[0],
            1
        );
        assert_eq!(
            geno_contrasts
                .var_key.get(&Parameter::Beta).unwrap()[0]
                .denominators.clone()[0],
            0
        );
        assert_eq!(
            strand_contrasts
                .var_key.get(&Parameter::Gamma).unwrap()[1]
                .denominators.clone()[0],
            0
        );
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
        assert_eq!(cfg.parse.vars[0], String::from("Alpha"));
        assert_eq!(cfg.parse.vars[3], String::from("non_sub"));
        assert_eq!(cfg.parse.threshold, 1.0);
        assert_eq!(cfg.geno_contrasts.clone().unwrap().Beta.unwrap()[0], String::from("vpsR-wt"));
        assert_eq!(cfg.geno_contrasts.clone().unwrap().Beta.unwrap()[2], String::from("vpsR - wt -vpsT + wt"));
        assert_eq!(cfg.geno_contrasts.clone().unwrap().Gamma.unwrap()[0], String::from("vpsR-wt"));
        assert_eq!(cfg.geno_contrasts.clone().unwrap().threshold, Some(1.0));
        assert_eq!(cfg.strand_contrasts.clone().unwrap().Gamma.unwrap()[0], String::from("plus-minus"));
        assert_eq!(cfg.strand_contrasts.clone().unwrap().Gamma.unwrap()[1], String::from("minus-plus"));
        assert!(cfg.strand_contrasts.clone().unwrap().threshold.is_none());
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

fn check_strand_error(
        strand_idx: usize,
        strand_names: &Vec<String>,
        substring: String,
) -> Result<(), StrandNotFoundError> {
    if strand_idx == strand_names.len() {
        return Err(StrandNotFoundError{ strand_name: substring.to_string() })
    }
    Ok(())
}

fn check_geno_error(
        geno_idx: usize,
        geno_vec: Vec<String>,
        substring: String,
        parameter: String,
) -> Result<(), GenotypeNotFoundError> {
    if geno_idx == geno_vec.len() {
        return Err(GenotypeNotFoundError{
            substring,
            parameter,
            geno_vec,
        })
    }
    Ok(())
}

// Define the custom error type
#[derive(Debug, Clone)]
pub struct GenotypeNotFoundError {
    substring: String,
    parameter: String,
    geno_vec: Vec<String>,
}

// Implement the Display trait for GenotypeNotFoundError for a user-friendly message
impl fmt::Display for GenotypeNotFoundError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "The genotype name '{}' in your genotype_contrast section of your config file is not valid. Allowed genotypes for parameter {:?} are {:?}.", self.substring, self.parameter, self.geno_vec)
    }
}

// Implement the Error trait for GenotypeNotFoundError
impl Error for GenotypeNotFoundError {}

// Define the custom error type
#[derive(Debug, Clone)]
pub struct StrandNotFoundError {
    strand_name: String,
}

// Implement the Display trait for StrandNotFoundError for a user-friendly message
impl fmt::Display for StrandNotFoundError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "The strand name '{}' in your strand_contrast section of your config file is not valid. Allowed strands are 'plus' and 'minus'.", self.strand_name)
    }
}

// Implement the Error trait for StrandNotFoundError
impl Error for StrandNotFoundError {}

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

impl ParserConfig {
    fn get_parse_instances(&self, expt_key: &ExptKey) -> BoxedResult<(HashSet<ParamID>, Vec<Parse>)> {
        let mut hs: HashSet<ParamID> = HashSet::new();
        let strand_vec = expt_key.get("strand_map").unwrap();
        // go through each parameter listed in expt_key to see whether we're keeping it.
        for (param_name, param_geno_vec) in ext_key.get("var_map").unwrap().iter() {
            // if the ParserConfig indicates we keep this parameter, go through its genotypes
            // and strands to set up the necessary Parse instances.
            if self.vars.contains(param_name) {
                let param = Parameter::from_str(param_name)?;
                for (geno_idx, geno_name) in self.vars.iter().enumerate() {
                    for (strand_idx, strand_name) in strand_vec.iter().enumerate() {
                        let param_id = ParamID::new(
                            param.clone(),
                            geno_idx.clone(),
                            strand_idx.clone(),
                        )
                        // a HashSet returns bool, true if inserted, false if not.
                        let _inserted = hs.insert(param_id);
                    }
                }
            }
        }
    }
}

#[derive(Debug)]
pub struct Parse {
    parameter: Parameter,
    string: String,
    param_effect_key: HashMap<>
    numerators: Vec<usize>,
    denominators: Vec<usize>,
    cl_array: Vec<Vec<f64>>,
    sum_array: Vec<f64>,
    sample_num: usize,
}

#[derive(Deserialize, Debug, Clone)]
struct ContrastsConfig {
    Alpha: Option<Vec<String>>,
    Beta: Option<Vec<String>>,
    Gamma: Option<Vec<String>>,
    non_sub: Option<Vec<String>>,
    threshold: Option<f64>,
}

#[derive(Hash, Eq, PartialEq, Debug, Clone, EnumString)]
enum Parameter {
    Alpha,
    Beta,
    Gamma,
    non_sub,
}

#[derive(Hash, Eq, PartialEq, Debug, Clone)]
struct ParamID {
    param: Parameter,
    geno_idx: usize,
    strand_idx: usize,
}

impl ParamID {
    fn new(param: Parameter, geno_idx: usize, strand_idx: usize) -> Self {
        Self {
            param,
            geno_idx,
            strand_idx,
        }
    }
}

/// A struct to allow tracking of each Parser or Genotype (anything implementing Tracker, really)
/// associated with a given set of parameters. Will be used during iteration over samples.
/// For a given cell in a draws.csv file, we'll check whether it's in the retained indices.
/// If so, look up all Trackers associated with this param and update their value
/// based on the action to take for this parameter.
///
/// Attributes:
///
/// key: HashMap<ParamID, Vec<Rc<impl Tracker>>>
///     Allows tracking of which Contrasts and Parsers make use of the ParamID.
///     Wrapping each Tracker in a Rc allows each Tracker
///     to be placed into potentially multiple ParamIDs.
#[derive(Debug, Clone)]
struct Trackers {
    var_key: HashMap<ParamID, Vec<Rc<impl Tracker>>>,
}

impl Trackers {
    fn new() -> Self {
        let hm: HashMap<ParamID, Vec<Rc<impl Tracker>>> = HashMap::new();
        Self{ var_key: hm }
    }

    fn insert(&mut self, param_id: ParamID, tracker: Rc<impl Tracker>) {
        if let Some(&mut current_tracker_list) = self.var_key.get_mut(&param_id) {
            current_tracker_list.push(Rc::clone(tracker));
        } else {
            let new_tracker_list: Vec<Rec<impl Tracker>> = vec![tracker];
            self.var_key.insert(param_id, new_tracker_list);
        }
    }

    fn contains_id(&self, param_id: &ParamID) -> bool {
        self.var_key.contains_key(param_id)
    }
}

trait Tracker {
    pub fn update_cl(&mut self) {}
    pub fn update_sum(&mut self) {}
}

/// A struct for building genotype contrast information
#[derive(Debug, Clone)]
struct GenoContrastsBuilder {}

impl GenoContrastsBuilder {
    fn new(cfg: &ContrastsConfig, expt_key: &ExptKey) -> BoxedResult<Self> {

        // instantiate a new mutable Trackers
        let mut trackers = Trackers::new();

        // one pass over parse, geno_contrasts, and strand_contrasts sections
        // to build set of ParamIDs.
        let mut param_id_set: HashSet<ParamID> = HashSet::new();
        

        // Build valid ParamID for each genotype/strand
        // if user wants alpha contrasts, build for alpha.
        if let Some(alpha) = &cfg.Alpha {
            // fetch alpha parameter's genotypes
            let param_genos = expt_key.get("var_map").unwrap()
            // for each genotype and strand, set up the ParamID, contrast, and place
            // ParamID into appropriate maps
            for (geno_idx, geno_name) in param_genos.iter().enumerate() {
                for (strand_idx, strand_name) in expt_key.get("strand_map").unwrap().iter().enumerate() {
                    let param_id = ParamID::new(Parameter::Alpha, *geno_idx, *strand_idx);
                    // is this param_id in those that we're using for any Parse or Contrast?
                    // if not, continue.
                    if param_id in 
                    // insert this tracker to associate it with its param_id.
                    trackers.insert(param_id, tracker);
                }
            }
        }




        let mut param_idx = 0;
        let pos_num = expt_key.pos_num[0];
        let sample_num = expt_key.sample_num[0];
        // if user wants alpha contrasts, build the HashMap for alpha.
        if let Some(alpha) = &cfg.Alpha {
            let mut inner_vec: Vec<GenoContrast> = Vec::new();
            for contrast_str in alpha {
                let mut contrast = GenoContrast::new(
                    contrast_str,
                    Parameter::Alpha,
                    pos_num.clone(),
                    sample_num.clone(),
                );
                contrast.classify_substrings("Alpha", &expt_key)?;
                inner_vec.push(contrast);
            }
            var_key.insert(Parameter::Alpha, inner_vec);
            param_idx += 1;
        }

        if let Some(beta) = &cfg.Beta {
            let mut inner_vec: Vec<GenoContrast> = Vec::new();
            for contrast_str in beta {
                let mut contrast = GenoContrast::new(
                    contrast_str,
                    Parameter::Beta,
                    pos_num.clone(),
                    sample_num.clone(),
                );
                contrast.classify_substrings("Beta", &expt_key)?;
                inner_vec.push(contrast);
            }
            var_key.insert(Parameter::Beta, inner_vec);
            param_idx += 1;
        }

        if let Some(gamma) = &cfg.Gamma {
            let mut inner_vec: Vec<GenoContrast> = Vec::new();
            for contrast_str in gamma {
                let mut contrast = GenoContrast::new(
                    contrast_str,
                    Parameter::Gamma,
                    pos_num.clone(),
                    sample_num.clone(),
                );
                contrast.classify_substrings("Gamma", &expt_key)?;
                inner_vec.push(contrast);
            }
            var_key.insert(Parameter::Gamma, inner_vec);
            param_idx += 1;

        }

        if let Some(non_sub) = &cfg.non_sub {
            let mut inner_vec: Vec<GenoContrast> = Vec::new();
            for contrast_str in non_sub {
                let mut contrast = GenoContrast::new(
                    contrast_str,
                    Parameter::non_sub,
                    pos_num.clone(),
                    sample_num.clone(),
                );
                contrast.classify_substrings("non_sub", &expt_key)?;
                inner_vec.push(contrast);
            }
            var_key.insert(Parameter::non_sub, inner_vec);
            param_idx += 1;
        }

        let threshold = if let Some(thresh) = &cfg.threshold {
            *thresh
        } else {
            0.0
        };
        Ok(Self{var_key, threshold})
    }
}

/// A struct for containing strand contrast information
///
/// Attributes:
///
/// var_key: HashMap<Parameter, Vec<StrandContrast>>
/// threshold: f64
///     Threshold for calculating evidence ratios. This threshold appies over
///     all contrasts
#[derive(Debug, Clone)]
struct StrandContrasts {
    var_key: HashMap<Parameter, Vec<StrandContrast>>,
    threshold: f64,
}

impl StrandContrasts {
    fn new(cfg: &ContrastsConfig, expt_key: &ExptKey) -> BoxedResult<Self> {
        let mut var_key: HashMap<Parameter, Vec<StrandContrast>>
            = HashMap::new();

        let mut param_idx = 0;
        let pos_num = expt_key.pos_num[0];
        let sample_num = expt_key.sample_num[0];
        // if user wants alpha contrasts, build the HashMap for alpha.
        if let Some(alpha) = &cfg.Alpha {
            let mut inner_vec: Vec<StrandContrast> = Vec::new();
            for contrast_str in alpha {
                let mut contrast = StrandContrast::new(
                    contrast_str,
                    Parameter::Alpha,
                    pos_num.clone(),
                    sample_num.clone(),
                );
                contrast.classify_substrings("Alpha", &expt_key)?;
                inner_vec.push(contrast);
            }
            var_key.insert(Parameter::Alpha, inner_vec);
            param_idx += 1;
        }

        if let Some(beta) = &cfg.Beta {
            let mut inner_vec: Vec<StrandContrast> = Vec::new();
            for contrast_str in beta {
                let mut contrast = StrandContrast::new(
                    contrast_str,
                    Parameter::Beta,
                    pos_num.clone(),
                    sample_num.clone(),
                );
                contrast.classify_substrings("Beta", &expt_key)?;
                inner_vec.push(contrast);
            }
            var_key.insert(Parameter::Beta, inner_vec);
            param_idx += 1;
        }

        if let Some(gamma) = &cfg.Gamma {
            let mut inner_vec: Vec<StrandContrast> = Vec::new();
            for contrast_str in gamma {
                let mut contrast = StrandContrast::new(
                    contrast_str,
                    Parameter::Gamma,
                    pos_num.clone(),
                    sample_num.clone(),
                );
                contrast.classify_substrings("Gamma", &expt_key)?;
                inner_vec.push(contrast);
            }
            var_key.insert(Parameter::Gamma, inner_vec);
            param_idx += 1;

        }

        if let Some(non_sub) = &cfg.non_sub {
            let mut inner_vec: Vec<StrandContrast> = Vec::new();
            for contrast_str in non_sub {
                let mut contrast = StrandContrast::new(
                    contrast_str,
                    Parameter::non_sub,
                    pos_num.clone(),
                    sample_num.clone(),
                );
                contrast.classify_substrings("non_sub", &expt_key)?;
                inner_vec.push(contrast);
            }
            var_key.insert(Parameter::non_sub, inner_vec);
            param_idx += 1;
        }

        let threshold = if let Some(thresh) = &cfg.threshold {
            *thresh
        } else {
            0.0
        };
        Ok(Self{var_key, threshold})
    }
}


/// The StrandContrast struct contains the information necessary to
/// retrieve proper values to perform the contrast of interest and
/// to write values to disk.
///
/// Attributes:
///
/// parameter: Parameter
///     Denotes whether this contrast is for Alpha, Beta, Gamma, or non_sub
/// text: String
///     The user's given text for the contrast, given in the config file.
/// name: String
///     The name of the contrast
/// numerators: Vec<usize>
///     Vec of usize. Each value is a strand index,
///     so that the proper values from this contrast's Parameter samples array
///     can be sliced and added to the result.
/// denominators: Vec<usize>
///     Same as numerators, but these values, when sliced from the Parameter samples array
///     will be subtracted from the result.
/// cl_array: Vec<Vec<f64>>
///     Outer vec is of length position_num, inners are of correct size to hold just the
///     samples that will be relevant to credible limit assignemnt plus one, so we can
///     insort new values and pop the middle one for every sample.
/// sum_array: Vec<f64>
///     Vec of length position_num to store cumulative sum of all samples for later
///     calculation of the mean of the samples.
/// sample_num: usize
///     Used to count number of samples in order to calculate mean.
#[derive(Debug, Clone)]
struct StrandContrast {
    parameter: Parameter,
    text: String,
    contrast_str: String,
    numerators: Vec<usize>,
    denominators: Vec<usize>,
    cl_array: Vec<Vec<f64>>,
    sum_array: Vec<f64>,
    sample_num: usize,
}

impl StrandContrast {
    fn new(text: &str, parameter: Parameter, position_num: usize, sample_num: usize) -> Self {
        let mut cl_outer: Vec<Vec<f64>> = Vec::with_capacity(position_num);
        let accumulator: Vec<f64> = Vec::with_capacity(position_num);
        let cl_num_needed = sample_num / 20 * 2 + 1;
        for p in 0..position_num {
            let cl_inner = Vec::with_capacity(cl_num_needed);
            cl_outer.push(cl_inner);
        }
        let numer_vec: Vec<usize> = Vec::new();
        let denom_vec: Vec<usize> = Vec::new();
        let cont_str = String::new();
        Self{
            parameter,
            text: text.to_string(),
            contrast_str: cont_str,
            numerators: numer_vec,
            denominators: denom_vec,
            cl_array: cl_outer,
            sum_array: accumulator,
            sample_num: sample_num,
        }
    }

    fn classify_substrings(
            &mut self,
            param_name: &str,
            expt_key: &ExptKey,
    ) -> BoxedResult<()> {

        let mut substring = "".to_string();
        let mut minus = false;
        let mut first_char = true;
        for c in self.text.chars() {
            if first_char {
                // if the first character is '-', set minus to true,
                // then first_char to false
                if c == '-' {
                    minus = true;
                // if the first character is space, continue without
                // switching first_char to false
                } else if c == ' ' {
                    continue
                // if the first character is '+' continue
                } else if c == '+' {
                    continue
                // if the first character is something else, push to string
                } else {
                    substring.push(c);
                }
                first_char = false;
                continue
            }
            // skip space characters
            if c == ' ' {
                continue
            // if we hit a '+' character, place current string into 
            // appropriate vec, then set minus to false and go to nect char
            } else if c == '+' || c == '-' {
                let mut strand_idx = 0;
                let mut unidentified_strand = false;

                let strand_names = &expt_key.strand_map;
                for strand_name in strand_names {
                    if *strand_name == substring {
                        break
                    }
                    strand_idx += 1;
                }
                // if we went through all strand names and didn't identify
                // the strand the user indicated, there's a problem, so error out
                check_strand_error(
                    strand_idx,
                    strand_names,
                    substring.clone(),
                )?;
                if minus {
                    self.denominators.push(strand_idx);
                    self.contrast_str.push('-');
                } else {
                    self.numerators.push(strand_idx);
                    if self.numerators.len() > 1 || self.denominators.len() > 0 {
                        self.contrast_str.push('+');
                    }
                }
                self.contrast_str.push_str(&substring);
                // after we've used substring, empty it.
                substring.clear();
                if c == '+' {
                    minus = false
                } else {
                    minus = true
                }
            } else {
                substring.push(c);
            }
        };

        let mut strand_idx = 0;
        let strand_names = &expt_key.strand_map;
        for strand_name in strand_names {
            if *strand_name == substring {
                break
            }
            strand_idx += 1;
        }
        // if we went through all strand names and didn't identify
        // the strand the user indicated, there's a problem, so error out
        check_strand_error(
            strand_idx,
            strand_names,
            substring.clone(),
        )?;

        if minus {
            self.denominators.push(strand_idx);
            self.contrast_str.push('-');
        } else {
            self.numerators.push(strand_idx);
            if self.numerators.len() > 1 || self.denominators.len() > 0 {
                self.contrast_str.push('+');
            }
        }
        self.contrast_str.push_str(&substring);

        Ok(())
    }
}

/// The Contrast struct contains the information necessary to
/// retrieve proper values to perform the contrast of interest and
/// to write values to disk.
///
/// Attributes:
///
/// parameter: Parameter
///     Denotes whether this contrast is for Alpha, Beta, Gamma, or non_sub
/// text: String
///     The user's given text for the contrast, given in the config file.
/// name: String
///     The name of the contrast
/// numerators: Vec<usize>
///     Each value gives a numerator's genotype index
///     so that the proper values from this contrast's Parameter samples array
///     can be sliced and added to the result.
/// denominators: Vec<usize>
///     Same as numerators, but these values, when sliced from the Parameter samples array
///     will be subtracted from the result.
/// cl_array: Vec<Vec<f64>>
///     Outer vec is of length position_num, inners are of correct size to hold just the
///     samples that will be relevant to credible limit assignemnt plus one, so we can
///     insort new values and pop the middle one for every sample.
/// sum_array: Vec<f64>
///     Vec of length position_num to store cumulative sum of all samples for later
///     calculation of the mean of the samples.
/// sample_num: usize
///     Used to count number of samples in order to calculate mean.
#[derive(Debug, Clone)]
struct GenoContrast {
    parameter: Parameter,
    text: String,
    contrast_str: String,
    numerators: Vec<usize>,
    denominators: Vec<usize>,
    cl_array: Vec<Vec<f64>>,
    sum_array: Vec<f64>,
    sample_num: usize,
}

enum VarType {
    Genotype,
    Strand,
}

impl GenoContrast {
    fn new(
            text: &str,
            parameter: Parameter,
            position_num: usize,
            sample_num: usize,
    ) -> Self {
        let mut cl_outer: Vec<Vec<f64>> = Vec::with_capacity(position_num);
        let accumulator: Vec<f64> = Vec::with_capacity(position_num);
        let cl_num_needed = sample_num / 20 * 2 + 1;
        for p in 0..position_num {
            let cl_inner = Vec::with_capacity(cl_num_needed);
            cl_outer.push(cl_inner);
        }
        let numer_vec: Vec<usize> = Vec::new();
        let denom_vec: Vec<usize> = Vec::new();
        let cont_str = String::new();
        Self{
            parameter,
            text: text.to_string(),
            contrast_str: cont_str,
            numerators: numer_vec,
            denominators: denom_vec,
            cl_array: cl_outer,
            sum_array: accumulator,
            sample_num: sample_num,
        }
    }

    fn classify_substrings(
            &mut self,
            param_name: &str,
            expt_key: &ExptKey,
    ) -> BoxedResult<()> {

        // Vectors to hold the substrings classified by their preceding character

        let mut substring = "".to_string();
        let mut minus = false;
        let mut first_char = true;
        for c in self.text.chars() {
            if first_char {
                // if the first character is '-', set minus to true,
                // then first_char to false
                if c == '-' {
                    minus = true;
                // if the first character is space, continue without
                // switching first_char to false
                } else if c == ' ' {
                    continue
                // if the first character is '+' continue
                } else if c == '+' {
                    continue
                // if the first character is something else, push to string
                } else {
                    substring.push(c);
                }
                first_char = false;
                continue
            }
            // skip space characters
            if c == ' ' {
                continue
            // if we hit a '+' character, place current string into 
            // appropriate vec, then set minus to false and go to nect char
            } else if c == '+' || c == '-' {
                let mut geno_idx = 0;
                let mut strand_idx = 0;

                // Will store each contrast as a Contrast struct
                let geno_names = expt_key.var_map.get(param_name).ok_or("No param_name")?;
                for geno_name in geno_names {
                    if *geno_name == substring {
                        break
                    }
                    geno_idx += 1;
                }
                check_geno_error(
                    geno_idx,
                    geno_names.to_vec(),
                    substring.clone(),
                    param_name.to_string(),
                )?;
                if minus {
                    self.denominators.push(geno_idx);
                    self.contrast_str.push('-');
                } else {
                    self.numerators.push(geno_idx);
                    if self.numerators.len() > 1 || self.denominators.len() > 0 {
                        self.contrast_str.push('+');
                    }
                }
                self.contrast_str.push_str(&substring);
                // after we've used substring, empty it.
                substring.clear();
                if c == '+' {
                    minus = false
                } else {
                    minus = true
                }
            } else {
                substring.push(c);
            }
        }

        let mut geno_idx = 0;
        let geno_names = expt_key.var_map.get(param_name).ok_or("No param_name")?;
        for geno_name in geno_names {
            if *geno_name == substring {
                break
            }
            geno_idx += 1;
        }
        check_geno_error(
            geno_idx,
            geno_names.to_vec(),
            substring.clone(),
            param_name.to_string(),
        )?;

        if minus {
            self.denominators.push(geno_idx);
            self.contrast_str.push('-');
        } else {
            self.numerators.push(geno_idx);
            if self.numerators.len() > 1 || self.denominators.len() > 0 {
                self.contrast_str.push('+');
            }
        }
        self.contrast_str.push_str(&substring);

        Ok(())
    }
}

#[derive(Deserialize, Debug, Clone)]
struct ConfigFile {
    input: InputConfig,
    output: OutputConfig,
    parse: ParserConfig,
    geno_contrasts: Option<ContrastsConfig>,
    strand_contrasts: Option<ContrastsConfig>,
}

impl ConfigFile {
    fn read(fname: &str) -> BoxedResult<ConfigFile> {
        let mut file = File::open(fname)?;
        let mut data = String::new();
        file.read_to_string(&mut data)?;
        let cfg_file: ConfigFile = toml::from_str(&data)?;
        Ok(cfg_file)
    }
}

#[derive(Debug)]
pub struct Config {
    pub input: InputConfig,
    pub output: OutputConfig,
}

impl Config {

    /// Config constructor
    ///
    /// Args:
    ///
    /// fname: &str
    ///     name of config file
    /// expt_key: &ExptKey
    ///     reference to ExptKey so that genotypes and strands can be matched
    ///     to their indices in the samples
    pub fn from_file(fname: &str, expt_key: &ExptKey) -> BoxedResult<Config> {
        let cfg_file = ConfigFile::read(&fname)?;

        let mut trackers = 

        // if there are geno_contrasts, fetch indices and construct Contrasts
        let geno_contrasts = if let Some(geno_cont) = cfg_file.geno_contrasts {
            let geno_con = GenoContrastsBuilder::new(&geno_cont, expt_key)?;
            Some(geno_con)
        } else {
            None
        };
        let strand_contrasts = if let Some(strand_cont) = cfg_file.strand_contrasts {
            let strand_con = StrandContrasts::new(&strand_cont, expt_key)?;
            Some(strand_con)
        } else {
            None
        };

        let parsers = Par ;

        let cfg = Config {
            input: cfg_file.input,
            output: cfg_file.output,
            parsers: parsers,
            geno_contrasts: geno_contrasts,
            strand_contrasts: strand_contrasts,
        };
        Ok((cfg, parsers, geno_contrasts, strand_contrasts))
    }
}

#[derive(Deserialize, Debug, Clone)]
struct InputConfig {
    expt_key_file: String,
    draws_file: String,
    pos_file_name: String,
}

#[derive(Deserialize, Debug, Clone)]
struct OutputConfig {
    out_direc: String,
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

