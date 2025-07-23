use clap::Parser;
use std::path::Path;
use serde::{Deserialize};
use serde_json;
use toml;
use std::fmt;
use std::fs::File;
use std::collections::{HashMap, HashSet};
use std::io::{self, Write, BufRead, BufReader, Read, BufWriter};
use std::error::Error;
use crate::{BoxedResult, Parameter, index::ParamID};

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
    fn test_prep_config() {
        let hash = set_up_tests();
        let expt_key = ExptKey::read(
            hash.get("expt_key").unwrap()
        ).unwrap();
        let cfg = Config::from_file(
            hash.get("config").unwrap(),
        ).unwrap();
        assert_eq!(cfg.input.expt_key_file, String::from("test_data/expt_key.json"));
        assert_eq!(cfg.input.draws_file, String::from("draws.csv"));
        assert_eq!(cfg.input.pos_file_name, String::from("position_map.tsv"));
        assert_eq!(cfg.output.out_direc, String::from("test_data/out_direc"));
        assert_eq!(cfg.parse.param_ids.len(), 24);
        assert_eq!(
            cfg.parse.param_ids[0].param_id,
            ParamID{
                parameter: "Alpha".to_string(),
                geno: "wt".to_string(),
                strand: Some("plus".to_string()),
            }
        );
        assert_eq!(
            cfg.parse.param_ids[3].param_id,
            ParamID{
                parameter: "Alpha".to_string(),
                geno: "vpsR".to_string(),
                strand: Some("minus".to_string()),
            }
        );
        assert_eq!(cfg.parse.threshold, 1.0);

        let geno_contrasts = cfg.geno_contrasts.unwrap();
        let strand_contrasts = cfg.strand_contrasts.unwrap();
        assert_eq!(geno_contrasts.threshold, 1.0);
        assert_eq!(strand_contrasts.threshold, 0.0);
        assert_eq!(
            geno_contrasts
                .contrasts[0]
                .numerators[0],
            ParamID{
                parameter: "Beta".to_string(),
                geno: "vpsR".to_string(),
                strand: Some("plus".to_string()),
            }
        );
        assert_eq!(
            geno_contrasts
                .contrasts[0]
                .denominators[0],
            ParamID{
                parameter: "Beta".to_string(),
                geno: "wt".to_string(),
                strand: Some("plus".to_string()),
            }
        );
        assert_eq!(
            geno_contrasts
                .contrasts[9]
                .numerators[0],
            ParamID{
                parameter: "Gamma".to_string(),
                geno: "vpsT".to_string(),
                strand: Some("minus".to_string()),
            }
        );
        assert_eq!(
            geno_contrasts
                .contrasts[9]
                .denominators[0],
            ParamID{
                parameter: "Gamma".to_string(),
                geno: "wt".to_string(),
                strand: Some("minus".to_string()),
            }
        );

    }

    #[test]
    fn test_read_config_file() {
        let hash = set_up_tests();
        let cfg = ConfigFile::read(
            hash.get("config").unwrap(),
        ).unwrap();
        assert_eq!(cfg.input.expt_key_file, String::from("test_data/expt_key.json"));
        assert_eq!(cfg.input.draws_file, String::from("draws.csv"));
        assert_eq!(cfg.input.pos_file_name, String::from("position_map.tsv"));
        assert_eq!(cfg.output.out_direc, String::from("test_data/out_direc"));
        assert_eq!(cfg.parse.vars[0], String::from("Alpha"));
        assert_eq!(cfg.parse.vars[3], String::from("non_sub"));
        assert_eq!(cfg.parse.threshold, Some(1.0));
        assert!(cfg.parse.rope.is_none());
        assert_eq!(cfg.geno_contrasts.clone().unwrap().Beta.unwrap()[0], String::from("vpsR-wt"));
        assert_eq!(cfg.geno_contrasts.clone().unwrap().Beta.unwrap()[2], String::from("vpsR - wt -vpsT + wt"));
        assert_eq!(cfg.geno_contrasts.clone().unwrap().Gamma.unwrap()[0], String::from("vpsR-wt"));
        assert_eq!(cfg.geno_contrasts.clone().unwrap().threshold, Some(1.0));
        assert_eq!(cfg.geno_contrasts.clone().unwrap().rope, Some(0.75));
        assert_eq!(cfg.strand_contrasts.clone().unwrap().Gamma.unwrap()[0], String::from("plus-minus"));
        assert_eq!(cfg.strand_contrasts.clone().unwrap().Gamma.unwrap()[1], String::from("minus-plus"));
        assert!(cfg.strand_contrasts.clone().unwrap().threshold.is_none());
        assert!(cfg.strand_contrasts.clone().unwrap().rope.is_none());
    }

    #[test]
    fn test_read_expt_key() {
        let hash = set_up_tests();
        let expt_key = ExptKey::read(
            hash.get("expt_key").unwrap()
        ).unwrap();
        assert_eq!(expt_key.var_map[0], "Alpha".to_string());
        assert_eq!(expt_key.strand_map[0], String::from("plus"));
        assert_eq!(expt_key.strand_map[1], String::from("minus"));
    }

    #[test]
    fn test_strand_contrast_parser() {
        let hash = set_up_tests();
        let expt_key = ExptKey::read(hash.get("expt_key").unwrap()).unwrap();
        let mut contrasts: Vec<Contrast> = Contrast::new(
            "plus - minus",
            Parameter::Alpha,
            "Alpha",
            &expt_key,
            ContrastType::Strand,
            "test_data/out_direc",
        ).unwrap();

        assert_eq!(contrasts[0].parameter, Parameter::Alpha);
        assert_eq!(contrasts[0].text, "plus - minus".to_string());
        assert_eq!(contrasts[0].contrast_str, "plus-minus".to_string());
        assert_eq!(
            ParamID{
                parameter: "Alpha".to_string(),
                geno: "wt".to_string(),
                strand: Some("plus".to_string()),
            },
            contrasts[0].numerators.clone()[0],
        );

        let mut contrasts = Contrast::new(
            "minus - plus",
            Parameter::Alpha,
            "Alpha",
            &expt_key,
            ContrastType::Strand,
            "test_data/out_direc",
        ).unwrap();

        assert_eq!(contrasts[0].parameter, Parameter::Alpha);
        assert_eq!(contrasts[0].text, "minus - plus");
        assert_eq!(contrasts[0].contrast_str, "minus-plus".to_string());
        assert_eq!(
            ParamID{
                parameter: "Alpha".to_string(),
                geno: "wt".to_string(),
                strand: Some("plus".to_string()),
            },
            contrasts[0].denominators.clone()[0],
        );
    }

    #[test]
    fn test_geno_contrast_parser() {
        let hash = set_up_tests();
        let expt_key = ExptKey::read(hash.get("expt_key").unwrap()).unwrap();
        let pos_num = expt_key.pos_num[0];
        let mut contrasts = Contrast::new(
            "-wt - vpsR+vpsT - wt",
            Parameter::Gamma,
            "Gamma",
            &expt_key,
            ContrastType::Genotype,
            "test_data/out_direc",
        ).unwrap();
        println!("contrasts: {:?}", &contrasts);

        assert_eq!(contrasts[0].parameter, Parameter::Gamma);
        assert_eq!(contrasts[0].text, "-wt - vpsR+vpsT - wt".to_string());
        assert_eq!(contrasts[0].contrast_str, "-wt-vpsR+vpsT-wt".to_string());
        assert_eq!(
            ParamID{
                parameter: "Gamma".to_string(),
                geno: "wt".to_string(),
                strand: Some("plus".to_string()),
            },
            contrasts[0].denominators.clone()[0],
        );
    }
}

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct Args {
    /// toml file containing information about variables to retain and contrasts to perform, made
    /// by user
    #[arg(long)]
    pub cfg_file: String,
}

#[derive(Deserialize)]
pub struct PositionMap {
    pub ctg: Vec<String>,
    pub start: Vec<usize>,
    pub end: Vec<usize>,
}

impl PositionMap {
    pub fn read(fname: &str) -> BoxedResult<Self> {
        let mut file = File::open(fname)?;
        let mut data = String::new();
        file.read_to_string(&mut data)?;
        let pos_map: Self = serde_json::from_str(&data)?;
        Ok(pos_map)
    }
}

#[derive(Deserialize, Debug, Clone)]
pub struct ConfigFile {
    pub input: InputConfig,
    pub output: OutputConfig,
    pub beta_sub_idx: Vec<usize>,
    pub parse: ParserConfig,
    pub geno_contrasts: Option<ContrastsConfig>,
    pub strand_contrasts: Option<ContrastsConfig>,
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

#[derive(Debug)]
pub struct Config {
    pub input: InputConfig,
    pub output: OutputConfig,
    pub beta_sub_idx: Vec<usize>,
    pub parse: ParamParsers,
    pub geno_contrasts: Option<Contrasts>,
    pub strand_contrasts: Option<Contrasts>,
}

#[derive(Debug)]
pub struct Contrasts {
    pub contrasts: Vec<Contrast>,
    pub threshold: f32,
    pub rope: f32,
}

#[derive(Debug)]
pub struct ParamParsers {
    pub param_ids: Vec<ParamParser>,
    pub threshold: f32,
    pub rope: f32,
}

#[derive(Debug)]
pub struct ParamParser {
    pub param_id: ParamID,
    pub low_file: BufWriter<File>,
    pub mean_file: BufWriter<File>,
    pub up_file: BufWriter<File>,
    pub Kge_file: BufWriter<File>,
    pub Klt_file: BufWriter<File>,
    pub K_in_rope_file: BufWriter<File>,
    pub K_not_in_rope_file: BufWriter<File>,
}

impl ParamParser {
    fn new(param_id: ParamID, out_dir: &str) -> BoxedResult<Self> {

        let mean_basename = format!("{}_mean.bedgraph", param_id);
        let mean_fname = Path::new(out_dir).join(mean_basename);
        let low_basename = format!("{}_lower.bedgraph", param_id);
        let low_fname = Path::new(out_dir).join(low_basename);
        let up_basename = format!("{}_upper.bedgraph", param_id);
        let up_fname = Path::new(out_dir).join(up_basename);
        let Kge_basename = format!("{}_Kge.bedgraph", param_id);
        let Kge_fname = Path::new(out_dir).join(Kge_basename);
        let Klt_basename = format!("{}_Klt.bedgraph", param_id);
        let Klt_fname = Path::new(out_dir).join(Klt_basename);
        let K_in_basename = format!("{}_K_in_rope.bedgraph", param_id);
        let K_in_fname = Path::new(out_dir).join(K_in_basename);
        let K_not_in_basename = format!("{}_K_not_in_rope.bedgraph", param_id);
        let K_not_in_fname = Path::new(out_dir).join(K_not_in_basename);

        let mean_f = File::create(mean_fname)?;
        let mean_file = BufWriter::new(mean_f);

        let low_f = File::create(low_fname)?;
        let low_file = BufWriter::new(low_f);

        let up_f = File::create(up_fname)?;
        let up_file = BufWriter::new(up_f);

        let Kge_f = File::create(Kge_fname)?;
        let Kge_file = BufWriter::new(Kge_f);

        let Klt_f = File::create(Klt_fname)?;
        let Klt_file = BufWriter::new(Klt_f);

        let K_in_f = File::create(K_in_fname)?;
        let K_in_rope_file = BufWriter::new(K_in_f);

        let K_not_in_f = File::create(K_not_in_fname)?;
        let K_not_in_rope_file = BufWriter::new(K_not_in_f);

        Ok(Self{
            param_id,
            low_file,
            mean_file,
            up_file,
            Kge_file,
            Klt_file,
            K_in_rope_file,
            K_not_in_rope_file,
        })
    }

    fn write_mean_to_file(
            &mut self,
            ctg: &str,
            start: &usize,
            end: &usize,
            val: f32,
    ) -> BoxedResult<()> {
        // Write the string content to the buffer
        let content = format!("{ctg}\t{start}\t{end}\t{val}\n");
        self.mean_file.write_all(content.as_bytes())?;
        Ok(())
    }

    fn flush_mean_buffer(&mut self) -> BoxedResult<()> {
        // Flush the buffer to ensure all data is written to the file
        self.mean_file.flush()?;
        Ok(())
    }

    fn write_lower_to_file(
            &mut self,
            ctg: &str,
            start: &usize,
            end: &usize,
            val: f32,
    ) -> BoxedResult<()> {
        // Write the string content to the buffer
        let content = format!("{ctg}\t{start}\t{end}\t{val}\n");
        self.low_file.write_all(content.as_bytes())?;
        Ok(())
    }

    fn flush_lower_buffer(&mut self) -> BoxedResult<()> {
        // Flush the buffer to ensure all data is written to the file
        self.low_file.flush()?;
        Ok(())
    }

    fn write_upper_to_file(
            &mut self,
            ctg: &str,
            start: &usize,
            end: &usize,
            val: f32,
    ) -> BoxedResult<()> {
        // Write the string content to the buffer
        let content = format!("{ctg}\t{start}\t{end}\t{val}\n");
        self.up_file.write_all(content.as_bytes())?;
        Ok(())
    }

    fn flush_upper_buffer(&mut self) -> BoxedResult<()> {
        // Flush the buffer to ensure all data is written to the file
        self.up_file.flush()?;
        Ok(())
    }

    fn write_Kge_to_file(
            &mut self,
            ctg: &str,
            start: &usize,
            end: &usize,
            val: f32,
    ) -> BoxedResult<()> {
        // Write the string content to the buffer
        let content = format!("{ctg}\t{start}\t{end}\t{val}\n");
        self.Kge_file.write_all(content.as_bytes())?;
        Ok(())
    }

    fn flush_Kge_buffer(&mut self) -> BoxedResult<()> {
        // Flush the buffer to ensure all data is written to the file
        self.Kge_file.flush()?;
        Ok(())
    }

    fn write_Klt_to_file(
            &mut self,
            ctg: &str,
            start: &usize,
            end: &usize,
            val: f32,
    ) -> BoxedResult<()> {
        // Write the string content to the buffer
        let content = format!("{ctg}\t{start}\t{end}\t{val}\n");
        self.Klt_file.write_all(content.as_bytes())?;
        Ok(())
    }

    fn flush_Klt_buffer(&mut self) -> BoxedResult<()> {
        // Flush the buffer to ensure all data is written to the file
        self.Klt_file.flush()?;
        Ok(())
    }

    fn write_K_in_to_file(
            &mut self,
            ctg: &str,
            start: &usize,
            end: &usize,
            val: f32,
    ) -> BoxedResult<()> {
        // Write the string content to the buffer
        let content = format!("{ctg}\t{start}\t{end}\t{val}\n");
        self.K_in_rope_file.write_all(content.as_bytes())?;
        Ok(())
    }

    fn flush_K_in_rope_buffer(&mut self) -> BoxedResult<()> {
        // Flush the buffer to ensure all data is written to the file
        self.K_in_rope_file.flush()?;
        Ok(())
    }

    fn write_K_not_in_to_file(
            &mut self,
            ctg: &str,
            start: &usize,
            end: &usize,
            val: f32,
    ) -> BoxedResult<()> {
        // Write the string content to the buffer
        let content = format!("{ctg}\t{start}\t{end}\t{val}\n");
        self.K_not_in_rope_file.write_all(content.as_bytes())?;
        Ok(())
    }

    fn flush_K_not_in_rope_buffer(&mut self) -> BoxedResult<()> {
        // Flush the buffer to ensure all data is written to the file
        self.K_not_in_rope_file.flush()?;
        Ok(())
    }

    pub fn flush_buffers(&mut self) -> BoxedResult<()> {
        self.flush_upper_buffer()?;
        self.flush_lower_buffer()?;
        self.flush_mean_buffer()?;
        self.flush_Kge_buffer()?;
        self.flush_Klt_buffer()?;
        self.flush_K_in_rope_buffer()?;
        self.flush_K_not_in_rope_buffer()?;
        Ok(())
    }

    pub fn write_vals_to_files(
            &mut self,
            ctg: &str,
            start: &usize,
            end: &usize,
            mean: f32,
            lower: f32,
            upper: f32,
            Kge: f32,
            Klt: f32,
            Kin: f32,
            Knotin: f32,
    ) -> BoxedResult<()> {
        self.write_mean_to_file(
            ctg,
            start,
            end,
            mean,
        )?;
        self.write_lower_to_file(
            ctg,
            start,
            end,
            lower,
        )?;
        self.write_upper_to_file(
            ctg,
            start,
            end,
            upper,
        )?;
        self.write_Kge_to_file(
            ctg,
            start,
            end,
            Kge,
        )?;
        self.write_Klt_to_file(
            ctg,
            start,
            end,
            Klt,
        )?;
        self.write_K_in_to_file(
            ctg,
            start,
            end,
            Kin,
        )?;
        self.write_K_not_in_to_file(
            ctg,
            start,
            end,
            Knotin,
        )?;
        Ok(())
    }
}

impl ParamParsers {
    fn new(parse_cfg: &ParserConfig, expt_key: &ExptKey, out_dir: &str) -> BoxedResult<Self> {
        let mut parsers: Vec<ParamParser> = Vec::new();
        for param_name in parse_cfg.vars.iter() {
            let param_idx = &expt_key.var_map.iter()
                .position(|s| s == param_name).ok_or("No such param_name")?;
            for geno_name in &expt_key.geno_map[*param_idx] {
                for strand_name in &expt_key.strand_map {
                    let param_id = ParamID{
                        parameter: param_name.to_string(),
                        geno: geno_name.to_string(),
                        strand: Some(strand_name.to_string()),
                    };
                    let par_parser = ParamParser::new(param_id, out_dir)?;
                    parsers.push(par_parser);
                }
            }
        }

        let threshold = if let Some(thresh) = parse_cfg.threshold {
            thresh
        } else {
            0.0
        };

        // set the default ROPE for parameters to log2-enrichment = abs(0.5)
        let rope = if let Some(rope_bound) = parse_cfg.rope {
            rope_bound
        } else {
            0.5
        };

        // Leave the rope as Option<value>, since some users won't set it.
        Ok(Self{param_ids: parsers, threshold: threshold, rope: rope})
    }
}

impl Config {
    pub fn from_file(fname: &str) -> BoxedResult<Config> {
        let cfg_file_res = ConfigFile::read(fname);
        let cfg_file = match cfg_file_res {
            Ok(file) => file,
            Err(error) => {
                println!("Error reading ConfigFile");
                println!("{:?}", error);
                return Err(error)
            }
        };
        let expt_key_res = ExptKey::read(&cfg_file.input.expt_key_file);
        let expt_key = match expt_key_res {
            Ok(file) => file,
            Err(error) => {
                println!("Error reading ExptKey");
                println!("{:?}", error);
                return Err(error)
            }
        };

        // if there are geno_contrasts, fetch indices and construct Contrasts
        let geno_contrasts = if let Some(geno_cont) = cfg_file.geno_contrasts {
            let geno_con = GenoContrastsBuilder::new(&geno_cont, &expt_key, &cfg_file.output.out_direc)?;
            Some(geno_con)
        } else {
            None
        };
        let strand_contrasts = if let Some(strand_cont) = cfg_file.strand_contrasts {
            let strand_con = StrandContrastsBuilder::new(&strand_cont, &expt_key, &cfg_file.output.out_direc)?;
            Some(strand_con)
        } else {
            None
        };

        let parsers = ParamParsers::new(&cfg_file.parse, &expt_key, &cfg_file.output.out_direc)?;

        let cfg = Config {
            input: cfg_file.input,
            output: cfg_file.output,
            beta_sub_idx: cfg_file.beta_sub_idx,
            parse: parsers,
            geno_contrasts: geno_contrasts,
            strand_contrasts: strand_contrasts,
        };
        Ok(cfg)
 
    }
}

/// A struct for building genotype contrast information
#[derive(Debug, Clone)]
struct StrandContrastsBuilder {}

impl StrandContrastsBuilder {
    fn new(cfg: &ContrastsConfig, expt_key: &ExptKey, out_dir: &str) -> BoxedResult<Contrasts> {

        let pos_num = expt_key.pos_num[0];
        let sample_num = expt_key.sample_num[0];
        let mut contrast_vec: Vec<Contrast> = Vec::new();
        // if user wants alpha contrasts, build the HashMap for alpha.
        if let Some(alpha) = &cfg.Alpha {
            for contrast_str in alpha {
                let mut contrasts = Contrast::new(
                    contrast_str,
                    Parameter::Alpha,
                    "Alpha",
                    expt_key,
                    ContrastType::Strand,
                    out_dir,
                )?;
                contrast_vec.append(&mut contrasts);
            }
        }

        if let Some(beta) = &cfg.Beta {
            for contrast_str in beta {
                let mut contrasts = Contrast::new(
                    contrast_str,
                    Parameter::Beta,
                    "Beta",
                    expt_key,
                    ContrastType::Strand,
                    out_dir,
                )?;
                contrast_vec.append(&mut contrasts);
            }
        }

        if let Some(gamma) = &cfg.Gamma {
            for contrast_str in gamma {
                let mut contrasts = Contrast::new(
                    contrast_str,
                    Parameter::Gamma,
                    "Gamma",
                    expt_key,
                    ContrastType::Strand,
                    out_dir,
                )?;
                contrast_vec.append(&mut contrasts);
            }

        }

        if let Some(non_sub) = &cfg.non_sub {
            for contrast_str in non_sub {
                let mut contrasts = Contrast::new(
                    contrast_str,
                    Parameter::non_sub,
                    "non_sub",
                    expt_key,
                    ContrastType::Strand,
                    out_dir,
                )?;
                contrast_vec.append(&mut contrasts);
            }
        }

        let threshold: f32 = if let Some(thresh) = &cfg.threshold {
            *thresh
        } else {
            0.0
        };
        
        // set the default ROPE for contrasts to log2(FC) = abs(0.5)
        let rope = if let Some(rope_bound) = &cfg.rope {
            *rope_bound
        } else {
            0.5
        };

        Ok(Contrasts{
            contrasts: contrast_vec,
            threshold: threshold,
            rope: rope,
        })
    }
}

/// A struct for building genotype contrast information
#[derive(Debug, Clone)]
struct GenoContrastsBuilder {}

impl GenoContrastsBuilder {
    fn new(cfg: &ContrastsConfig, expt_key: &ExptKey, out_dir: &str) -> BoxedResult<Contrasts> {

        let pos_num = expt_key.pos_num[0];
        let sample_num = expt_key.sample_num[0];
        let mut contrast_vec: Vec<Contrast> = Vec::new();
        // if user wants alpha contrasts, build the HashMap for alpha.
        if let Some(alpha) = &cfg.Alpha {
            for contrast_str in alpha {
                let mut contrasts = Contrast::new(
                    contrast_str,
                    Parameter::Alpha,
                    "Alpha",
                    expt_key,
                    ContrastType::Genotype,
                    out_dir,
                )?;
                contrast_vec.append(&mut contrasts);
            }
        }

        if let Some(beta) = &cfg.Beta {
            for contrast_str in beta {
                let mut contrasts = Contrast::new(
                    contrast_str,
                    Parameter::Beta,
                    "Beta",
                    expt_key,
                    ContrastType::Genotype,
                    out_dir,
                )?;
                contrast_vec.append(&mut contrasts);
            }
        }

        if let Some(gamma) = &cfg.Gamma {
            for contrast_str in gamma {
                let mut contrasts = Contrast::new(
                    contrast_str,
                    Parameter::Gamma,
                    "Gamma",
                    expt_key,
                    ContrastType::Genotype,
                    out_dir,
                )?;
                contrast_vec.append(&mut contrasts);
            }

        }

        if let Some(non_sub) = &cfg.non_sub {
            for contrast_str in non_sub {
                let mut contrasts = Contrast::new(
                    contrast_str,
                    Parameter::non_sub,
                    "non_sub",
                    expt_key,
                    ContrastType::Genotype,
                    out_dir,
                )?;
                contrast_vec.append(&mut contrasts);
            }
        }

        let threshold = if let Some(thresh) = &cfg.threshold {
            *thresh
        } else {
            0.0
        };

        // set the default ROPE for contrasts to log2(FC) = abs(0.5)
        let rope = if let Some(rope_bound) = &cfg.rope {
            *rope_bound
        } else {
            0.5
        };

        Ok(Contrasts{
            contrasts: contrast_vec,
            threshold: threshold,
            rope: rope,
        })
    }
}

#[derive(Deserialize, Debug, Clone)]
struct ParserConfig {
    vars: Vec<String>,
    threshold: Option<f32>,
    rope: Option<f32>,
}

#[derive(Deserialize, Debug, Clone)]
struct ContrastsConfig {
    Alpha: Option<Vec<String>>,
    Beta: Option<Vec<String>>,
    Gamma: Option<Vec<String>>,
    non_sub: Option<Vec<String>>,
    threshold: Option<f32>,
    rope: Option<f32>,
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
/// numerators: Vec<ParamID>
///     Each value gives a numerator's ParamID so that the proper values from the blob
///     can be read and added to the result.
/// denominators: Vec<ParamID>
///     Same as numerators, but these values, when read from the blob,
///     will be subtracted from the result.
#[derive(Debug)]
pub struct Contrast {
    pub parameter: Parameter,
    pub text: String,
    pub contrast_str: String,
    pub numerators: Vec<ParamID>,
    pub denominators: Vec<ParamID>,
    pub low_file: BufWriter<File>,
    pub mean_file: BufWriter<File>,
    pub up_file: BufWriter<File>,
    pub Kge_file: BufWriter<File>,
    pub Klt_file: BufWriter<File>,
    pub K_in_rope_file: BufWriter<File>,
    pub K_not_in_rope_file: BufWriter<File>,
}

enum ContrastType {
    Genotype,
    Strand,
}

fn build_tmp_id_vec(
        contrast_type: &ContrastType,
        geno_names: &Vec<String>,
        strand_names: &Vec<String>,
        param_name: &str,
        substring: &str,
) -> BoxedResult<Vec<ParamID>> {
    let mut tmp_ids: Vec<ParamID> = Vec::new();
    match contrast_type {
        // if we're doing genotype contrast, get the geno_idx
        ContrastType::Genotype => {
            // just for checking whether substring is in geno_names
            let geno_idx = geno_names.iter()
                .position(|s| *s == substring).ok_or("No such geno_name")?;
            for strand_name in strand_names {
                let param_id = ParamID{
                    parameter: param_name.to_string(),
                    geno: substring.to_string(),
                    strand: Some(strand_name.to_string()),
                };
                tmp_ids.push(param_id);
            }
        },
        // if we're doing strand contrast, get the strand_idx
        ContrastType::Strand => {
            // just for checking whether substring is in strand_names
            let strand_idx = strand_names.iter()
                .position(|s| *s == substring).ok_or("No such strand name")?;
            for geno_name in geno_names {
                let param_id = ParamID{
                    parameter: param_name.to_string(),
                    geno: geno_name.to_string(),
                    strand: Some(substring.to_string()),
                };
                tmp_ids.push(param_id);
            }
        },
    };
    Ok(tmp_ids)
}

fn append_to_vecs(
        substring: &str,
        minus: &bool,
        contrast_type: &ContrastType,
        strand_names: &Vec<String>,
        geno_names: &Vec<String>,
        numerator_vec: &mut Vec<Vec<ParamID>>,
        denominator_vec: &mut Vec<Vec<ParamID>>,
        contrast_str_vec: &mut Vec<String>,
        tmp_ids: &Vec<ParamID>,
) {
    if *minus {
        match contrast_type {
            ContrastType::Genotype => {
                for i in 0..strand_names.len() {
                    denominator_vec[i].push(tmp_ids[i].clone());
                    contrast_str_vec[i].push('-');
                    contrast_str_vec[i].push_str(substring);
                }
            },
            ContrastType::Strand => {
                for i in 0..geno_names.len() {
                    denominator_vec[i].push(tmp_ids[i].clone());
                    contrast_str_vec[i].push('-');
                    contrast_str_vec[i].push_str(substring);
                }
            }
        };
    } else {
        match contrast_type {
            ContrastType::Genotype => {
                for i in 0..strand_names.len() {
                    numerator_vec[i].push(tmp_ids[i].clone());
                    if numerator_vec[i].len() > 1 || denominator_vec[i].len() > 0 {
                        contrast_str_vec[i].push('+');
                    }
                    contrast_str_vec[i].push_str(substring);
                }
            },
            ContrastType::Strand => {
                for i in 0..geno_names.len() {
                    numerator_vec[i].push(tmp_ids[i].clone());
                    if numerator_vec[i].len() > 1 || denominator_vec[i].len() > 0 {
                        contrast_str_vec[i].push('+');
                    }
                    contrast_str_vec[i].push_str(substring);
                }
            }
        };
    }
}

impl Contrast {
    fn new(
            text: &str,
            parameter: Parameter,
            param_name: &str,
            expt_key: &ExptKey,
            contrast_type: ContrastType,
            out_direc: &str,
    ) -> BoxedResult<Vec<Contrast>> {

        // The trick here is to say, if we're doing genotype contrasts,
        // collect all the genotype information, then afterward, go through
        // all strand names to produce the genotype contrasts for each strand.
        //
        // However, if we're doing strand contrasts, collect the strand information,
        // then afterward go through all genotypes to produce strand contrasts
        // for each genotype.
        let strand_names = &expt_key.strand_map;
        let param_idx = expt_key.var_map.iter()
            .position(|s| *s == param_name).ok_or("No such param name")?;
        let geno_names = &expt_key.geno_map[param_idx];

        // If this is a genotype contrast, outer vec is for each strand.
        // If this is a strand contrast, outer vec is for each genotype.
        let mut numerator_vec: Vec<Vec<ParamID>> = Vec::new();
        let mut denominator_vec: Vec<Vec<ParamID>> = Vec::new();
        let mut contrast_str_vec: Vec<String> = Vec::new();
        // instantiate inner vecs
        match contrast_type {
            ContrastType::Genotype => {
                for i in 0..strand_names.len() {
                    denominator_vec.push(Vec::new());
                    numerator_vec.push(Vec::new());
                    contrast_str_vec.push(String::new());
                }
            },
            ContrastType::Strand => {
                for i in 0..geno_names.len() {
                    denominator_vec.push(Vec::new());
                    numerator_vec.push(Vec::new());
                    contrast_str_vec.push(String::new());
                }
            }
        };

        let mut substring = "".to_string();
        let mut minus = false;
        let mut first_char = true;
        for c in text.chars() {
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

                // Will store each contrast as a Contrast struct
                // need to know whether this is a genotype or a strand contrast so that we
                // can iterate over either genotypes (for strand contrast) or strands (for genotype
                // contrast) below to get all proper contrasts for this group of genotypes or
                // strands.
                let tmp_ids = build_tmp_id_vec(
                    &contrast_type,
                    &geno_names,
                    &strand_names,
                    param_name,
                    &substring,
                )?;

                append_to_vecs(
                    &substring,
                    &minus,
                    &contrast_type,
                    &strand_names,
                    &geno_names,
                    &mut numerator_vec,
                    &mut denominator_vec,
                    &mut contrast_str_vec,
                    &tmp_ids,
                );

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

        let tmp_ids = build_tmp_id_vec(
            &contrast_type,
            &geno_names,
            &strand_names,
            param_name,
            &substring,
        )?;
        append_to_vecs(
            &substring,
            &minus,
            &contrast_type,
            &strand_names,
            &geno_names,
            &mut numerator_vec,
            &mut denominator_vec,
            &mut contrast_str_vec,
            &tmp_ids,
        );

        let mut contrasts: Vec<Contrast> = Vec::new();
        // outer vecs are for either strand or genotype, depending on type of contrast
        match contrast_type {
            ContrastType::Genotype => {
                for i in 0..strand_names.len() {
                    let low_base = format!("{}_{}_lower.bedgraph", &parameter, &contrast_str_vec[i]);
                    let low_fname = Path::new(out_direc).join(low_base);
                    let low_f = File::create(low_fname)?;
                    let low_file = BufWriter::new(low_f);

                    let mean_base = format!("{}_{}_mean.bedgraph", &parameter, &contrast_str_vec[i]);
                    let mean_fname = Path::new(out_direc).join(mean_base);
                    let mean_f = File::create(mean_fname)?;
                    let mean_file = BufWriter::new(mean_f);

                    let up_base = format!("{}_{}_upper.bedgraph", &parameter, &contrast_str_vec[i]);
                    let up_fname = Path::new(out_direc).join(up_base);
                    let up_f = File::create(up_fname)?;
                    let up_file = BufWriter::new(up_f);

                    let Kge_base = format!("{}_{}_Kge.bedgraph", &parameter, &contrast_str_vec[i]);
                    let Kge_fname = Path::new(out_direc).join(Kge_base);
                    let Kge_f = File::create(Kge_fname)?;
                    let Kge_file = BufWriter::new(Kge_f);

                    let Klt_base = format!("{}_{}_Klt.bedgraph", &parameter, &contrast_str_vec[i]);
                    let Klt_fname = Path::new(out_direc).join(Klt_base);
                    let Klt_f = File::create(Klt_fname)?;
                    let Klt_file = BufWriter::new(Klt_f);

                    let K_in_rope_base = format!("{}_{}_K_in_rope.bedgraph", &parameter, &contrast_str_vec[i]);
                    let K_in_rope_fname = Path::new(out_direc).join(K_in_rope_base);
                    let K_in_rope_f = File::create(K_in_rope_fname)?;
                    let K_in_rope_file = BufWriter::new(K_in_rope_f);

                    let K_not_in_rope_base = format!("{}_{}_K_not_in_rope.bedgraph", &parameter, &contrast_str_vec[i]);
                    let K_not_in_rope_fname = Path::new(out_direc).join(K_not_in_rope_base);
                    let K_not_in_rope_f = File::create(K_not_in_rope_fname)?;
                    let K_not_in_rope_file = BufWriter::new(K_not_in_rope_f);

                    let cont = Contrast {
                        parameter: parameter.clone(),
                        text: text.to_string(),
                        contrast_str: contrast_str_vec[i].clone(),
                        numerators: numerator_vec[i].clone(),
                        denominators: denominator_vec[i].clone(),
                        low_file: low_file,
                        mean_file: mean_file,
                        up_file: up_file,
                        Kge_file: Kge_file,
                        Klt_file: Klt_file,
                        K_in_rope_file: K_in_rope_file,
                        K_not_in_rope_file: K_not_in_rope_file,
                    };
                    contrasts.push(cont);
                }
            },
            ContrastType::Strand => {
                for i in 0..geno_names.len() {
                    let low_base = format!("{}_{}_lower.bedgraph", &parameter, &contrast_str_vec[i]);
                    let low_fname = Path::new(out_direc).join(low_base);
                    let low_f = File::create(low_fname)?;
                    let low_file = BufWriter::new(low_f);
                    let mean_base = format!("{}_{}_mean.bedgraph", &parameter, &contrast_str_vec[i]);
                    let mean_fname = Path::new(out_direc).join(mean_base);
                    let mean_f = File::create(mean_fname)?;
                    let mean_file = BufWriter::new(mean_f);
                    let up_base = format!("{}_{}_upper.bedgraph", &parameter, &contrast_str_vec[i]);
                    let up_fname = Path::new(out_direc).join(up_base);
                    let up_f = File::create(up_fname)?;
                    let up_file = BufWriter::new(up_f);
                    let Kge_base = format!("{}_{}_Kge.bedgraph", &parameter, &contrast_str_vec[i]);
                    let Kge_fname = Path::new(out_direc).join(Kge_base);
                    let Kge_f = File::create(Kge_fname)?;
                    let Kge_file = BufWriter::new(Kge_f);

                    let Klt_base = format!("{}_{}_Klt.bedgraph", &parameter, &contrast_str_vec[i]);
                    let Klt_fname = Path::new(out_direc).join(Klt_base);
                    let Klt_f = File::create(Klt_fname)?;
                    let Klt_file = BufWriter::new(Klt_f);

                    let K_in_rope_base = format!("{}_{}_K_in_rope.bedgraph", &parameter, &contrast_str_vec[i]);
                    let K_in_rope_fname = Path::new(out_direc).join(K_in_rope_base);
                    let K_in_rope_f = File::create(K_in_rope_fname)?;
                    let K_in_rope_file = BufWriter::new(K_in_rope_f);

                    let K_not_in_rope_base = format!("{}_{}_K_not_in_rope.bedgraph", &parameter, &contrast_str_vec[i]);
                    let K_not_in_rope_fname = Path::new(out_direc).join(K_not_in_rope_base);
                    let K_not_in_rope_f = File::create(K_not_in_rope_fname)?;
                    let K_not_in_rope_file = BufWriter::new(K_not_in_rope_f);

                    let cont = Contrast {
                        parameter: parameter.clone(),
                        text: text.to_string(),
                        contrast_str: contrast_str_vec[i].clone(),
                        numerators: numerator_vec[i].clone(),
                        denominators: denominator_vec[i].clone(),
                        low_file: low_file,
                        mean_file: mean_file,
                        up_file: up_file,
                        Kge_file: Kge_file,
                        Klt_file: Klt_file,
                        K_in_rope_file: K_in_rope_file,
                        K_not_in_rope_file: K_not_in_rope_file,
                    };
                    contrasts.push(cont);
                }
            }
        };
        Ok(contrasts)
    }

    fn write_mean_to_file(
            &mut self,
            ctg: &str,
            start: &usize,
            end: &usize,
            val: f32,
    ) -> BoxedResult<()> {
        // Write the string content to the buffer
        let content = format!("{ctg}\t{start}\t{end}\t{val}\n");
        self.mean_file.write_all(content.as_bytes())?;
        Ok(())
    }

    fn flush_mean_buffer(&mut self) -> BoxedResult<()> {
        // Flush the buffer to ensure all data is written to the file
        self.mean_file.flush()?;
        Ok(())
    }

    fn write_lower_to_file(
            &mut self,
            ctg: &str,
            start: &usize,
            end: &usize,
            val: f32,
    ) -> BoxedResult<()> {
        // Write the string content to the buffer
        let content = format!("{ctg}\t{start}\t{end}\t{val}\n");
        self.low_file.write_all(content.as_bytes())?;
        Ok(())
    }

    fn flush_lower_buffer(&mut self) -> BoxedResult<()> {
        // Flush the buffer to ensure all data is written to the file
        self.low_file.flush()?;
        Ok(())
    }

    fn write_upper_to_file(
            &mut self,
            ctg: &str,
            start: &usize,
            end: &usize,
            val: f32,
    ) -> BoxedResult<()> {
        // Write the string content to the buffer
        let content = format!("{ctg}\t{start}\t{end}\t{val}\n");
        self.up_file.write_all(content.as_bytes())?;
        Ok(())
    }

    fn flush_upper_buffer(&mut self) -> BoxedResult<()> {
        // Flush the buffer to ensure all data is written to the file
        self.up_file.flush()?;
        Ok(())
    }

    fn write_Kge_to_file(
            &mut self,
            ctg: &str,
            start: &usize,
            end: &usize,
            val: f32,
    ) -> BoxedResult<()> {
        // Write the string content to the buffer
        let content = format!("{ctg}\t{start}\t{end}\t{val}\n");
        self.Kge_file.write_all(content.as_bytes())?;
        Ok(())
    }

    fn flush_Kge_buffer(&mut self) -> BoxedResult<()> {
        // Flush the buffer to ensure all data is written to the file
        self.Kge_file.flush()?;
        Ok(())
    }

    fn write_Klt_to_file(
            &mut self,
            ctg: &str,
            start: &usize,
            end: &usize,
            val: f32,
    ) -> BoxedResult<()> {
        // Write the string content to the buffer
        let content = format!("{ctg}\t{start}\t{end}\t{val}\n");
        self.Klt_file.write_all(content.as_bytes())?;
        Ok(())
    }

    fn flush_Klt_buffer(&mut self) -> BoxedResult<()> {
        // Flush the buffer to ensure all data is written to the file
        self.Klt_file.flush()?;
        Ok(())
    }

    fn write_K_in_to_file(
            &mut self,
            ctg: &str,
            start: &usize,
            end: &usize,
            val: f32,
    ) -> BoxedResult<()> {
        // Write the string content to the buffer
        let content = format!("{ctg}\t{start}\t{end}\t{val}\n");
        self.K_in_rope_file.write_all(content.as_bytes())?;
        Ok(())
    }

    fn flush_K_in_rope_buffer(&mut self) -> BoxedResult<()> {
        // Flush the buffer to ensure all data is written to the file
        self.K_in_rope_file.flush()?;
        Ok(())
    }

    fn write_K_not_in_to_file(
            &mut self,
            ctg: &str,
            start: &usize,
            end: &usize,
            val: f32,
    ) -> BoxedResult<()> {
        // Write the string content to the buffer
        let content = format!("{ctg}\t{start}\t{end}\t{val}\n");
        self.K_not_in_rope_file.write_all(content.as_bytes())?;
        Ok(())
    }

    fn flush_K_not_in_rope_buffer(&mut self) -> BoxedResult<()> {
        // Flush the buffer to ensure all data is written to the file
        self.K_not_in_rope_file.flush()?;
        Ok(())
    }

    pub fn flush_buffers(&mut self) -> BoxedResult<()> {
        self.flush_upper_buffer()?;
        self.flush_lower_buffer()?;
        self.flush_mean_buffer()?;
        self.flush_Kge_buffer()?;
        self.flush_Klt_buffer()?;
        self.flush_K_in_rope_buffer()?;
        self.flush_K_not_in_rope_buffer()?;
        Ok(())
    }

    pub fn write_vals_to_files(
            &mut self,
            ctg: &str,
            start: &usize,
            end: &usize,
            mean: f32,
            lower: f32,
            upper: f32,
            Kge: f32,
            Klt: f32,
            Kin: f32,
            Knotin: f32,
    ) -> BoxedResult<()> {
        self.write_mean_to_file(
            ctg,
            start,
            end,
            mean,
        )?;
        self.write_lower_to_file(
            ctg,
            start,
            end,
            lower,
        )?;
        self.write_upper_to_file(
            ctg,
            start,
            end,
            upper,
        )?;
        self.write_Kge_to_file(
            ctg,
            start,
            end,
            Kge,
        )?;
        self.write_Klt_to_file(
            ctg,
            start,
            end,
            Klt,
        )?;
        self.write_K_in_to_file(
            ctg,
            start,
            end,
            Kin,
        )?;
        self.write_K_not_in_to_file(
            ctg,
            start,
            end,
            Knotin,
        )?;
        Ok(())
    }

    pub fn to_string(&self) -> String {
        format!("{}_{}", &self.parameter, &self.contrast_str)
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
    pub var_map: Vec<String>,
    pub geno_map: Vec<Vec<String>>,
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

