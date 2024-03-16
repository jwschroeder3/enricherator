use clap::Parser;
use std::{fmt, thread, time};
use std::sync::mpsc;
use std::fs::File;
use std::collections::HashMap;
//use tokio::io::{AsyncReadExt};
use std::io::{self, Write, BufRead, BufReader, Read};
//use bytes::{BytesMut, BufMut};
use std::error::Error;
use std::str;

use ndarray::{Array5, Axis};
use ndarray::prelude::*;
//use ordered_float::OrderedFloat;
use serde::{Deserialize};
use serde_json;

//use bio_anno_rs::BEDGraphData;
//use bio_anno_rs::BEDGraphRecord;

const BUF_SIZE: usize = 1064;
const SAMP_NUM: usize = 501;

type BoxedResult<T> = std::result::Result<T, Box<dyn Error>>;


#[allow(dead_code)]
#[cfg(test)]
mod tests {
    use super::*;
    //use assert_approx_eq::*;

    fn set_up_tests() -> HashMap<String, &'static str> {
        let into_colname_text = "#header text \nlp__,log_p__,log_g__,sig_noise.1,sub_Alpha.1.1.1,sub_Beta.1.1.1,Alpha_Beta.1.1,Alpha_Beta.2.1,Alpha_Beta.3.1,Alpha_Beta.1.2,Alpha_Beta.2.2,Alpha_Beta.3.2\n";
        let into_sample_text = "#header text \nlp__,log_p__,log_g__,sig_noise.1,sub_Alpha.1.1.1,sub_Beta.1.1.1,Alpha_Beta.1.1,Alpha_Beta.2.1,Alpha_Beta.3.1,Alpha_Beta.1.2,Alpha_Beta.2.2,Alpha_Beta.3.2\n#next header \n0.0,0.00,0.0,0.99999,1.0,1.2,1.0,2.0,3.0,4.0,5.0,6.0\n";
        let first_batch = "#header text \nlp__,log_p__,log_g__,sig_noise.1,sub_Alpha.1.1.1,sub_Beta.1";
        let leftover = "sub_Beta.1";
        let the_rest = ".1.1,Alpha_Beta.1.1,Alpha_Beta.2.1,Alpha_Beta.3.1,Alpha_Beta.1.2,Alpha_Beta.2.2,Alpha_Beta.3.2\n";
        let colname_text = "lp__,log_p__,log_g__,sig_noise.1,sub_Alpha.1.1.1,sub_Beta.1.1.1,Alpha_Beta.1.1,Alpha_Beta.2.1,Alpha_Beta.3.1,Alpha_Beta.1.2,Alpha_Beta.2.2,Alpha_Beta.3.2\n";
        let fname = "test_data/test_parse.csv";
        let covar_fname = "test_data/covar_key.json";
        let contrast_covar_fname = "test_data/contrast_covar_key.json";

        let contrast_text = "#header text \nlp__,log_p__,log_g__,sig_noise.1,sub_Alpha.1.1.1,sub_Beta.1.1.1,Alpha_Beta.1.1,Alpha_Beta.2.1,Alpha_Beta.3.1,Alpha_Beta.1.2,Alpha_Beta.2.2,Alpha_Beta.3.2,Alpha_Beta.1.3,Alpha_Beta.2.3,Alpha_Beta.3.3,Alpha_Beta.1.4,Alpha_Beta.2.4,Alpha_Beta.3.4\n#next header \n0.0,0.00,0.0,0.99999,1.0,1.2,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0\n";

        let mut hash = HashMap::new();
        hash.insert(String::from("into_colname"), into_colname_text);
        hash.insert(String::from("into_sample"), into_sample_text);
        hash.insert(String::from("first_batch"), first_batch);
        hash.insert(String::from("leftover"), leftover);
        hash.insert(String::from("the_rest"), the_rest);
        hash.insert(String::from("colname_text"), colname_text);
        hash.insert(String::from("fname"), fname);
        hash.insert(String::from("covar_fname"), covar_fname);
        hash.insert(String::from("contrast_txt"), contrast_text);
        hash.insert(String::from("contrast_covar_fname"), contrast_covar_fname);
        hash
    }

    #[test]
    fn assign_to_samp_arr() {
        let idx = (0,0,0,0,0);

        let mut samp_arr = SamplesArray::zeros_by_shape(
            (2,2,2,2,2)
        );
        assert_eq!(samp_arr.samples[[0,0,0,0,0]], 0.0);
        samp_arr.assign_at(1.0, idx);
        assert_eq!(samp_arr.samples[[0,0,0,0,0]], 1.0);
        assert_eq!(samp_arr.samples[[0,0,0,0,1]], 0.0);
    }

    #[test]
    fn test_get_geno_contrast_keys() {
        let covar_key = CovarKey::read_covar_key("test_data/contrast_covar_key.json").unwrap();
        let vars = vec![String::from("Alpha"), String::from("Beta")];
        let geno_contrasts = Some(vec![String::from("mut-wt")]);
        let (
            idx_key,
            fname_key,
            geno_contrast_key_opt,
            geno_contrast_file_map_opt,
            _,
            _,
        ) = get_index_key(
            &covar_key,
            &vars,
            "test_data/test_outs",
            geno_contrasts,
            None,
        ).unwrap();
        let target_contrasts = vec![1,3];
        let result_contrasts_key = geno_contrast_key_opt.unwrap();
        assert_eq!(
            &target_contrasts[0],
            result_contrasts_key[1].get("numerator").unwrap()
        );
        assert_eq!(
            &target_contrasts[1],
            result_contrasts_key[1].get("denominator").unwrap()
        );
        let target_contrasts = vec![0,2];
        assert_eq!(
            &target_contrasts[0],
            result_contrasts_key[0].get("numerator").unwrap()
        );
        assert_eq!(
            &target_contrasts[1],
            result_contrasts_key[0].get("denominator").unwrap()
        );

    }

    #[test]
    fn test_get_index_key() {
        let covar_key = CovarKey::read_covar_key("test_data/covar_key.json").unwrap();
        let vars = vec![String::from("Alpha"),String::from("Beta")];
        let (idx_key,fname_key,contrast_key_opt,contrast_file_map_opt,_,_) = get_index_key(
            &covar_key,
            &vars,
            "test_data/test_outs",
            None,
            None,
        ).unwrap();
        assert_eq!(idx_key.get(&0).unwrap().0, 0);
        assert_eq!(idx_key.get(&0).unwrap().1, 0);
        assert_eq!(idx_key.get(&0).unwrap().2, 0);
        assert_eq!(idx_key.get(&1).unwrap().0, 0);
        assert_eq!(idx_key.get(&1).unwrap().1, 0);
        assert_eq!(idx_key.get(&1).unwrap().2, 1);
        assert_eq!(idx_key.get(&2).unwrap().0, 1);
        assert_eq!(idx_key.get(&2).unwrap().1, 0);
        assert_eq!(idx_key.get(&2).unwrap().2, 0);
        assert_eq!(idx_key.get(&3).unwrap().0, 1);
        assert_eq!(idx_key.get(&3).unwrap().1, 0);
        assert_eq!(idx_key.get(&3).unwrap().2, 1);
    }

    #[test]
    fn test_read_covar_key() {
        let hash = set_up_tests();
        let covar_key = CovarKey::read_covar_key("test_data/covar_key.json").unwrap();
        assert_eq!(covar_key.map[0].genotype_name, String::from("mut"));
        assert_eq!(covar_key.map[1].genotype_name, String::from("mut"));
        assert_eq!(covar_key.map[0].strand_name, String::from("both"));
        assert_eq!(covar_key.map[1].strand_name, String::from("both"));
        assert_eq!(covar_key.map[0].var_name, String::from("Alpha"));
        assert_eq!(covar_key.map[1].var_name, String::from("Beta"));
    }

    #[test]
    fn update_byter() {
        let new_bytes = vec![0,0,0,0,0];
        let bytes = Vec::new();
        let leftover = vec![5,5,5,5,5];
        let linetype = PriorByteIn::Header;
        let sample_idx = 0;
        let mut byter = BytesIterator{bytes, leftover, linetype, sample_idx};

        // test that byter.bytes is empty
        let target: Vec<u8> = Vec::new();
        assert_eq!(target, byter.bytes);
        // test that byter.leftover is what we expect
        let target: Vec<u8> = vec![5,5,5,5,5];
        assert_eq!(target, byter.leftover);

        // update byter
        byter.update(new_bytes);

        // test that byter.bytes has what we expect
        let target: Vec<u8> = vec![5,5,5,5,5,0,0,0,0,0];
        assert_eq!(target, byter.bytes);
        // test that byter.leftover is now empty
        let target: Vec<u8> = Vec::new();
        assert_eq!(target, byter.leftover);
    }
    
    #[test]
    fn parse_field_names() {
        let hash = set_up_tests();
        let text = hash.get("first_batch").unwrap();
        let target_leftover = hash.get("leftover").unwrap();
        let bytes = text.as_bytes();
        let mut byter = BytesIterator::new(bytes.to_vec());
        let mut col_idx = 0;
        //println!("{:?}", byter);
        assert_eq!(byter.bytes[0], b'#');
        let fields = byter.field_scan(&mut col_idx).unwrap();
        assert_eq!(target_leftover.as_bytes(), byter.leftover);
        //println!("{:?}", fields);
        assert_eq!("lp__", fields.fields[0].string);
        let text = hash.get("the_rest").unwrap();
        byter.update(text.as_bytes().to_vec());
        let fields = byter.field_scan(&mut col_idx).unwrap();
        //println!("{:?}", fields);
        assert_eq!("sub_Beta.1.1.1", fields.fields[0].string);
    }

    #[test]
    fn enter_colnames() {
        let hash = set_up_tests();
        let mut col_idx = 0;
        let bytes = hash.get("into_colname").unwrap().as_bytes();
        let mut byter = BytesIterator::new(bytes.to_vec());
        let fields = byter.field_scan(&mut col_idx).unwrap();
        assert_eq!("lp__", fields.fields[0].string);
        assert_eq!("Alpha_Beta.3.2", fields.fields[fields.len()-1].string);
    }

    #[test]
    fn filter_fields() {
        let hash = set_up_tests();
        let mut col_idx = 0;
        let bytes = hash.get("into_colname").unwrap().as_bytes();
        let mut byter = BytesIterator::new(bytes.to_vec());
        let mut fields = byter.field_scan(&mut col_idx).unwrap();
        let mut retained_fields: Vec<(usize, (usize, usize, usize, usize))> = Vec::new();
        let mut retained_idx: usize = 0;

        let covar_key = CovarKey::read_covar_key("test_data/covar_key.json").unwrap();
        let vars = vec![String::from("Alpha"),String::from("Beta")];
        let (idx_key,fname_key,_,_,_,_) = get_index_key(&covar_key, &vars, "test_data/test_outs", None, None).unwrap();

    // samples array is shape (geno, strand, params, position, samples)
        let mut samples_arr = SamplesArray::zeros_by_shape(
            (covar_key.geno_num, covar_key.strand_num, vars.len(), covar_key.pos_num, 2)
        );

        fields.parse_fields(
            &mut retained_fields,
            &mut retained_idx,
            &vars,
            &covar_key,
            &idx_key,
            &mut samples_arr,
        ).unwrap();

        let target_cols = vec![6,7,8,9,10,11];
        let target_genos = vec![0,0,0,0,0,0];
        let target_strands = vec![0,0,0,0,0,0];
        let target_vars = vec![0,0,0,1,1,1];
        for (i,col) in target_cols.iter().enumerate() {
            assert_eq!(col, &retained_fields[i].0);
            assert_eq!(target_genos[i], retained_fields[i].1.0);
            assert_eq!(target_strands[i], retained_fields[i].1.1);
            assert_eq!(target_vars[i], retained_fields[i].1.2);
        }
        // should not be equal, so assert that
        let target_cols = vec![7,8,9,10,11,12];
        for (i,col) in target_cols.iter().enumerate() {
            assert_ne!(col, &retained_fields[i].0);
        }

    }

    #[test]
    fn test_fields_into_samples() {
        let hash = set_up_tests();
        let mut col_idx = 0;
        let bytes = hash.get("into_sample").unwrap().as_bytes();
        let mut byter = BytesIterator::new(bytes.to_vec());
        let mut fields = byter.field_scan(&mut col_idx).unwrap();
        assert_eq!(fields.fields[fields.fields.len()-1].col_idx, 11);
        assert_eq!(fields.fields[fields.fields.len()-1].linetype, PriorByteIn::Samples(0));
    }

    #[test]
    fn test_into_samples() {
        let hash = set_up_tests();
        let bytes = hash.get("into_sample").unwrap().as_bytes();
        let mut col_idx = 0;
        let mut byter = BytesIterator::new(bytes.to_vec());

        let mut fields = byter.field_scan(&mut col_idx).unwrap();

        let mut retained_fields: Vec<(usize, (usize, usize, usize, usize))> = Vec::new();
        let mut retained_idx: usize = 0;

        let covar_key = CovarKey::read_covar_key("test_data/covar_key.json").unwrap();
        let vars = vec![String::from("Alpha"),String::from("Beta")];
        let (idx_key,fname_key,_,_,_,_) = get_index_key(&covar_key, &vars, "test_data/test_outs", None, None).unwrap();

        let mut samples_arr = SamplesArray::zeros_by_shape(
            (covar_key.geno_num, covar_key.strand_num, vars.len(), covar_key.pos_num, 2)
        );

        fields.parse_fields(
            &mut retained_fields,
            &mut retained_idx,
            &vars,
            &covar_key,
            &idx_key,
            &mut samples_arr,
        ).unwrap();

        assert_eq!(retained_fields[0].0, 6);
        assert_ne!(retained_fields[0].0, 7);
        assert_eq!(retained_fields[retained_fields.len()-1].0, 11);
        assert_eq!(samples_arr.samples[[0,0,0,0,0]], 1.0);
        assert_eq!(samples_arr.samples[[0,0,0,1,0]], 2.0);
        assert_eq!(samples_arr.samples[[0,0,0,2,0]], 3.0);
        assert_eq!(samples_arr.samples[[0,0,1,0,0]], 4.0);
        assert_eq!(samples_arr.samples[[0,0,1,1,0]], 5.0);
        assert_eq!(samples_arr.samples[[0,0,1,2,0]], 6.0);
    }
}

#[derive(Deserialize, Debug)]
struct CovarKey {
    map: Vec<CovarInfo>,
    resolution: u8,
    pos_num: usize,
    geno_num: usize,
    strand_num: usize,
}

#[derive(Deserialize, Debug)]
struct CovarInfo {
    var_name: String,
    genotype_name: String,
    strand_name: String,
}

impl CovarKey {
    fn read_covar_key(fname: &str) -> BoxedResult<CovarKey> {
        let mut file = File::open(fname)?;
        let mut data = String::new();
        file.read_to_string(&mut data)?;
        let covar_key: CovarKey = serde_json::from_str(&data)?;
        Ok(covar_key)
    }
}

/// Enumeration to handle information appropriately, depending on what we know about a line
#[derive(Debug,Copy,Clone,PartialEq)]
enum PriorByteIn {
    Header,
    EndOfHeader,
    HeaderAfterColnames,
    EndOfHeaderAfterColnames,
    ColNames,
    Samples(usize),
}

#[derive(Debug)]
enum HandleSamplesError {
    Placeholder,
}

impl fmt::Display for HandleSamplesError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            HandleSamplesError::Placeholder => write!(f, "nothing, really"),
        }
    }
}

impl Error for HandleSamplesError {}

#[derive(Debug)]
enum HandleFieldsError {
    CheckColnamesError,
    FilterFieldsError,
    ParseFloatError,
}

impl fmt::Display for HandleFieldsError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            HandleFieldsError::CheckColnamesError => write!(f, "nothing, really"),
            HandleFieldsError::FilterFieldsError => write!(f, "nothing, really"),
            HandleFieldsError::ParseFloatError => write!(f, "Unable to parse string to float"),
        }
    }
}

impl Error for HandleFieldsError {}

#[derive(Debug)]
enum ByteIterError {
    HashAfterSampleError,
    ReceiverError,
}

impl fmt::Display for ByteIterError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            ByteIterError::HashAfterSampleError =>
                write!(f, "Reached a '#' symbol after samples. This should not occur in cmdstan output files. Check your draws file."),
            ByteIterError::ReceiverError => write!(f, "Receiver error"),
        }
    }
}

impl Error for ByteIterError {}

#[derive(Debug)]
struct BytesIterator {
    bytes: Vec<u8>,
    leftover: Vec<u8>,
    linetype: PriorByteIn,
    sample_idx: usize,
}

#[derive(Debug)]
struct ByteField {
    bytes: Vec<u8>,
}

#[derive(Debug, PartialEq, Clone)]
struct StringField {
    string: String,
    linetype: PriorByteIn,
    col_idx: usize,
}

impl StringField {
    fn new(string: String, linetype: PriorByteIn, col_idx: usize) -> StringField {
        StringField { string, linetype, col_idx }
    }
}

impl ByteField {
    fn new() -> ByteField {
        let bytes = Vec::new();
        ByteField{bytes}
    }

    fn from_vec(bytes: Vec<u8>) -> ByteField {
        ByteField { bytes }
    }

    fn push(&mut self, byte: u8) {
        self.bytes.push(byte);
    }

    fn clear(&mut self) {
        self.bytes.clear();
    }
}

#[derive(Debug)]
struct SamplesArray {
    // samples array is shape (geno, strand, params, position, samples)
    // so that for each genotype, strand, parameter combination, there is a 2D sub-array
    // of shape (positions, samples). I think that placing samples last will promote
    // faster sorting if I'm not mistaken
    pub samples: Array5<f32>,
}

impl SamplesArray {

    fn zeros_by_shape(shape: (usize,usize,usize,usize,usize)) -> SamplesArray {
        let samples = Array5::zeros(shape);
        SamplesArray{ samples }
    }

    fn assign_at(&mut self, val: f32, idx: (usize,usize,usize,usize,usize)) {
        //println!("self.samples[idx]: {:?}", &self.samples[idx]);
        //println!("val: {}", &val);
        self.samples[idx] = val;
        //println!("self.samples[idx]: {:?}", &self.samples[idx]);
    }

    fn get_subarray_view(&self, idx: (usize,usize,usize)) -> BoxedResult<ArrayView<f32, Ix2>> {
        let subarr_view = self.samples.slice(s![idx.0, idx.1, idx.2, .., ..]);
        Ok(subarr_view)
    }

    /// Method to summarize samples for each genotype/strand/parameter/position.
    ///
    /// Args:
    /// -----
    /// subarr_map: just maps covariate indices to the first three indices of
    ///     self.samples
    /// writer_map: maps a covariate to its writers for easy grabbing of
    ///     correct writers for each genotype/strand/parameter
    /// position_map: for each index in the second-to-last axis of self.samples,
    ///     has "ctg\tstart\tend" ready and waiting to be written to files
    ///     in file_map
    /// threshold: defines value for gathering evidence ratios
    fn fetch_summaries(
            &mut self,
            subarr_map: &HashMap<usize, (usize, usize, usize)>,
            mut file_map: HashMap<usize, Vec<File>>,
            position_map: &Vec<String>,
            threshold: f32,
    ) -> BoxedResult<()> {

        let lower_idx = 24;
        let upper_idx = 474;

        // allocate the appropriate-size vec for reuse
        let mut vals: Vec<f32> = vec![0.0; SAMP_NUM];

        for (covar_idx, subarr_idx) in subarr_map.iter() {

            let covar_subarr_view = self.get_subarray_view(*subarr_idx)?;
            let files: &mut Vec<File> = file_map.get_mut(&covar_idx).unwrap();

            for (r_idx,row) in covar_subarr_view.axis_iter(Axis(0)).enumerate() {
                for i in 0..SAMP_NUM {
                    unsafe {
                        vals[i] = *row.uget(i);
                    }
                }
                vals.sort_unstable_by(|a,b| a.total_cmp(b));
                let lower = vals[lower_idx];
                let upper = vals[upper_idx];
                let median = (vals[249]+vals[250]) / 2.0;
                let mean = vals.iter().map(|x| *x).sum::<f32>() / 501.0;

                let mut negative_numbers: f32 = 0.0;
                for val in &vals {
                    if *val < 0.0 {
                        negative_numbers += 1.0;
                    } else {
                        break
                    }
                }
                let positive_numbers = 501.0 - negative_numbers;
                let mut K_gt = positive_numbers / negative_numbers;
                let mut K_lt = negative_numbers / positive_numbers;
                if K_gt > 501.0 {
                    K_gt = 501.0;
                }
                if K_lt > 501.0 {
                    K_lt = 501.0
                }
                let three_cols = &position_map[r_idx];
                files[0].write_all(
                    format!("{}\t{}\n", three_cols, lower).as_bytes()
                )?;
                files[1].write_all(
                    format!("{}\t{}\n", three_cols, upper).as_bytes()
                )?;
                files[2].write_all(
                    format!("{}\t{}\n", three_cols, median).as_bytes()
                )?;
                files[3].write_all(
                    format!("{}\t{}\n", three_cols, mean).as_bytes()
                )?;
                files[4].write_all(
                    format!("{}\t{}\n", three_cols, K_gt).as_bytes()
                )?;
                files[5].write_all(
                    format!("{}\t{}\n", three_cols, K_lt).as_bytes()
                )?;
            }
        }
        Ok(())
    }

    /// Method to summarize contrasts for each contrasts provided at cli.
    ///
    /// Args:
    /// -----
    /// subarr_map: just maps covariate indices to the first three indices of
    ///     self.samples
    /// contrast_map: A vec of 2-tuples, (numerator,denominator) pairings
    ///     that will be used to look up subarrays using subarr_map
    /// file_map: maps a covariate to its writers for easy grabbing of
    ///     correct writers for each contrast/strand/parameter
    /// position_map: for each index in the second-to-last axis of self.samples,
    ///     has "ctg\tstart\tend" ready and waiting to be written to files
    ///     in file_map
    /// threshold: defines value for gathering evidence ratios
    fn fetch_contrasts(
            &mut self,
            subarr_map: &HashMap<usize, (usize, usize, usize)>,
            genotype_contrast_map: Option<Vec<HashMap<String,usize>>>,
            mut genotype_file_map: Option<Vec<Vec<File>>>,
            strand_contrast_map: Option<Vec<HashMap<String,usize>>>,
            mut strand_file_map: Option<Vec<Vec<File>>>,
            position_map: &Vec<String>,
            genotype_threshold: f32,
            strand_threshold: f32,
    ) -> BoxedResult<()> {

        let lower_idx = 24;
        let upper_idx = 474;

        // allocate the appropriate-size vec for reuse
        let mut vals: Vec<f32> = vec![0.0; SAMP_NUM];

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// The below regions really should be wrapped into functions;
// they are nearly identical to each other
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
        if let Some(contrast_map) = strand_contrast_map {
            let mut file_map = strand_file_map.unwrap();
            for (cont_idx, cont_idx_map) in contrast_map.iter().enumerate() {

                let numer_key = cont_idx_map.get("numerator").unwrap();
                let denom_key = cont_idx_map.get("denominator").unwrap();

                let numer_idx = subarr_map.get(numer_key).unwrap();
                let denom_idx = subarr_map.get(denom_key).unwrap();

                let numer_subarr_view = self.get_subarray_view(*numer_idx)?;
                let denom_subarr_view = self.get_subarray_view(*denom_idx)?;
                //let contrast = numer_subarr_view - denom_subarr_view;
                let files: &mut Vec<File> = &mut file_map[cont_idx];

                for (r_idx,numer_row) in numer_subarr_view.axis_iter(Axis(0)).enumerate() {
                    let denom_row = denom_subarr_view.slice(s![r_idx,..]);
                    for i in 0..SAMP_NUM {
                        unsafe {
                            vals[i] = *numer_row.uget(i) - *denom_row.uget(i);
                        }
                    }
                    vals.sort_unstable_by(|a,b| a.total_cmp(b));
                    let lower = vals[lower_idx];
                    let upper = vals[upper_idx];
                    let median = (vals[249]+vals[250]) / 2.0;
                    let mean = vals.iter().map(|x| *x).sum::<f32>() / 501.0;

                    let mut negative_numbers: f32 = 0.0;
                    for val in &vals {
                        if *val < 0.0 {
                            negative_numbers += 1.0;
                        } else {
                            break
                        }
                    }
                    let positive_numbers = 501.0 - negative_numbers;
                    let mut K_gt = positive_numbers / negative_numbers;
                    let mut K_lt = negative_numbers / positive_numbers;
                    if K_gt > 501.0 {
                        K_gt = 501.0;
                    }
                    if K_lt > 501.0 {
                        K_lt = 501.0
                    }
                    let three_cols = &position_map[r_idx];
                    files[0].write_all(
                        format!("{}\t{}\n", three_cols, lower).as_bytes()
                    )?;
                    files[1].write_all(
                        format!("{}\t{}\n", three_cols, upper).as_bytes()
                    )?;
                    files[2].write_all(
                        format!("{}\t{}\n", three_cols, median).as_bytes()
                    )?;
                    files[3].write_all(
                        format!("{}\t{}\n", three_cols, mean).as_bytes()
                    )?;
                    files[4].write_all(
                        format!("{}\t{}\n", three_cols, K_gt).as_bytes()
                    )?;
                    files[5].write_all(
                        format!("{}\t{}\n", three_cols, K_lt).as_bytes()
                    )?;
                }
            }
        }

        if let Some(contrast_map) = genotype_contrast_map {
            let mut file_map = genotype_file_map.unwrap();
            for (cont_idx, cont_idx_map) in contrast_map.iter().enumerate() {

                let numer_key = cont_idx_map.get("numerator").unwrap();
                let denom_key = cont_idx_map.get("denominator").unwrap();

                let numer_idx = subarr_map.get(numer_key).unwrap();
                let denom_idx = subarr_map.get(denom_key).unwrap();

                let numer_subarr_view = self.get_subarray_view(*numer_idx)?;
                let denom_subarr_view = self.get_subarray_view(*denom_idx)?;
                //let contrast = numer_subarr_view - denom_subarr_view;
                let files: &mut Vec<File> = &mut file_map[cont_idx];

                for (r_idx,numer_row) in numer_subarr_view.axis_iter(Axis(0)).enumerate() {
                    let denom_row = denom_subarr_view.slice(s![r_idx,..]);
                    for i in 0..SAMP_NUM {
                        unsafe {
                            vals[i] = *numer_row.uget(i) - *denom_row.uget(i);
                        }
                    }
                    vals.sort_unstable_by(|a,b| a.total_cmp(b));
                    let lower = vals[lower_idx];
                    let upper = vals[upper_idx];
                    let median = (vals[249]+vals[250]) / 2.0;
                    let mean = vals.iter().map(|x| *x).sum::<f32>() / 501.0;

                    let mut negative_numbers: f32 = 0.0;
                    for val in &vals {
                        if *val < 0.0 {
                            negative_numbers += 1.0;
                        } else {
                            break
                        }
                    }
                    let positive_numbers = 501.0 - negative_numbers;
                    let mut K_gt = positive_numbers / negative_numbers;
                    let mut K_lt = negative_numbers / positive_numbers;
                    if K_gt > 501.0 {
                        K_gt = 501.0;
                    }
                    if K_lt > 501.0 {
                        K_lt = 501.0
                    }
                    let three_cols = &position_map[r_idx];
                    files[0].write_all(
                        format!("{}\t{}\n", three_cols, lower).as_bytes()
                    )?;
                    files[1].write_all(
                        format!("{}\t{}\n", three_cols, upper).as_bytes()
                    )?;
                    files[2].write_all(
                        format!("{}\t{}\n", three_cols, median).as_bytes()
                    )?;
                    files[3].write_all(
                        format!("{}\t{}\n", three_cols, mean).as_bytes()
                    )?;
                    files[4].write_all(
                        format!("{}\t{}\n", three_cols, K_gt).as_bytes()
                    )?;
                    files[5].write_all(
                        format!("{}\t{}\n", three_cols, K_lt).as_bytes()
                    )?;
                }
            }
        }
        Ok(())
    }

}

#[derive(Debug)]
struct FieldVec {
    fields: Vec<StringField>,
}

impl FieldVec {
    fn new() -> FieldVec {
        FieldVec{fields: Vec::<StringField>::new()}
    }

    unsafe fn push_bytes(&mut self, field: &mut ByteField, linetype: &PriorByteIn, col_idx: usize) {
        //println!("Input bytes: {:?}", &field);
        let field_string = unsafe {
            String::from_utf8_unchecked(field.bytes.to_vec())
        };
        //println!("String from the bytes: {:?}", &field_string);
        // if pushing to a fieldvec, always want to clear the field
        field.clear();
        //println!("Should be empty: {:?}", &field);
        let string_field = StringField::new(field_string, *linetype, col_idx);
        //println!("StringField: {:?}", &string_field);
        // push this field to fields
        self.fields.push(string_field);
    }

    fn len(&self) -> usize {
        self.fields.len()
    }

    // retained fields is used to map column indices to their proper position in the final array
    fn parse_fields(
            &mut self,
            retained_fields: &mut Vec<(usize, (usize, usize, usize, usize))>,
            retained_idx: &mut usize,
            vars: &Vec<String>,
            covar_key: &CovarKey,
            idx_key: &HashMap<usize, (usize, usize, usize)>,
            samples_arr: &mut SamplesArray,
    ) -> Result<(), HandleFieldsError> {
        //println!("Parsing fields passed from thread");
        for field in &self.fields {
            //println!("{:?}", field.linetype);
            match field.linetype {
                // if this field is a column name, do the following
                PriorByteIn::ColNames => {
                    // check whether the col name starts with "Alpha_Beta."
                    if field.string.starts_with("Alpha_Beta") {
                        // split string on ".", unpack to position and covar_idx
                        let field_split: Vec<&str> = field.string.split('.').collect();
                        // toss out the first element, keeping final two, converting to usize
                        let field_usize: Vec<usize> = field_split[1..3].iter().map(|x| {
                            x.parse::<usize>().unwrap()
                        }).collect();
                        let covar_idx = field_usize[1] - 1;
                        let position_idx = field_usize[0] - 1;
                        let this_covar_info = &covar_key.map[covar_idx];
                        // if any value in vars matches this var_name, do the following
                        if vars.iter().any(|x| x == &this_covar_info.var_name) {
                            // we're keeping values in this column, so we need to
                            // mark this col_idx as retained, and track which position in the
                            // samples array to place this column's values into once we get
                            // to its samples
                            let (geno_idx,strand_idx,var_idx) = idx_key.get(&covar_idx).unwrap();
                            retained_fields.push(
                                (field.col_idx, (*geno_idx,*strand_idx,*var_idx, position_idx))
                            );
                        }
                    }
                },
                PriorByteIn::Samples(_) => {
                    // set retained_idx to zero if required
                    if retained_idx == &retained_fields.len() {
                        *retained_idx = 0;
                    }
                    // now, if the column index == retained_fields[retained_idx.0, assign idx in
                    // samples_arr to val.
                    if field.col_idx == retained_fields[*retained_idx].0 {
                        let info = retained_fields[*retained_idx].1;

                        let val_result = field.string.parse();
                        let val: f32 = match val_result {
                            Ok(val) => val,
                            Err(_error) => return Err(HandleFieldsError::ParseFloatError)
                        };
                        
                        if let PriorByteIn::Samples(samp_idx) = field.linetype {
                            let idx = (info.0, info.1, info.2, info.3, samp_idx);
                            //println!("val: {}", &val);
                            //println!("idx: {:?}", &idx);
                            samples_arr.assign_at(
                                val,
                                idx,
                            )
                        }
                        *retained_idx += 1;
                    }
                },
                _ => continue
            }
        }
        Ok(())
    }
}

impl BytesIterator {
    fn new(bytes: Vec<u8>) -> BytesIterator {
        let leftover: Vec<u8> = Vec::new();
        let linetype = PriorByteIn::Header;
        let sample_idx: usize = 0;
        BytesIterator{bytes, leftover, linetype, sample_idx}
    }

    fn update(&mut self, mut bytes: Vec<u8>) {
        let mut complete_data = Vec::new();
        // if there is something leftover from prior buffer read, prepend to data here
        if !self.leftover.is_empty() {
            complete_data.append(&mut self.leftover);
            complete_data.append(&mut bytes);
            self.bytes = complete_data;
            self.leftover.clear();
        } else {
            self.bytes = bytes;
        }
        self.leftover.clear();
    }

    fn field_scan(&mut self, col_idx: &mut usize) -> Result<FieldVec, ByteIterError> {

        let mut fields = FieldVec::new();
        let mut field = ByteField::new();
        // look at every byte to scan for fields of interest
        for byte in &self.bytes {
            match *byte {
                // if byte is '#', entering a header
                b'#' => {
                    match self.linetype {
                        PriorByteIn::Header => {continue},
                        PriorByteIn::EndOfHeader => {
                            self.linetype = PriorByteIn::Header;
                        },
                        PriorByteIn::HeaderAfterColnames => {continue},
                        PriorByteIn::EndOfHeaderAfterColnames => {
                            self.linetype = PriorByteIn::HeaderAfterColnames;
                        },
                        PriorByteIn::ColNames => {
                            self.linetype = PriorByteIn::HeaderAfterColnames;
                            println!("Finished parsing column header line");
                        },
                        PriorByteIn::Samples(_) => {
                            return Err(ByteIterError::HashAfterSampleError)
                        },
                    };
                },
                // if the current byte is an EOL character
                b'\n' => {
                    match self.linetype {
                        // if the current byte is EOL and we're coming
                        // from a header, continue
                        PriorByteIn::Header => {
                            self.linetype = PriorByteIn::EndOfHeader;
                        },
                        PriorByteIn::EndOfHeader => {continue},
                        PriorByteIn::HeaderAfterColnames => {
                            self.linetype = PriorByteIn::EndOfHeaderAfterColnames;
                            *col_idx = 0;
                        },
                        PriorByteIn::EndOfHeaderAfterColnames => {continue},
                        // if the current byte is EOL and we're coming from colnames,
                        // gather the final field in the row
                        // and set linetype to headeraftercolnames
                        PriorByteIn::ColNames => {
                            unsafe { fields.push_bytes(&mut field, &self.linetype, *col_idx) };
                            *col_idx = 0;
                        },
                        PriorByteIn::Samples(samp_idx) => {
                            unsafe { fields.push_bytes(&mut field, &self.linetype, *col_idx) };
                            let this_samp_num = samp_idx + 1;
                            println!("Finished with sample {}", &this_samp_num);
                            // increment sample index by one and re-set col index to 0
                            self.linetype = PriorByteIn::Samples(this_samp_num);
                            *col_idx = 0;
                        },
                    };
                },
                b',' => {
                    match self.linetype {
                        // if the current byte is ',' and we're coming from a
                        // header, continue
                        PriorByteIn::Header => {continue},
                        PriorByteIn::EndOfHeader => {continue},
                        PriorByteIn::HeaderAfterColnames => {continue},
                        PriorByteIn::EndOfHeaderAfterColnames => {continue},
                        // if the current byte is ',' and we're in colnames,
                        // push field to fieldvec
                        PriorByteIn::ColNames => {
                            unsafe {
                                fields.push_bytes(&mut field, &self.linetype, *col_idx)
                            };
                            // increment col idx by one
                            *col_idx += 1;
                        },
                        // if current byte is ',' and we're in samples,
                        // push field to fieldvec
                        PriorByteIn::Samples(_) => {
                            unsafe {
                                fields.push_bytes(&mut field, &self.linetype, *col_idx)
                            };
                            // increment col idx by one
                            *col_idx += 1;
                        },
                    };
                },
                _ => {
                    match self.linetype {
                        PriorByteIn::Header => {continue},
                        PriorByteIn::EndOfHeader => {
                            self.linetype = PriorByteIn::ColNames;
                            println!("Beginning parsing column header line");
                            field.push(*byte);
                        },
                        PriorByteIn::HeaderAfterColnames => {continue},
                        PriorByteIn::EndOfHeaderAfterColnames => {
                            self.linetype = PriorByteIn::Samples(0);
                            field.push(*byte);
                        },
                        PriorByteIn::ColNames => {
                            field.push(*byte);
                        },
                        PriorByteIn::Samples(_) => {
                            field.push(*byte);
                        },
                    };
                },
            };
        };

        // place remaining bytes in field into leftover
        self.leftover = field.bytes;

        Ok(fields)
    }
}

/// returns the hashmap that will be used to map a covariate index to a position
/// in the first three axes of the final samples array and to files to be written
fn get_index_key(
        covar_key: &CovarKey,
        vars: &Vec<String>,
        out_direc: &str,
        genotype_contrasts: Option<Vec<String>>,
        strand_contrasts: Option<Vec<String>>,
) -> BoxedResult<(
        HashMap<usize, (usize, usize, usize)>,
        HashMap<usize, Vec<File>>,
        Option<Vec<HashMap<String, usize>>>,
        Option<Vec<Vec<File>>>,
        Option<Vec<HashMap<String, usize>>>,
        Option<Vec<Vec<File>>>,
)> {
    // samples array is shape (geno, strand, params, position, samples)
    let mut geno_vec: Vec<&str> = Vec::new();
    let mut strand_vec: Vec<&str> = Vec::new();

    let mut arr_map: HashMap<usize, (usize, usize, usize)> = HashMap::new();
    let mut file_map: HashMap<usize, Vec<File>> = HashMap::new();
    let summaries: Vec<&str> = vec![
        "lower", "upper", "median", "mean", "K_gt", "K_lt"
    ];

    for (i,covar_info) in covar_key.map.iter().enumerate() {
        let mut file_vec: Vec<File> = Vec::new();
        // if the genotype name isn't currently in geno_vec, append it
        if !geno_vec.iter().any(|x| x == &covar_info.genotype_name) {
            geno_vec.push(&covar_info.genotype_name);
        }
        // if the strand name isn't currently in strand_vec, append it
        if !strand_vec.iter().any(|x| x == &covar_info.strand_name) {
            strand_vec.push(&covar_info.strand_name);
        }
        // get the genotype index and strand index from the position of their name match in the
        // vecs
        let geno_idx = geno_vec.iter().position(|x| x == &covar_info.genotype_name).unwrap();
        let strand_idx = strand_vec.iter().position(|x| x == &covar_info.strand_name).unwrap();
        
        if vars.iter().any(|x| x == &covar_info.var_name) {
            // get the var index in the samples array that corresponds to this var_name
            let var_idx = vars.iter().position(|x| x == &covar_info.var_name).unwrap();
            arr_map.insert(i, (geno_idx, strand_idx, var_idx));
            for summary in &summaries {
                let fname = format!(
                    "{}/{}_{}_{}_{}.bedgraph",
                    out_direc,
                    &covar_info.genotype_name,
                    &covar_info.strand_name,
                    &covar_info.var_name,
                    summary,
                );
                let mut file = File::create(&fname)?;
                file_vec.push(file);
            }
            file_map.insert(i, file_vec);
        }
    }
    println!("geno_vec: {:?}", &geno_vec);

    let (geno_contrast_key,geno_contrast_file_map) = if let Some(contrasts) = genotype_contrasts {

        let mut geno_contrast_file_map: Vec<Vec<File>> = Vec::new();
        let mut geno_contrast_key: Vec<HashMap<String, usize>> = Vec::new();

        // for each contrast at cli, collect the two genotypes being
        // contrasted, get their indices in geno_vec
        // then iterate over values in arr_map above to store the keys
        // in arr_map that must be stored in geno_contrast_key to map
        // numerator and denominator to sub-arrays in final samples array
        for (i,contrast) in contrasts.iter().enumerate() {
            
            // get numerator and denominator genotype names
            let cont_vec: Vec<&str> = contrast.split("-").collect::<Vec<&str>>()
                .iter().map(|x| x.trim()).collect();
            let numer_geno = cont_vec[0];
            let denom_geno = cont_vec[1];

            // get indices for each of numerator and denominator
            let numer_geno_idx = geno_vec.iter()
                .position(|x| x == &numer_geno).unwrap();
            let denom_geno_idx = geno_vec.iter()
                .position(|x| x == &denom_geno).unwrap();

            for (j,strand_name) in strand_vec.iter().enumerate() {
                for (k,var_name) in vars.iter().enumerate() {

                    let mut this_contrast_key: HashMap<String, usize> = HashMap::new();
                    let mut file_vec: Vec<File> = Vec::new();

                    for summary in &summaries {
                        let fname = format!(
                            "{}/{}_{}_{}_{}.bedgraph",
                            out_direc,
                            contrast,
                            strand_name,
                            var_name,
                            summary,
                        );
                        let mut file = File::create(&fname)?;
                        file_vec.push(file);
                    }
                    geno_contrast_file_map.push(file_vec);

                    // now loop over the var key above to grab keys that match the strand, var,
                    // geno icx here
                    for (l, (query_geno_idx, query_strand_idx, query_var_idx)) in arr_map.iter() {
                        // if both strand and var indices match j and k, check if this is numerator
                        // or denominator and insert into the map
                        if query_strand_idx == &j && query_var_idx == &k {
                            if query_geno_idx == &numer_geno_idx {
                                this_contrast_key.insert(String::from("numerator"), *l);
                            } else if query_geno_idx == &denom_geno_idx {
                                this_contrast_key.insert(String::from("denominator"), *l);
                            }
                        }
                    }
                    geno_contrast_key.push(this_contrast_key);
                }
            }
        }
        (Some(geno_contrast_key),Some(geno_contrast_file_map))
    } else {
        (None,None)
    };

    let (strand_contrast_key,strand_contrast_file_map) = if let Some(contrasts) = strand_contrasts {

        let mut strand_contrast_file_map: Vec<Vec<File>> = Vec::new();
        let mut strand_contrast_key: Vec<HashMap<String, usize>> = Vec::new();

        for (i,contrast) in contrasts.iter().enumerate() {
            
            // get numerator and denominator strand names
            let cont_vec: Vec<&str> = contrast.split("-").collect::<Vec<&str>>()
                .iter().map(|x| x.trim()).collect();
            let numer_strand = cont_vec[0];
            let denom_strand = cont_vec[1];

            // get indices for each of numerator and denominator
            let numer_strand_idx = strand_vec.iter()
                .position(|x| x == &numer_strand).unwrap();
            let denom_strand_idx = strand_vec.iter()
                .position(|x| x == &denom_strand).unwrap();

            for (j,geno_name) in geno_vec.iter().enumerate() {
                for (k,var_name) in vars.iter().enumerate() {

                    let mut this_contrast_key: HashMap<String, usize> = HashMap::new();
                    let mut file_vec: Vec<File> = Vec::new();

                    for summary in &summaries {
                        let fname = format!(
                            "{}/{}_{}_{}_{}.bedgraph",
                            out_direc,
                            geno_name,
                            contrast,
                            var_name,
                            summary,
                        );
                        let mut file = File::create(&fname)?;
                        file_vec.push(file);
                    }
                    strand_contrast_file_map.push(file_vec);

                    // now loop over the var key above to grab keys that match the strand, var,
                    // geno idx here
                    for (l, (query_geno_idx, query_strand_idx, query_var_idx)) in arr_map.iter() {
                        // if both strand and var indices match j and k, check if this is numerator
                        // or denominator and insert into the map
                        if query_geno_idx == &j && query_var_idx == &k {
                            if query_strand_idx == &numer_strand_idx {
                                this_contrast_key.insert(String::from("numerator"), *l);
                            } else if query_strand_idx == &denom_strand_idx {
                                this_contrast_key.insert(String::from("denominator"), *l);
                            }
                        }
                    }
                    strand_contrast_key.push(this_contrast_key);
                }
            }
        }
        (Some(strand_contrast_key),Some(strand_contrast_file_map))
    } else {
        (None,None)
    };

    Ok((arr_map, file_map, geno_contrast_key, geno_contrast_file_map, strand_contrast_key, strand_contrast_file_map))
}

fn make_position_map(pos_file_name: &str) -> BoxedResult<Vec<String>> {
    let file = File::open(pos_file_name)?;
    let reader = BufReader::new(file);
    let mut line_map: Vec<String> = Vec::new();
    // Iterate over each line in the file.
    for line_result in reader.lines() {
        let line = line_result?;
        let stripped = line.trim_end();
        let split = stripped.rsplit_once('\t').unwrap();
        line_map.push(split.0.into());
    }
    Ok(line_map)
}

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Name of the person to greet
    #[arg(long)]
    draws_file: String,
    /// Comma-separated list of variables to keep
    #[arg(long, value_parser, use_value_delimiter=true, value_delimiter=',')]
    vars: Vec<String>,
    /// Comma-separated list of genotype-level contrasts to perform
    #[arg(long, value_parser, use_value_delimiter=true, value_delimiter=',')]
    genotype_contrasts: Option<String>,
    /// Comma-separated list of strand-level contrasts to perform
    #[arg(long, value_parser, use_value_delimiter=true, value_delimiter=',')]
    strand_contrasts: Option<String>,
    /// File containing covariate key
    #[arg(long)]
    covar_key_file: String,
    /// File containing positions to be mapped back to fit positions
    #[arg(long)]
    pos_file_name: String,
    /// Threshold for gathering evidence ratios for enrichment inferences
    #[arg(long, default_value_t = 1.0)]
    threshold: f32,
    /// Threshold for gathering evidence ratios for genotype contrasts
    #[arg(long, default_value_t = 0.0)]
    genotype_contrast_threshold: f32,
    /// Threshold for gathering evidence ratios for strand contrasts
    #[arg(long, default_value_t = 0.0)]
    strand_contrast_threshold: f32,
    /// Directory into which to write output files
    #[arg(long)]
    out_direc: String,
}

fn main() -> BoxedResult<()> {

    let args = Args::parse();
    let infname = &args.draws_file;
    let vars = args.vars;
    let covar_key_file = &args.covar_key_file;
    let pos_file_name = &args.pos_file_name;
    let threshold: f32 = args.threshold;
    let genotype_threshold: f32 = args.genotype_contrast_threshold;
    let strand_threshold: f32 = args.genotype_contrast_threshold;
    let out_direc: &str = &args.out_direc;
    let genotype_contrasts: Option<String> = args.genotype_contrasts;
    let strand_contrasts: Option<String> = args.strand_contrasts;

    let position_map = make_position_map(pos_file_name)?;
    let covar_key = CovarKey::read_covar_key(covar_key_file)?;

    let mut genotype_contrast_vec: Option<Vec<String>>  = Some(Vec::new());
    if let Some(cont) = genotype_contrasts {
        let split: Vec<String> = cont.split(',').map(|x| x.to_string()).collect();
        genotype_contrast_vec = Some(split);
    } else {
        genotype_contrast_vec = None;
    };

    let mut strand_contrast_vec: Option<Vec<String>>  = Some(Vec::new());
    if let Some(cont) = strand_contrasts {
        let split: Vec<String> = cont.split(',').map(|x| x.to_string()).collect();
        strand_contrast_vec = Some(split);
    } else {
        strand_contrast_vec = None;
    };

    let (
        idx_key,
        file_map,
        geno_contrast_key_opt,
        geno_contrast_file_map_opt,
        strand_contrast_key_opt,
        strand_contrast_file_map_opt,
    ) = get_index_key(
        &covar_key,
        &vars,
        out_direc,
        genotype_contrast_vec,
        strand_contrast_vec,
    )?;

    // channel for threads to communicate
    let (reader_tx, reader_rx) = mpsc::sync_channel::<Vec<u8>>(1);
    //let (reader_tx, reader_rx) = mpsc::channel::<Vec<u8>>();
    // If parser channel is sync_channel, hangs for a long time, maybe forever
    let (parser_tx, parser_rx) = mpsc::sync_channel(1);
    //let (parser_tx, parser_rx) = mpsc::channel();
    let recd: Vec<u8> = Vec::new();
    let mut byte_iterator = BytesIterator::new(recd);

    let file = File::open(infname)?;

    let mut reader = BufReader::new(file);
    let mut buffer = [0; BUF_SIZE];

    // spawn reader thread
    // this thread is simply responsible for reading BUF_SIZE bytes at a time
    // and passing those bytes to the next thread.
    let _reader_thread = thread::spawn(move || -> io::Result<()> {

        loop {
            let num_bytes_read = reader.read(&mut buffer)?;
            // zero bytes read indicates EOF
            if num_bytes_read == 0 {
                break
            }

            // send the slice of buffer as an owned vec with data to parser
            //println!("\nnum_bytes_read: {}", &num_bytes_read);
            let to_send = buffer[..num_bytes_read].to_vec();
            //println!("{:?}", &to_send);
            let send_result = reader_tx.send(to_send);
            //println!("send_result: {:?}", &send_result);
            send_result.unwrap();
        }

        Ok(())
    });

    // parser thread
    // this thread is repsonsible for splitting fields and determining what
    // type of row the fields are in, i.e., column header vs. samples
    // Information to pass to the handler thread is, the actual data in the field,
    // the type of data contained (col head vs samples), and
    // the sample number if we're in samples
    //
    // fields is a FieldVec, which has a single attribute, `fields`, of type Vec<StringField>
    let mut col_idx: usize = 0;
    let parser_thread = thread::spawn(move || -> Result<(), ByteIterError> {
        loop {
            let read_result = reader_rx.recv_timeout(time::Duration::from_secs(5));
            //let read_result = reader_rx.recv();
            let bytes_read: Vec<u8> = match read_result {
                Ok(bytes) => bytes,
                Err(error) => {
                    println!("Encountered a ReceiverError or hit timeout");
                    println!("{:?}", error);
                    return Err(ByteIterError::ReceiverError)
                }
            };
 
            //let bytes_read = read_result.unwrap();
            //println!("bytes read by reader_tx: {:?}", &bytes_read);
            byte_iterator.update(bytes_read);
            // fields is a FieldVec
            let fields = byte_iterator.field_scan(&mut col_idx)?;
            if !fields.fields.is_empty() {
                //println!("fields sent by parser tx: {:?}", &fields);
                //let copy = fields.fields.to_vec();
                let parser_send_result = parser_tx.send(fields);
                let _ = match parser_send_result {
                    Ok(send) => send,
                    Err(error) => {
                        println!("parser_send_error: {:?}", error);
                        break
                    }
                };
            }
        }
        Ok(())
    });

    //parser_thread.join()?;
    let mut samples_arr = SamplesArray::zeros_by_shape(
        (covar_key.geno_num, covar_key.strand_num, vars.len(), covar_key.pos_num, SAMP_NUM)
    );
    let mut retained_fields: Vec<(usize, (usize, usize, usize, usize))> = Vec::new();
    let mut retained_idx: usize = 0;

    // field handler
    // Responsible for accepting the FieldVec passed by parser thread
    // either assigning each field as retained or not, or assigning each value to
    // its proper position in the samples array.
    //
    // Tasks it must perform
    // 1) When in colname line, get vec of column indices and info
    //      Place the col index into retained_fields
    // 2) When in data, determine whether we're keeping each datum
    //      Place into appropriate cell in array.
    loop {
        let parse_result = parser_rx.recv_timeout(time::Duration::from_secs(5));
        let _ = match parse_result {
            Ok(mut fields) => {
                fields.parse_fields(
                    &mut retained_fields,
                    &mut retained_idx,
                    &vars,
                    &covar_key,
                    &idx_key,
                    &mut samples_arr,
                ).unwrap();
            },
            Err(error) => {
                println!("{:?}", error);
                break
            }
        };

    }

    println!("Finished parsing file, summarizing samples and writing output files.");
    samples_arr.fetch_summaries(
        &idx_key,
        file_map,
        &position_map,
        threshold,
    )?;

    samples_arr.fetch_contrasts(
        &idx_key,
        geno_contrast_key_opt,
        geno_contrast_file_map_opt,
        strand_contrast_key_opt,
        strand_contrast_file_map_opt,
        &position_map,
        genotype_threshold,
        strand_threshold,
    )?;

    Ok(())
}
