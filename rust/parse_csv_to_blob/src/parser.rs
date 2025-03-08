use clap::Parser;
use std::path::Path;
use std::fs::File;
use std::error::Error;
use std::collections::HashMap;
use std::fmt;

use crate::index::{Index, ParamID};
use crate::config::ExptKey;
use crate::blob::Blob;

const BUF_SIZE: usize = 1064;
const SAMP_NUM: usize = 501;

type BoxedResult<T> = std::result::Result<T, Box<dyn Error>>;

#[allow(dead_code)]
#[cfg(test)]
mod tests {
    use super::*;
    //use assert_approx_eq::*;

    fn set_up_tests() -> HashMap<String, &'static str> {
        let into_colname_text = "#header text \nlp__,log_p__,log_g__,sig_noise.1,sub_Alpha.1.1.1,sub_Beta.1.1.1,sub_Gamma.1.1.1,Alpha.1.1.1,Alpha.2.1.1,Alpha.3.1.1,Alpha.1.1.2,Alpha.2.1.2,Alpha.3.1.2,Gamma.3.2.2\n";
        let into_sample_text = "#header text \nlp__,log_p__,log_g__,sig_noise.1,sub_Alpha.1.1.1,sub_Beta.1.1.1,sub_Gamma.1.1.1,Alpha.1.1.1,Alpha.2.1.1,Alpha.3.1.1,Alpha.1.1.2,Alpha.2.1.2,Alpha.3.1.2,Gamma.3.2.2\n#next header \n0.0,0.00,0.0,0.99999,1.0,1.2,1.2,1.0,2.0,3.0,4.0,5.0,6.0,7.0\n";
        let first_batch = "#header text \nlp__,log_p__,log_g__,sig_noise.1,sub_Alpha.1.1.1,sub_Beta.1";
        let leftover = "sub_Beta.1";
        let the_rest = ".1.1,sub_Gamma.1.1.1,Alpha.1.1.1,Alpha.2.1.1,Alpha.3.1.1,Alpha.1.1.2,Alpha.2.1.2,Alpha.3.1.2\n";
        let colname_text = "lp__,log_p__,log_g__,sig_noise.1,sub_Alpha.1.1.1,sub_Beta.1.1.1,sub_Gamma.1.1.1,Alpha.1.1.1,Alpha.2.1.1,Alpha.3.1.1,Alpha.1.1.2,Alpha.2.1.2,Alpha.3.1.2\n";
        let fname = "test_data/test_parse.csv";
        let expt_key = "test_data/expt_key.json";
        let parse_key = "test_data/parse_key.toml";

        let contrast_text = "#header text \nlp__,log_p__,log_g__,sig_noise.1,sub_Alpha.1.1.1,sub_Beta.1.1.1,Alpha.1.1.1,Alpha.2.1.1,Alpha.3.1.1,Alpha.1.1.2,Alpha.2.1.2,Alpha.3.1.2,Beta.1.1.1,Beta.2.1.1,Beta.3.1.1,Beta.1.1.2,Beta.2.1.2,Beta.3.1.2,Gamma.1.1.1,Gamma.2.1.1,Gamma.3.1.1,Gamma.1.1.2,Gamma.2.1.2,Gamma.3.1.2\n#next header \n0.0,0.00,0.0,0.99999,1.0,1.2,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0\n0.0,0.00,0.0,0.99999,1.0,1.2,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5\n";

        let mut hash = HashMap::new();
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
        assert_eq!("Gamma.3.2.2", fields.fields[fields.len()-1].string);
    }

    #[test]
    fn var_names_contains() {
        let vars = [
            "Alpha",
            "Beta",
            "Gamma",
            "non_sub",
            "a_cubic",
            "b_cubic",
            "c_cubic",
        ];
        assert!(vars.contains(&"Beta"));
        assert!(!vars.contains(&"Betaa"));
    }

    #[test]
    fn filter_fields() {

        let hash = set_up_tests();
        let expt_key = ExptKey::read("test_data/expt_key.json").unwrap();
        let index = Index::new(&expt_key).unwrap();
        //index.to_file("test_data/out_direc/filter_fields").unwrap();

        let mut col_idx = 0;
        let bytes = hash.get("into_colname").unwrap().as_bytes();
        let mut byter = BytesIterator::new(bytes.to_vec());

        let mut fields = byter.field_scan(&mut col_idx).unwrap();
        let mut retained_fields: Vec<(usize, u64)> = Vec::new();
        let mut retained_idx: usize = 0;

        let mut blob_file = Blob::create(
            Path::new("test_data/out_direc/filter_fields").join("filter_fields.blob"),
            index.total_bytes,
            index.buffer_size,
        ).unwrap();

        fields.parse_fields(
            &index,
            &mut blob_file,
            &500,
            &mut retained_fields,
            &mut retained_idx,
        ).unwrap();

        let target_cols = vec![7,8,9,10,11,12,13];
        let result_cols: Vec<usize> = retained_fields.iter().map(|x| x.0).collect();
        assert_eq!(target_cols, result_cols);

        // multipliers are param, geno loops, strand loops, pos, sampnum, bytes per val
        let target_offsets = vec![
            0, // Alpha.1.1.1; wt, plus, 1
            200000, // Alpha.2.1.1; vpsR, plus, 1
            400000, // Alpha.3.1.1; vpsT, plus, 1
            500 * (2 - 1) * 4, // Alpha.1.1.2; wt, plus, 2
            200000 + 500 * (2 - 1) * 4, // Alpha.2.1.2; vpsR, plus, 2
            400000 + 500 * (2 - 1) * 4, // Alpha.3.1.2; vpsT, plus, 2
            1700000 + 500 * (2 - 1) * 4, // Gamma.3.2.2; vpsT, minus, 2
        ];
        let result_offsets: Vec<u64> = retained_fields.iter().map(|x| x.1).collect();
        assert_eq!(target_offsets, result_offsets);
    }

    #[test]
    fn test_fields_into_samples() {
        let hash = set_up_tests();
        let mut col_idx = 0;
        let bytes = hash.get("into_sample").unwrap().as_bytes();
        let mut byter = BytesIterator::new(bytes.to_vec());
        let mut fields = byter.field_scan(&mut col_idx).unwrap();
        assert_eq!(fields.fields[fields.fields.len()-1].col_idx, 13);
        assert_eq!(fields.fields[fields.fields.len()-1].linetype, PriorByteIn::Samples(0));
    }

    #[test]
    fn test_into_samples() {

        let hash = set_up_tests();
        let bytes = hash.get("into_sample").unwrap().as_bytes();

        let expt_key = ExptKey::read("test_data/expt_key.json").unwrap();
        let index = Index::new(&expt_key).unwrap();

        let mut col_idx = 0;
        let mut byter = BytesIterator::new(bytes.to_vec());

        let mut fields = byter.field_scan(&mut col_idx).unwrap();

        let mut retained_fields: Vec<(usize, u64)> = Vec::new();
        let mut retained_idx: usize = 0;

        let mut blob_file = Blob::create(
            Path::new("test_data/out_direc/filter_fields").join("filter_fields.blob"),
            index.total_bytes,
            index.buffer_size,
        ).unwrap();

        fields.parse_fields(
            &index,
            &mut blob_file,
            &500,
            &mut retained_fields,
            &mut retained_idx,
        ).unwrap();

        let target_offsets = vec![
            0, // Alpha.1.1.1; wt, plus, 1
            200000, // Alpha.2.1.1; vpsR, plus, 1
            400000, // Alpha.3.1.1; vpsT, plus, 1
            500 * (2 - 1) * 4, // Alpha.1.1.2; wt, plus, 2
            200000 + 500 * (2 - 1) * 4, // Alpha.2.1.2; vpsR, plus, 2
            400000 + 500 * (2 - 1) * 4, // Alpha.3.1.2; vpsT, plus, 2
            1700000 + 500 * (2 - 1) * 4, // Gamma.3.2.2; vpsT, minus, 2
        ];

        // The following f32 vals should have been written to blob at the positions in
        // target_offsets
        let values_in_blob = vec![1.0,2.0,3.0,4.0,5.0,6.0,7.0];

        for (i, (targ_off, targ_val)) in target_offsets.iter().zip(values_in_blob).enumerate() {
            // read_from_position reads to the Blob's `buffer` attribute
            let blob_bytes = blob_file.read_from_position(*targ_off);
            // convert first four bytes of buffer into [u8;4] to then cast to f32
            let mut val_bytes: [u8; 4] =  blob_file.buffer[0..4].try_into().unwrap();
            let blob_val = f32::from_be_bytes(val_bytes);
            // assert the float values are what we expect
            assert_eq!(targ_val, blob_val);
        }
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
    Samples(u64),
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
    NoMatchingParamIDError,
}

impl fmt::Display for HandleFieldsError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            HandleFieldsError::CheckColnamesError => write!(f, "nothing, really"),
            HandleFieldsError::FilterFieldsError => write!(f, "nothing, really"),
            HandleFieldsError::ParseFloatError => write!(f, "Unable to parse string to float"),
            HandleFieldsError::NoMatchingParamIDError => write!(f, "Searched for ParamID not found in index. Check your genotype and strand names carefully and try again."),
        }
    }
}

impl Error for HandleFieldsError {}

#[derive(Debug)]
pub enum ByteIterError {
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
pub struct BytesIterator {
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
pub struct StringField {
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
pub struct FieldVec {
    pub fields: Vec<StringField>,
}

//fn get_offset(
//        param_idx: usize,
//        param_dist: usize,
//        geno_idx: usize,
//        geno_dist: usize,
//        strand_idx: usize,
//        strand_dist: usize
//        pos_idx: usize,
//        pos_dist: usize,
//        samp_idx: usize,
//        samp_dist: usize,
//) -> usize {
//    let offset = (param_idx + 1) * param_dist
//        + (geno_idx + 1) * geno_dist
//        + (strand_idx + 1) * strand_dist
//        + (pos_idx + 1) * pos_dist
//        + (samp_idx + 1) * samp_dist
//        - 1;
//    offset
//}

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

    pub fn parse_fields(
            &mut self,
            index: &Index,
            blob: &mut Blob,
            sample_num: &u64,
            retained_fields: &mut Vec<(usize, u64)>,
            retained_idx: &mut usize,
    ) -> BoxedResult<()> {
        let vars = vec!["Alpha", "Beta", "Gamma", "non_sub"];
        //println!("Parsing fields passed from thread");
        for field in &self.fields {
            //println!("{:?}", field.linetype);
            match field.linetype {
                // if this field is a column name, do the following
                PriorByteIn::ColNames => {
                    // check whether the col name starts with "Alpha, Beta, Gamma, or
                    // non_sub"
                    if vars.iter().any(|var| field.string.starts_with(var)) {
                        // split string on ".", unpack to position and covar_idx
                        let field_split: Vec<&str> = field.string.split('.').collect();
                        // toss out the first element, keeping final two, converting to usize
                        let field_usize: Vec<usize> = field_split[1..4].iter().map(|x| {
                            x.parse::<usize>().unwrap()
                        }).collect();
                        let geno_idx = field_usize[0] - 1;
                        let strand_idx = field_usize[1] - 1;
                        let position_idx: u64 = (field_usize[2] - 1).try_into()?;
                        //let mut var_idx = 0;
                        let mut param_name = String::new();
                        for (var_i,var) in vars.iter().enumerate() {
                            if field.string.starts_with(var) {
                                //var_idx = var_i;
                                param_name = var.to_string();
                                break
                            }
                        };
                        // set the ParamID for this parameter
                        let param_id = ParamID {
                            parameter: param_name,
                            geno_idx: geno_idx,
                            strand_idx: Some(strand_idx),
                        };
                        println!("param_id {:?}", param_id);
                        let base_offset_opt = index.offsets.get(&param_id.to_string());
                        if let Some(base_offset) = base_offset_opt {
                            println!("base_offset: {:?}", &base_offset.offset);
                            println!("position_idx: {:?}", &position_idx);
                            println!("sample_num: {:?}", &sample_num);
                            // add position_idx * 4 bytes to the base offset for this param
                            let this_offset: u64 = base_offset.offset
                                + position_idx * sample_num * 4;
                            // we're keeping values in this column, so we need to
                            // mark this col_idx as retained, and track the offset at which
                            // we must write this field's values to blob
                            retained_fields.push(
                                (field.col_idx, this_offset)
                            );
                        } else {
                            return Err(Box::new(HandleFieldsError::NoMatchingParamIDError))
                        }
                    }
                },
                PriorByteIn::Samples(_) => {
                    // set retained_idx to zero if we reached the end of the the
                    // retained_fields vec in a prior field.
                    if retained_idx == &retained_fields.len() {
                        *retained_idx = 0;
                    }
                    // now, if the column index == retained_fields[retained_idx.0, 
                    // write the f32 val to blob at proper position
                    if field.col_idx == retained_fields[*retained_idx].0 {
                        let mut offset: u64 = retained_fields[*retained_idx].1;

                        let val_result = field.string.parse();
                        let val: f32 = match val_result {
                            Ok(val) => val,
                            Err(_error) => return Err(Box::new(HandleFieldsError::ParseFloatError))
                        };
                        
                        // get the sample_index from the linetype attribute
                        if let PriorByteIn::Samples(samp_idx) = field.linetype {
                            // add four bytes per sample
                            offset += samp_idx * 4;
                            blob.write_f32_to_position(offset, val)?;

                            //let idx = (info.0, info.1, info.2, info.3, samp_idx);
                            //println!("val: {}", &val);
                            //println!("idx: {:?}", &idx);
                            //samples_arr.assign_at(
                            //    val,
                            //    idx,
                            //)
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
    pub fn new(bytes: Vec<u8>) -> BytesIterator {
        let leftover: Vec<u8> = Vec::new();
        let linetype = PriorByteIn::Header;
        let sample_idx: usize = 0;
        BytesIterator{bytes, leftover, linetype, sample_idx}
    }

    pub fn update(&mut self, mut bytes: Vec<u8>) {
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

    pub fn field_scan(&mut self, col_idx: &mut usize) -> Result<FieldVec, ByteIterError> {

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

