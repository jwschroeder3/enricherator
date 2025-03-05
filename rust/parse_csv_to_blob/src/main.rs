use clap::Parser;
use std::fs::File;
use std::error::Error;
use std::collections::HashMap;
use std::fmt;

mod config;
use config::{Args, ExptKey, ConfigFile};

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
        let into_sample_text = "#header text \nlp__,log_p__,log_g__,sig_noise.1,sub_Alpha.1.1.1,sub_Beta.1.1.1,sub_Gamma.1.1.1,Alpha.1.1.1,Alpha.2.1.1,Alpha.3.1.1,Alpha.1.1.2,Alpha.2.1.2,Alpha.3.1.2,Gamma.3.1.2\n#next header \n0.0,0.00,0.0,0.99999,1.0,1.2,1.2,1.0,2.0,3.0,4.0,5.0,6.0,7.0\n";
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

    //#[test]
    //fn filter_fields() {
    //    let hash = set_up_tests();
    //    let mut col_idx = 0;
    //    let bytes = hash.get("into_colname").unwrap().as_bytes();
    //    let mut byter = BytesIterator::new(bytes.to_vec());
    //    let mut fields = byter.field_scan(&mut col_idx).unwrap();
    //    let mut retained_fields: Vec<(usize, (usize, usize, usize, usize))> = Vec::new();
    //    let mut retained_idx: usize = 0;

    //    let expt_key = ExptKey::read("test_data/expt_key.json").unwrap();
    //    let vars = vec![String::from("Alpha"),String::from("Beta"),String::from("Gamma")];

    //    fields.parse_fields(
    //        &mut retained_fields,
    //        &mut retained_idx,
    //        &vars,
    //        &expt_key,
    //        &mut samples_arr,
    //    ).unwrap();
    //    let target_cols = vec![7,8,9,10,11,12,13];
    //    let target_genos = vec![0,1,2,0,1,2,2];
    //    let target_strands = vec![0,0,0,0,0,0,1];
    //    let target_vars = vec![0,0,0,0,0,0,2];
    //    for (i,col) in target_cols.iter().enumerate() {
    //        assert_eq!(col, &retained_fields[i].0);
    //        assert_eq!(target_genos[i], retained_fields[i].1.0);
    //        assert_eq!(target_strands[i], retained_fields[i].1.1);
    //        assert_eq!(target_vars[i], retained_fields[i].1.2);
    //    }
    //    // should not be equal, so assert that
    //    let target_cols = vec![6,7,8,9,10,11];
    //    for (i,col) in target_cols.iter().enumerate() {
    //        assert_ne!(col, &retained_fields[i].0);
    //    }

    //}

    //#[test]
    //fn test_fields_into_samples() {
    //    let hash = set_up_tests();
    //    let mut col_idx = 0;
    //    let bytes = hash.get("into_sample").unwrap().as_bytes();
    //    let mut byter = BytesIterator::new(bytes.to_vec());
    //    let mut fields = byter.field_scan(&mut col_idx).unwrap();
    //    assert_eq!(fields.fields[fields.fields.len()-1].col_idx, 13);
    //    assert_eq!(fields.fields[fields.fields.len()-1].linetype, PriorByteIn::Samples(0));
    //}

    //#[test]
    //fn test_into_samples() {
    //    let hash = set_up_tests();
    //    let bytes = hash.get("into_sample").unwrap().as_bytes();
    //    let mut col_idx = 0;
    //    let mut byter = BytesIterator::new(bytes.to_vec());

    //    let mut fields = byter.field_scan(&mut col_idx).unwrap();

    //    let mut retained_fields: Vec<(usize, (usize, usize, usize, usize))> = Vec::new();
    //    let mut retained_idx: usize = 0;

    //    let expt_key = ExptKey::read("test_data/expt_key.json").unwrap();
    //    let vars = vec![String::from("Alpha"),String::from("Beta"),String::from("Gamma")];
    //    let (idx_key,fname_key,_,_,_,_) = get_index_key(&expt_key, &vars, "test_data/test_outs", None, None).unwrap();

    //    let mut samples_arr = SamplesArray::zeros_by_shape(
    //        (expt_key.geno_num[0], expt_key.strand_num[0], vars.len(), expt_key.pos_num[0], 2)
    //    );

    //    fields.parse_fields(
    //        &mut retained_fields,
    //        &mut retained_idx,
    //        &vars,
    //        &expt_key,
    //        &mut samples_arr,
    //    ).unwrap();

    //    assert_eq!(retained_fields[0].0, 7);
    //    assert_eq!(retained_fields[retained_fields.len()-1].0, 13);
    //    assert_eq!(samples_arr.samples[[0,0,0,0,0]], 1.0);
    //    assert_eq!(samples_arr.samples[[1,0,0,0,0]], 2.0);
    //    assert_eq!(samples_arr.samples[[2,0,0,0,0]], 3.0);
    //    assert_eq!(samples_arr.samples[[0,0,0,1,0]], 4.0);
    //    assert_eq!(samples_arr.samples[[1,0,0,1,0]], 5.0);
    //    assert_eq!(samples_arr.samples[[2,0,0,1,0]], 6.0);
    //    assert_eq!(samples_arr.samples[[2,0,2,1,0]], 7.0);
    //}
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
            expt_key: &ExptKey,
    ) -> Result<(), HandleFieldsError> {
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
                        let position_idx = field_usize[2] - 1;
                        let mut var_idx = 0;
                        for (var_i,var) in vars.iter().enumerate() {
                            if field.string.starts_with(var) {
                                var_idx = var_i;
                                break
                            }
                        };
                        // we're keeping values in this column, so we need to
                        // mark this col_idx as retained, and track which position in the
                        // samples array to place this column's values into once we get
                        // to its samples
                        retained_fields.push(
                            (field.col_idx, (geno_idx,strand_idx,var_idx,position_idx))
                        );
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

fn main() -> BoxedResult<()> {

    let args = Args::parse();
    let expt_key = ExptKey::read(&args.expt_key)?;
    let cfg = ConfigFile::read(&args.cfg_file)?;
    println!("Config: {:?}", &cfg);
    println!("expt_key: {:?}", &expt_key);

    let vars = [
        "Alpha",
        "Beta",
        "Gamma",
        "non_sub",
        "a_cubic",
        "b_cubic",
        "c_cubic",
    ];
    let enricherator_params = &vars[0..4];
    let cubic_params = &vars[4..7];

    let infname = &cfg.input.draws_file;
    let out_direc: &str = &cfg.output.out_direc;

    ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////
    // At this point I need to: 
    // 1. figure out how many bytes my blob is,
    // 2. write zeros to that blob,
    // 3. write my index file with the offsets necessary to do the math to calculate pointer
    //     positions
    // 4. parse csv column names to assign each column index to its pointer position in blob
    ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////

    //let (
    //    idx_key,
    //    file_map,
    //    geno_contrast_key_opt,
    //    geno_contrast_file_map_opt,
    //    strand_contrast_key_opt,
    //    strand_contrast_file_map_opt,
    //) = get_index_key(
    //    &expt_key,
    //    &vars,
    //    out_direc,
    //    genotype_contrast_vec,
    //    strand_contrast_vec,
    //)?;

    //// channel for threads to communicate
    //let (reader_tx, reader_rx) = mpsc::sync_channel::<Vec<u8>>(1);
    ////let (reader_tx, reader_rx) = mpsc::channel::<Vec<u8>>();
    //// If parser channel is sync_channel, hangs for a long time, maybe forever
    //let (parser_tx, parser_rx) = mpsc::sync_channel(1);
    ////let (parser_tx, parser_rx) = mpsc::channel();
    //let recd: Vec<u8> = Vec::new();
    //let mut byte_iterator = BytesIterator::new(recd);

    //let file = File::open(infname)?;

    //let mut reader = BufReader::new(file);
    //let mut buffer = [0; BUF_SIZE];

    //// spawn reader thread
    //// this thread is simply responsible for reading BUF_SIZE bytes at a time
    //// and passing those bytes to the next thread.
    //let _reader_thread = thread::spawn(move || -> io::Result<()> {

    //    loop {
    //        let num_bytes_read = reader.read(&mut buffer)?;
    //        // zero bytes read indicates EOF
    //        if num_bytes_read == 0 {
    //            break
    //        }

    //        // send the slice of buffer as an owned vec with data to parser
    //        //println!("\nnum_bytes_read: {}", &num_bytes_read);
    //        let to_send = buffer[..num_bytes_read].to_vec();
    //        //println!("{:?}", &to_send);
    //        let send_result = reader_tx.send(to_send);
    //        //println!("send_result: {:?}", &send_result);
    //        send_result.unwrap();
    //    }

    //    Ok(())
    //});

    //// parser thread
    //// this thread is repsonsible for splitting fields and determining what
    //// type of row the fields are in, i.e., column header vs. samples
    //// Information to pass to the handler thread is, the actual data in the field,
    //// the type of data contained (col head vs samples), and
    //// the sample number if we're in samples
    ////
    //// fields is a FieldVec, which has a single attribute, `fields`, of type Vec<StringField>
    //let mut col_idx: usize = 0;
    //let parser_thread = thread::spawn(move || -> Result<(), ByteIterError> {
    //    loop {
    //        let read_result = reader_rx.recv_timeout(time::Duration::from_secs(5));
    //        //let read_result = reader_rx.recv();
    //        let bytes_read: Vec<u8> = match read_result {
    //            Ok(bytes) => bytes,
    //            Err(error) => {
    //                println!("Encountered a ReceiverError or hit timeout");
    //                println!("{:?}", error);
    //                return Err(ByteIterError::ReceiverError)
    //            }
    //        };
 
    //        //let bytes_read = read_result.unwrap();
    //        //println!("bytes read by reader_tx: {:?}", &bytes_read);
    //        byte_iterator.update(bytes_read);
    //        // fields is a FieldVec
    //        let fields = byte_iterator.field_scan(&mut col_idx)?;
    //        if !fields.fields.is_empty() {
    //            //println!("fields sent by parser tx: {:?}", &fields);
    //            //let copy = fields.fields.to_vec();
    //            let parser_send_result = parser_tx.send(fields);
    //            let _ = match parser_send_result {
    //                Ok(send) => send,
    //                Err(error) => {
    //                    println!("parser_send_error: {:?}", error);
    //                    break
    //                }
    //            };
    //        }
    //    }
    //    Ok(())
    //});

    //let mut retained_fields: Vec<(usize, (usize, usize, usize, usize))> = Vec::new();
    //let mut retained_idx: usize = 0;

    //// field handler
    //// Responsible for accepting the FieldVec passed by parser thread
    //// either assigning each field as retained or not, or assigning each value to
    //// its proper position in the samples array.
    ////
    //// Tasks it must perform
    //// 1) When in colname line, get vec of column indices and info
    ////      Place the col index into retained_fields
    //// 2) When in data, determine whether we're keeping each datum
    ////      Place into appropriate cell in array.
    //loop {
    //    let parse_result = parser_rx.recv_timeout(time::Duration::from_secs(5));
    //    let _ = match parse_result {
    //        Ok(mut fields) => {
    //            fields.parse_fields(
    //                &mut retained_fields,
    //                &mut retained_idx,
    //                &vars,
    //                &expt_key,
    //            ).unwrap();
    //        },
    //        Err(error) => {
    //            println!("{:?}", error);
    //            break
    //        }
    //    };

    //}

    Ok(())
}
