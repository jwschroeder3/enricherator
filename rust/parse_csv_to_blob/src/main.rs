use clap::Parser;
use std::path::Path;
use std::fs::File;
use std::io::{self, Read, BufReader};
use std::time;
use std::error::Error;
use std::collections::HashMap;
use std::fmt;
use std::sync::mpsc;
use std::thread;

use parse_csv_to_blob::{
    config::{Args, ExptKey, ConfigFile},
    index::Index,
    blob::Blob,
    parser::{
        ByteIterError,
        BytesIterator,
    },
    BoxedResult,
};

const BUF_SIZE: usize = 1064;
const SAMP_NUM: usize = 501;

fn main() -> BoxedResult<()> {

    let args = Args::parse();
    let cfg = ConfigFile::read(&args.cfg_file)?;
    let expt_key = ExptKey::read(&cfg.input.expt_key_file)?;
    println!("Config: {:?}", &cfg);
    println!("expt_key: {:?}", &expt_key);

    let vars = [
        "Alpha",
        "Beta",
        "Gamma",
        "non_sub",
        "intercept",
        "slope",
    ];
    let enricherator_params = &vars[0..4];
    let sub_params = &vars[4..6];

    let infname = &cfg.input.draws_file;
    let out_direc: &str = &cfg.output.out_direc;

    let offset_idx = Index::new(&expt_key)?;
    offset_idx.to_file(out_direc)?;

    let blob_fname = Path::new(out_direc).join("params.blob");
    println!("Creating binary file {:?}", &blob_fname);
    // instantiate the Blob to manage IO to binary file
    let mut blob_file = Blob::create(
        blob_fname,
        offset_idx.total_bytes,
        offset_idx.buffer_size,
    )?;
    println!("Successfully created binary file");

    // channel for threads to communicate
    let (reader_tx, reader_rx) = mpsc::sync_channel::<Vec<u8>>(1);
    //let (reader_tx, reader_rx) = mpsc::channel::<Vec<u8>>();
    // If parser channel is sync_channel, hangs for a long time, maybe forever
    let (parser_tx, parser_rx) = mpsc::sync_channel(1);
    //let (parser_tx, parser_rx) = mpsc::channel();
    let recd: Vec<u8> = Vec::new();
    let mut byte_iterator = BytesIterator::new(recd);

    println!("Opening connection to draws file {:?}", &infname);
    let file = File::open(infname)?;
    println!("Successfully opened connection to draws file {:?}", &infname);

    let mut reader = BufReader::new(file);
    let mut read_buffer = [0; BUF_SIZE];

    // spawn reader thread
    // this thread is simply responsible for reading BUF_SIZE bytes at a time
    // and passing those bytes to the next thread.
    let _reader_thread = thread::spawn(move || -> io::Result<()> {

        loop {
            let num_bytes_read = reader.read(&mut read_buffer)?;
            // zero bytes read indicates EOF
            if num_bytes_read == 0 {
                break
            }

            // send the slice of buffer as an owned vec with data to parser
            //println!("\nnum_bytes_read: {}", &num_bytes_read);
            let to_send = read_buffer[..num_bytes_read].to_vec();
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

    let mut retained_fields: Vec<(usize, u64)> = Vec::new();
    let mut retained_idx: usize = 0;
    let sample_num: u64 = expt_key.sample_num[0].try_into()?;

    // Responsible for accepting the FieldVec passed by parser thread
    // either assigning each field as retained or not, or assigning each value to
    // its proper position in the samples array.
    //
    // Tasks it must perform
    // 1) When in colname line, get vec of column indices and offsets
    //      Place the col index into retained_fields
    // 2) When in data, determine whether we're keeping each datum
    //      Place bytes into appropriate cell in blob.
    loop {
        let parse_result = parser_rx.recv_timeout(time::Duration::from_secs(5));
        let _ = match parse_result {
            Ok(mut fields) => {
                fields.parse_fields(
                    &offset_idx,
                    &expt_key,
                    &mut blob_file,
                    &sample_num,
                    &mut retained_fields,
                    &mut retained_idx,
                ).unwrap();
            },
            Err(error) => {
                println!("{:?}", error);
                break
            }
        };

    }

    Ok(())
}

