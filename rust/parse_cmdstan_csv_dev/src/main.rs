use ordered_float::OrderedFloat;

use std::error::Error;
use std::iter::zip;

use tokio::fs::File;
use tokio::prelude::*;
use bytes::{BytesMut, BufMut};
use std::str;

async fn parse_lines<'a>(
        name: &str,
        nsamps: usize,
        vars: Vec<&'a str>,
        samp_out_fname: &str,
) -> Result<(Vec<(usize,(usize,usize,usize,String))>, Vec<Vec<OrderedFloat<f32>>>), Box<dyn Error>> {

    let mut file = File::open(name).await?;
    let mut buffer = BytesMut::with_capacity(1000);
    let mut line = String::new();
    let mut column_indices = Vec::new();
    let mut column_index = 0;
    let mut in_header = false;
    let mut in_samples = false;

    let target_word = "target word";  // put the sub-string you are looking for here

    loop {
        let bytes_read = file.read_buf(&mut buffer).await?;

        if bytes_read == 0 {
            break;
        }

        for byte in buffer.iter() {
            match *byte {
                b'#' => {
                    in_header = true;
                },
                b'\n' => {
                    in_header = false;
                    // check if line matches target word before a new line
                    if !line.is_empty() && line.trim() == target_word {
                        column_indices.push(column_index);
                    }
                    line.clear();
                    column_index = 0;
                },
                b',' if !in_header => {
                    if !in_samples {
                        // check if line matches target word
                        if line == target_word {
                            column_indices.push(column_index);
                        }
                        line.clear();
                    } else {
                    }
                    column_index += 1;
                },
                other if !in_header => {
                    line.push(other as char);
                },
                _ => continue,
            }
        }

        // reset the inner buffer to free memory space
        let remaining = buffer.split_to(bytes_read);
        buffer.clear();
        buffer.put(remaining);
    }
    
    println!("Target word found at columns: {:?}", column_indices);


}
 
#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {

    let args: Vec<String> = env::args().collect();
    let infname = &args[1];
    let summary_outfname = &args[2];
    let sample_outfname = &args[3];
    let samp_num: usize = args[4].parse().unwrap();
    let var_arg = &args[5]; // Alpha,Beta
    println!("Reading {:?} from {}", var_arg, infname);
    let vars: Vec<&str> = var_arg.split(",").collect();
 

    Ok(())
}
