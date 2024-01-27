use std::env;
use tokio::io::{BufReader, BufWriter, AsyncWriteExt, AsyncBufReadExt, Lines};
use tokio::fs::File;
//use std::io::{BufReader, BufRead, BufWriter};
//use std::fs::File;
use std::error::Error;
use std::iter::zip;
//use std::io::prelude::*;
use rayon::prelude::*;

use ordered_float::OrderedFloat;

 
fn map_line_to_vec(
        line: String,
        arr: &mut Vec<OrderedFloat<f32>>,
) -> Result<(), Box<dyn Error>> {
    *arr = line.split(',')
        .map(|x| x.parse().unwrap())
        .collect::<Vec<OrderedFloat<f32>>>();
    Ok(())
}

async fn parse_lines<'a>(
        name: &str,
        nsamps: usize,
        vars: Vec<&'a str>,
        samp_out_fname: &str,
) -> Result<(Vec<(usize,(usize,usize,usize,String))>, Vec<Vec<OrderedFloat<f32>>>), Box<dyn Error>> {
    let file = File::open(name).await?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();
    // skip header
    let mut col_name_line = String::new();
    while let Some(line_str) = lines.next_line().await? {
    //let mut lines = lines.skip_while(|line| {
        //let line_str: &str = line.as_ref().unwrap();
        if line_str.starts_with("#") {
            continue
        } else {
            col_name_line = line_str;
            break
        }
    };

    let samp_outfile = File::create(samp_out_fname).await?;
    let mut samp_out_writer = BufWriter::new(samp_outfile);

    // get column names
    //let col_name_line = line_str;//lines.next().ok_or("unable to get next line")??;
    let split: Vec<&str> = col_name_line.split(',').collect();
    println!("total number of columns: {}", split.len());

    let mut row_data = Vec::<OrderedFloat<f32>>::with_capacity(split.len());

    let col_names: Vec<(usize,(usize,usize,usize,String))> = split.par_iter()
        .enumerate()
        .filter(|(_,x)| {
            vars.iter().any(|var_name| x.starts_with(var_name))
        }).map(|(i,x)| {
            let splitx: Vec<&str> = x.split('.').collect();
            let var_name = splitx[0];
            let geno: usize = splitx[1].parse::<usize>().unwrap()-1;
            let strand: usize = splitx[2].parse::<usize>().unwrap()-1;
            let pos: usize = splitx[3].parse::<usize>().unwrap()-1;

            (i, (geno,strand,pos,var_name.to_string()))
        }).collect();

    // write header to smaller samples file
    let col_num = col_names.len();
    for i in 0..col_num-1 {
        let col_name = &col_names[i];
        let var = format!(
            "{}.{}.{}.{},",
            col_name.1.3, col_name.1.0+1, col_name.1.1+1, col_name.1.2+1
        );
        samp_out_writer.write(var.as_bytes()).await?;
    }
    let last_var = format!(
        "{}.{}.{}.{}\n",
        col_names[col_num-1].1.3, col_names[col_num-1].1.0+1, col_names[col_num-1].1.1+1, col_names[col_num-1].1.2+1
    );
    samp_out_writer.write(last_var.as_bytes()).await?;


    // samples is a vector of vectors.
    // outer vec is number of parameters long
    // inner vecs are number of samples long
    // So, outer is columns, inners are rows for each column
    let mut samples: Vec<Vec<OrderedFloat<f32>>> = Vec::with_capacity(col_num);
    for _ in 0..col_num {
        samples.push(Vec::with_capacity(nsamps));
    }

    println!("First var column: {:?}", col_names[0]);
    println!("Final var column: {:?}", col_names[col_names.len()-1]);
    println!("Number of vars: {}", col_names.len());

    // skip the cmdstan output's header
    while let Some(line_str) = lines.next_line().await? {
    //let lines = lines.skip_while(|line| {
        //let line_str: &str = line.as_ref().unwrap();
        if line_str.starts_with("#") {
            continue
        } else {
            println!("Getting values for sample 1");
            map_line_to_vec(line_str, &mut row_data)?;
            for j in 0..(col_num-1) {
                let col_info = &col_names[j];
                let k = col_info.0;
                //let col_idx = col_info.1.0;
                let this_samp = &mut samples[j];
                // convert from ln to log2, and place in vector of samples for this parameter
                let this_val = row_data[k] / (2.0_f32).ln();
                this_samp.push(this_val);
                samp_out_writer.write(format!("{},",this_val).as_bytes()).await?;
            }

            // gather and write final number to line
            let col_info = &col_names[col_num-1];
            let k = col_info.0;
            //let col_idx = col_info.1.0;
            let this_samp = &mut samples[col_num-1];
            // convert from ln to log2, and place in vector of samples for this parameter
            let this_val = row_data[k] / (2.0_f32).ln();
            this_samp.push(this_val);
            samp_out_writer.write(format!("{}\n",this_val).as_bytes()).await?;
            break
        }
    };

    let mut i = 1;
    //for (i,line) in lines.enumerate() {
    while let Some(line_str) = lines.next_line().await? {
        if i+1 > nsamps {
            break
        }
        println!("Getting values for sample {}", i+1);

        //let line_str = line?;
        map_line_to_vec(line_str, &mut row_data)?;

        for j in 0..(col_num-1) {
            let col_info = &col_names[j];
            let k = col_info.0;
            //let col_idx = col_info.1.0;
            let this_samp = &mut samples[j];
            // convert from ln to log2, and place in vector of samples for this parameter
            let this_val = row_data[k] / (2.0_f32).ln();
            this_samp.push(this_val);
            samp_out_writer.write(format!("{},",this_val).as_bytes()).await?;
        }

        // gather and write final number to line
        let col_info = &col_names[col_num-1];
        let k = col_info.0;
        //let col_idx = col_info.1.0;
        let this_samp = &mut samples[col_num-1];
        // convert from ln to log2, and place in vector of samples for this parameter
        let this_val = row_data[k] / (2.0_f32).ln();
        this_samp.push(this_val);
        samp_out_writer.write(format!("{}\n",this_val).as_bytes()).await?;
        i += 1;
    }
    samp_out_writer.flush().await?;

    // sort the samples
    println!("Sorting samples for quantile summaries");
    for samp_vec in samples.iter_mut() {
        samp_vec.sort();
    }

    Ok((col_names, samples))
}


fn fetch_summaries(data: &Vec<Vec<OrderedFloat<f32>>>) -> Vec<Vec<f32>> {

    let lower_idx = 24;
    //let lower_idx = 1;
    let upper_idx = 474;
    //let upper_idx = 8;

    println!("Getting summary statistics for each parameter");
    let mut summaries: Vec<Vec<f32>> = Vec::with_capacity(data.len());
    for samples in data.iter() {
        let mut summary: Vec<f32> = Vec::with_capacity(4);
        summary.push(f32::from(samples[lower_idx]));
        summary.push(f32::from(samples[249]+samples[250])/2.0);
        //summary.push(f32::from((samples[4]+samples[5])/2.0));
        summary.push(f32::from(samples[upper_idx]));
        summary.push(
            f32::from(samples.iter().sum::<OrderedFloat<f32>>()
                / samples.len() as f32)
        );
        let positive_numbers: f32 = samples.iter().filter(|&&x| x > OrderedFloat(0.0)).count() as f32;
        let negative_numbers: f32 = samples.len() as f32 - positive_numbers;
        let mut K_gt = positive_numbers / negative_numbers;
        let mut K_lt = negative_numbers / positive_numbers;
        if K_gt > samples.len() as f32 {
            K_gt = samples.len() as f32;
        }
        if K_lt > samples.len() as f32 {
            K_lt = samples.len() as f32;
        }
        summary.push(K_gt);
        summary.push(K_lt);
        summaries.push(summary);
    }
    summaries
}

async fn write_summaries(
        fname: &str,
        summaries: Vec<Vec<f32>>,
        info: Vec<(usize,(usize,usize,usize,String))>,
) -> Result<(), Box<dyn Error>> {
    let file = File::create(fname).await?;
    let mut writer = BufWriter::new(file);

    println!("Writing summary statistics to file {}", fname);
    writer.write(b"var_name\tgenotype\tstrand\tposition\tlower\tmedian\tupper\tmean\tK_gt\tK_lt\n").await?;
    let iterator = zip(summaries, info);
    for dat in iterator {
        let line = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            dat.1.1.3, dat.1.1.0+1, dat.1.1.1+1, dat.1.1.2+1, dat.0[0], dat.0[1], dat.0[2], dat.0[3], dat.0[4], dat.0[5]
        );

        writer.write(line.as_bytes()).await?;
    }
    writer.flush().await?;
    Ok(())
}

#[tokio::main]
async fn main() {
    let args: Vec<String> = env::args().collect();
    let infname = &args[1];
    let summary_outfname = &args[2];
    let sample_outfname = &args[3];
    let samp_num: usize = args[4].parse().unwrap();
    let var_arg = &args[5]; // Alpha,Beta
    println!("Reading {:?} from {}", var_arg, infname);
    let vars: Vec<&str> = var_arg.split(",").collect();
    let (cols,samps) = parse_lines(infname, samp_num, vars, sample_outfname).await.unwrap();
    //println!("{:?}", &samps[0]);
    let summaries = fetch_summaries(&samps);
    //println!("{:?}", &summaries[0..2]);
    write_summaries(summary_outfname, summaries, cols).await.unwrap();
}
