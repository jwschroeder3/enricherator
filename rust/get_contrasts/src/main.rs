use std::env;
use tokio::io::{BufReader, BufWriter, AsyncWriteExt, AsyncBufReadExt, Lines};
use tokio::fs::File;
use std::error::Error;
use std::collections::HashMap;
//use std::io::prelude::*;
use ndarray::prelude::*;
use rayon::prelude::*;

use ordered_float::OrderedFloat;

#[cfg(test)]
mod tests {
    use super::*;
    use approx::*;

    #[test]
    fn test_parse_header() {
        let target_col_names: Vec<(usize,(usize,usize,usize,String))> = vec![
            (0, (0,1,49,String::from("Beta"))),
            (1, (0,0,50,String::from("Beta"))),
            (2, (0,1,50,String::from("Beta"))),
            //(3, (0,0,49,String::from("Alpha"))),
            //(4, (0,1,49,String::from("Alpha"))),
            //(5, (0,1,0,String::from("sub_Alpha"))),
            //(6, (0,1,0,String::from("sub_Beta"))),
        ];
        let col_name_line = String::from("Beta.1.2.50,Beta.1.1.51,Beta.1.2.51,Alpha.1.1.50,Alpha.1.2.50,sub_Alpha.1.2.1,sub_Beta.1.2.1");
        let split: Vec<&str> = col_name_line.split(',').collect();
        let col_names = get_col_names(split, String::from("Beta"));
        assert_eq!(target_col_names.len(), col_names.len());
        for (i,col) in col_names.iter().enumerate() {
            let target_col = &target_col_names[i];
            // assert indices match
            assert_eq!(col.0, target_col.0);
            // assert genotypes match
            assert_eq!(col.1.0, target_col.1.0);
            // assert strands match
            assert_eq!(col.1.1, target_col.1.1);
            // assert positions match
            assert_eq!(col.1.2, target_col.1.2);
            // assert vars match
            assert_eq!(col.1.3, target_col.1.3);
        }

        let target_col_names: Vec<(usize,(usize,usize,usize,String))> = vec![
            //(0, (0,1,49,String::from("Beta"))),
            //(1, (0,0,50,String::from("Beta"))),
            //(2, (0,1,50,String::from("Beta"))),
            (3, (0,0,49,String::from("Alpha"))),
            (4, (0,1,49,String::from("Alpha"))),
            //(5, (0,1,0,String::from("sub_Alpha"))),
            //(6, (0,1,0,String::from("sub_Beta"))),
        ];

        let col_name_line = String::from("Beta.1.2.50,Beta.1.1.51,Beta.1.2.51,Alpha.1.1.50,Alpha.1.2.50,sub_Alpha.1.2.1,sub_Beta.1.2.1");
        let split: Vec<&str> = col_name_line.split(',').collect();
        let col_names = get_col_names(split, String::from("Alpha"));
        assert_eq!(target_col_names.len(), col_names.len());
        for (i,col) in col_names.iter().enumerate() {
            let target_col = &target_col_names[i];
            // assert indices match
            assert_eq!(col.0, target_col.0);
            // assert genotypes match
            assert_eq!(col.1.0, target_col.1.0);
            // assert strands match
            assert_eq!(col.1.1, target_col.1.1);
            // assert positions match
            assert_eq!(col.1.2, target_col.1.2);
            // assert vars match
            assert_eq!(col.1.3, target_col.1.3);
        }
    }

    #[test]
    fn test_get_arr_shape() {
        let target_arr_shape = (2,2,2);
        // 12 columns, two beta genotypes
        let col_name_line = String::from("Beta.2.1.1,Beta.1.1.1,Beta.2.2.1,Beta.1.2.1,Beta.2.1.2,Beta.1.1.2,Beta.2.2.2,Beta.1.2.2,Alpha.1.1.1,Alpha.1.2.1,sub_Alpha.1.2.1,sub_Beta.1.2.1");
        let split: Vec<&str> = col_name_line.split(',').collect();
        let col_names = get_col_names(split, String::from("Beta"));
        let arr_shape = get_array_shape(&col_names);
        assert_eq!(target_arr_shape, arr_shape);
        let target_arr_shape = (2,2,1);
        assert_ne!(target_arr_shape, arr_shape);
    }

    fn add_to_slice(mut arr: Array4<f32>, a: usize, b: usize, c: usize, d: usize, val: f32) {
    }

    #[tokio::test]
    async fn test_fill_array() {
        let file = File::open("src/test_data/test.csv").await.unwrap();
        let reader = BufReader::new(file);
        let mut lines = reader.lines();
        let col_name_line = lines.next_line().await.unwrap().ok_or("unable to get next line").unwrap();
        let split: Vec<&str> = col_name_line.split(',').collect();
        let row_data = Vec::<OrderedFloat<f32>>::with_capacity(split.len());
        let col_names = get_col_names(split, String::from("Beta"));
        let arr_shape = get_array_shape(&col_names);


        let mut target_array = Array4::zeros((2,2,2,2));
        let mut sl = target_array.slice_mut(s![0,0,0,0]);
        sl += 2.0;
        let mut sl = target_array.slice_mut(s![0,0,0,1]);
        sl += 3.0;
        let mut sl = target_array.slice_mut(s![0,0,1,0]);
        sl += 4.0;
        let mut sl = target_array.slice_mut(s![0,0,1,1]);
        sl += 5.0;
        let mut sl = target_array.slice_mut(s![0,1,0,0]);
        sl += 2.5;
        let mut sl = target_array.slice_mut(s![0,1,0,1]);
        sl += 3.5;
        let mut sl = target_array.slice_mut(s![0,1,1,0]);
        sl += 2.0;
        let mut sl = target_array.slice_mut(s![0,1,1,1]);
        sl += 3.0;
        let mut sl = target_array.slice_mut(s![1,0,0,0]);
        sl += 1.0;
        let mut sl = target_array.slice_mut(s![1,0,0,1]);
        sl += 2.0;
        let mut sl = target_array.slice_mut(s![1,0,1,0]);
        sl += 3.0;
        let mut sl = target_array.slice_mut(s![1,0,1,1]);
        sl += 4.0;
        let mut sl = target_array.slice_mut(s![1,1,0,0]);
        sl += 1.5;
        let mut sl = target_array.slice_mut(s![1,1,0,1]);
        sl += 2.5;
        let mut sl = target_array.slice_mut(s![1,1,1,0]);
        sl += 1.0;
        let mut sl = target_array.slice_mut(s![1,1,1,1]);
        //sl += 3.0;
        sl += 2.0;
        println!("{:?}", target_array);
        for item in target_array.iter_mut() {
            *item = *item / (2.0_f32).ln();
        }
        println!("{:?}", target_array);
        let nsamps: usize = 2;
        let res = fill_array(arr_shape, nsamps, lines, &col_names, row_data).await.unwrap();
        println!("{:?}", res);
        
        for geno in 0..2 {
            for strand in 0..2 {
                for pos in 0..2 {
                    for samp in 0..2 {
                        //println!("{:?}", res.slice(s![geno,strand,pos,samp]));
                        //println!("{:?}", target_array.slice(s![geno,strand,pos,samp]));
                        assert_eq!(res.slice(s![geno,strand,pos,samp]), target_array.slice(s![geno,strand,pos,samp]));
                    }
                }
            }
        }
    }

    #[tokio::test]
    async fn test_genotype_contrast() {
        let file = File::open("src/test_data/test.csv").await.unwrap();
        let reader = BufReader::new(file);
        let mut lines = reader.lines();
        let col_name_line = lines.next_line().await.unwrap().ok_or("unable to get next line").unwrap();
        let split: Vec<&str> = col_name_line.split(',').collect();
        let row_data = Vec::<OrderedFloat<f32>>::with_capacity(split.len());
        let col_names = get_col_names(split, String::from("Beta"));
        let arr_shape = get_array_shape(&col_names);

        let nsamps: usize = 2;
        let data_arr = fill_array(arr_shape, nsamps, lines, &col_names, row_data).await.unwrap();
        let contrasts = vec!["genoA-genoB"];
        let mut target_array = Array4::zeros((
            contrasts.len(), data_arr.dim().1, data_arr.dim().2, data_arr.dim().3
        ));
        let mut sl = target_array.slice_mut(s![0,0,0,0]);
        sl += 1.0 / 2.0_f32.ln();
        let mut sl = target_array.slice_mut(s![0,0,0,1]);
        sl += 1.0 / 2.0_f32.ln();
        let mut sl = target_array.slice_mut(s![0,0,1,0]);
        sl += 1.0 / 2.0_f32.ln();
        let mut sl = target_array.slice_mut(s![0,0,1,1]);
        sl += 1.0 / 2.0_f32.ln();
        let mut sl = target_array.slice_mut(s![0,1,0,0]);
        sl += 1.0 / 2.0_f32.ln();
        let mut sl = target_array.slice_mut(s![0,1,0,1]);
        sl += 1.0 / 2.0_f32.ln();
        let mut sl = target_array.slice_mut(s![0,1,1,0]);
        sl += 1.0 / 2.0_f32.ln();
        let mut sl = target_array.slice_mut(s![0,1,1,1]);
        sl += 1.0 / 2.0_f32.ln();
        //sl += 11.0 / 2.0_f32.ln();

        let mut lut = HashMap::new();
        lut.insert(String::from("genoA"), 0);
        lut.insert(String::from("genoB"), 1);
        let result_array = do_contrasts(&contrasts, data_arr, lut, String::from("genotype")).unwrap();
        for cont in 0..1 {
            for strand in 0..2 {
                for pos in 0..2 {
                    for samp in 0..2 {
                        assert_abs_diff_eq!(result_array[(cont as usize, strand as usize, pos as usize, samp as usize)].into_inner(), target_array[(cont as usize, strand as usize, pos as usize, samp as usize)], epsilon=0.00001);
                    }
                }
            }
        }
    }

    #[tokio::test]
    async fn test_strand_contrast() {
        let file = File::open("src/test_data/test.csv").await.unwrap();
        let reader = BufReader::new(file);
        let mut lines = reader.lines();
        let col_name_line = lines.next_line().await.unwrap().ok_or("unable to get next line").unwrap();
        let split: Vec<&str> = col_name_line.split(',').collect();
        let row_data = Vec::<OrderedFloat<f32>>::with_capacity(split.len());
        let col_names = get_col_names(split, String::from("Beta"));
        let arr_shape = get_array_shape(&col_names);

        let nsamps: usize = 2;
        let data_arr = fill_array(arr_shape, nsamps, lines, &col_names, row_data).await.unwrap();
        let contrasts = vec!["plus-minus"];
        let mut target_array = Array4::<f32>::zeros((
            data_arr.dim().0, contrasts.len(), data_arr.dim().2, data_arr.dim().3
        ));

        // 2.0 - 2.5
        let mut sl = target_array.slice_mut(s![0,0,0,0]);
        sl += -0.5 / 2.0_f32.ln();
        // 3.0 - 3.5
        let mut sl = target_array.slice_mut(s![0,0,0,1]);
        sl += -0.5 / 2.0_f32.ln();
        // 4.0-2.0
        let mut sl = target_array.slice_mut(s![0,0,1,0]);
        sl += 2.0 / 2.0_f32.ln();
        // 5.0-3.0
        let mut sl = target_array.slice_mut(s![0,0,1,1]);
        sl += 2.0 / 2.0_f32.ln();
        // 1.0-1.5
        let mut sl = target_array.slice_mut(s![1,0,0,0]);
        sl += -0.5 / 2.0_f32.ln();
        // 2.0-2.5
        let mut sl = target_array.slice_mut(s![1,0,0,1]);
        sl += -0.5 / 2.0_f32.ln();
        // 3.0-1.0
        let mut sl = target_array.slice_mut(s![1,0,1,0]);
        sl += 2.0 / 2.0_f32.ln();
        // 4.0-2.0
        let mut sl = target_array.slice_mut(s![1,0,1,1]);
        sl += 2.0 / 2.0_f32.ln();

        let mut lut = HashMap::new();
        lut.insert(String::from("plus"), 0);
        lut.insert(String::from("minus"), 1);
        let result_array = do_contrasts(&contrasts, data_arr, lut, String::from("strand")).unwrap();
        for geno in 0..2 {
            for cont in 0..1 {
                for pos in 0..2 {
                    for samp in 0..2 {
                        let val1 = result_array[(geno as usize, cont as usize, pos as usize, samp as usize)].into_inner();
                        let val2 = target_array[(geno as usize, cont as usize, pos as usize, samp as usize)];
                        assert_abs_diff_eq!(val1, val2, epsilon=0.0001);
                    }
                }
            }
        }

    }
}

fn map_line_to_vec(
        line: String,
        arr: &mut Vec<OrderedFloat<f32>>,
) -> Result<(), Box<dyn Error>> {
    // mutate arr in place
    *arr = line.split(',')
        .map(|x| x.parse().unwrap())
        .collect::<Vec<OrderedFloat<f32>>>();
    Ok(())
}

async fn read_genotype_lut(name: &str) -> Result<HashMap<String, usize>, Box<dyn Error>> {
    let file = File::open(name).await?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    // get column names
    let _col_name_line = lines.next_line().await?.ok_or("unable to get next line")?;

    let mut geno_lut: HashMap<String, usize> = HashMap::new();

    //for line in lines {
    while let Some(line_str) = lines.next_line().await? {
        //let line_str = line?;
        let line_vec: Vec<&str> = line_str.split(',').collect();
        let geno_name = line_vec[0].to_string();
        let geno_num: usize = line_vec[1].parse()?;
        let geno_idx = geno_num - 1;
        geno_lut.insert(geno_name, geno_idx);
    }
    Ok(geno_lut)
}

async fn read_strand_lut(name: &str) -> Result<HashMap<String, usize>, Box<dyn Error>> {
    let file = File::open(name).await?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    // get column names
    let _col_name_line = lines.next_line().await?.ok_or("unable to get next line")?;

    let mut strand_lut: HashMap<String, usize> = HashMap::new();

    while let Some(line_str) = lines.next_line().await? {
    //for line in lines {
        //let line_str = line?;
        let line_vec: Vec<&str> = line_str.split(',').collect();
        let strand_name = line_vec[0].to_string();
        let strand_num: usize = line_vec[1].parse()?;
        let strand_idx = strand_num - 1;
        strand_lut.insert(strand_name, strand_idx);
    }
    Ok(strand_lut)
}

fn get_col_names(split: Vec<&str>, var: String) -> Vec<(usize,(usize,usize,usize,String))> {
    let col_names: Vec<(usize,(usize,usize,usize,String))> = split.par_iter()
        .enumerate()
        .filter(|(_,x)| x.starts_with(&var) )
        .map(|(i,x)| {
            let splitx: Vec<&str> = x.split('.').collect();
            let param = splitx[0];
            let geno: usize = splitx[1].parse::<usize>().unwrap()-1;
            let strand: usize = splitx[2].parse::<usize>().unwrap()-1;
            let pos: usize = splitx[3].parse::<usize>().unwrap()-1;
            (i, (geno,strand,pos,param.to_string()))
        }).collect();
    col_names
}

fn get_array_shape(col_names: &Vec<(usize,(usize,usize,usize,String))>) -> (usize,usize,usize) {

    let mut max_geno = 0;
    let mut max_strand = 0;
    let mut max_pos = 0;

    for col_info in col_names.iter() {
        let geno = col_info.1.0;
        let strand = col_info.1.1;
        let pos = col_info.1.2;
        // track this information to later set size of data array we'll fill with values
        if geno > max_geno {
            max_geno = geno
        }
        if strand > max_strand {
            max_strand = strand
        }
        if pos > max_pos {
            max_pos = pos
        }
    }
    (max_geno+1, max_strand+1, max_pos+1)
}

async fn fill_array(
        arr_shape: (usize,usize,usize),
        nsamps: usize,
        mut lines: Lines<BufReader<File>>,
        col_names: &Vec<(usize,(usize,usize,usize,String))>,
        mut row_data: Vec::<OrderedFloat<f32>>
) -> Result<Array<OrderedFloat<f32>,Ix4>, Box<dyn Error>> {
    // instantiate array of proper size
    let mut data_arr = Array4::zeros((arr_shape.0, arr_shape.1, arr_shape.2, nsamps));
    println!("data_arr dims: {:?}", data_arr.dim());

    let col_num = col_names.len();
    let mut i = 0;
    //for (i,line) in lines.enumerate() {
    while let Some(line_str) = lines.next_line().await? {
        if i+1 > nsamps {
            break
        }
        println!("Getting vals for sample {}", i+1);

        // row_data updated in place here
        map_line_to_vec(line_str, &mut row_data)?;
        for j in 0..col_num {
            let col_info = &col_names[j];
            let k = col_info.0;
            // convert from ln to log2, and place in vector of samples for this parameter
            let this_val = row_data[k] / (2.0_f32).ln();

            let (geno_idx,strand_idx,pos_idx,_) = col_info.1;
            let mut cell = data_arr.slice_mut(s![geno_idx, strand_idx, pos_idx, i]);
            cell.fill(this_val);
        }
        i += 1;
    }

    Ok(data_arr)
}

async fn parse_lines<'a>(
        name: &str,
        nsamps: usize,
        var: &str,
) -> Result<Array<OrderedFloat<f32>,Ix4>, Box<dyn Error>> {
    let mut var_name: String = var.into();
    let file_result = File::open(name).await;
    let file = match file_result {
        Ok(file) => file,
        Err(e) => return Err(Box::new(e)),
    };
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    // get column names
    let col_name_line = lines.next_line().await?.ok_or("unable to get next line")?;
    let split: Vec<&str> = col_name_line.split(',').collect();
    println!("total number of columns: {}", split.len());

    let row_data = Vec::<OrderedFloat<f32>>::with_capacity(split.len());

    var_name.push('.');
    let col_names = get_col_names(split, var_name);
    let arr_shape = get_array_shape(&col_names);
    let data_arr = fill_array(arr_shape, nsamps, lines, &col_names, row_data).await?;
    println!("Finished parsing input file.");
    Ok(data_arr)
}

fn do_contrasts(
        contrasts: &Vec<&str>,
        data: Array<OrderedFloat<f32>,Ix4>,
        lut: HashMap<String, usize>,
        contrast_type: String,
) -> Result<Array<OrderedFloat<f32>, Ix4>, Box<dyn Error>> {

    let mut result_arr = if contrast_type == String::from("genotype") {
        Array4::zeros((
            contrasts.len(), data.dim().1, data.dim().2, data.dim().3
        ))
    } else {
        Array4::zeros((
            data.dim().0, contrasts.len(), data.dim().2, data.dim().3
        ))
    };
    println!("Cont arr dims: {:?}", result_arr.dim());
    println!("Contrasts: {:?}", contrasts);

    for (i,contrast) in contrasts.iter().enumerate() {
        println!("Running contrast: {}", contrast);
        let cont_vec: Vec<&str> = contrast.split('-').collect::<Vec<&str>>()
            .iter().map(|x| x.trim()).collect();
        let numer_geno = cont_vec[0];
        let denom_geno = cont_vec[1];
        let numer_idx = lut.get(numer_geno).unwrap();
        let denom_idx = lut.get(denom_geno).unwrap();
        let numer_vals = if contrast_type == String::from("genotype") {
            data.slice(s![*numer_idx, .., .., ..]).to_owned()
        } else {
            data.slice(s![.., *numer_idx, .., ..]).to_owned()
        };
        let denom_vals = if contrast_type == String::from("genotype") {
            data.slice(s![*denom_idx, .., .., ..])
        } else {
            data.slice(s![.., *denom_idx, .., ..])
        };
        let mut contrast_slice = if contrast_type == String::from("genotype") {
            result_arr.slice_mut(s![i, .., .., ..])
        } else {
            result_arr.slice_mut(s![.., i, .., ..])
        };
        let contrast_vals = numer_vals - denom_vals;
        contrast_slice.assign(&contrast_vals);
    }
    println!("Finished calculating contrasts");
    Ok(result_arr)
}

async fn write_summaries(
        fname: &str,
        var_name: &str,
        cont_arr: Array<OrderedFloat<f32>, Ix4>,
        contrasts: Vec<&str>,
        cont_type: String,
) -> Result<(), Box<dyn Error>> {

    let file = File::create(fname).await?;
    let mut writer = BufWriter::new(file);

    println!("Summarizing contrasts and writing results to {}.", fname);
    writer.write(b"var_name\tgenotype\tstrand\tposition\tlower\tmedian\tupper\tmean\tK_gt\tK_lt\n").await?;

    let lower_idx = 24;
    //let lower_idx = 1;
    let upper_idx = 474;
    //let upper_idx = 8;
    let arr_dims = cont_arr.dim();

    // iterate over the contrasts
    if cont_type == String::from("genotype") {
        for i in 0..contrasts.len() {
            for strand_idx in 0..arr_dims.1 {
                for pos_idx in 0..arr_dims.2 {
                    writer.write(format!(
                        "{}\t{}\t{}\t{}\t",
                        var_name,i+1,strand_idx+1,pos_idx+1
                    ).as_bytes()).await?;
                    let mut samples: Vec<OrderedFloat<f32>> = cont_arr.slice(
                        s![i, strand_idx, pos_idx, ..]
                    ).to_vec();
                    samples.sort();
                    let lower = f32::from(samples[lower_idx]);
                    let median = f32::from(samples[249]+samples[250])/2.0;
                    let upper = f32::from(samples[upper_idx]);
                    let mean = f32::from(samples.iter().sum::<OrderedFloat<f32>>()
                        / samples.len() as f32);

                    let positive_numbers: f32 = samples.iter().filter(|&&x| x > OrderedFloat(0.0)).count() as f32;
                    let negative_numbers: f32 = samples.len() as f32 - positive_numbers;
                    let K_gt = positive_numbers / negative_numbers;
                    let K_lt = negative_numbers / positive_numbers;
                    writer.write(format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\n",
                        lower,median,upper,mean,K_gt,K_lt
                    ).as_bytes()).await?;
                }
            }
        }
    } else {
        for geno_idx in 0..arr_dims.0 {
            for strand_idx in 0..contrasts.len() {
                for pos_idx in 0..arr_dims.2 {
                    writer.write(format!(
                        "{}\t{}\t{}\t{}\t",
                        var_name,geno_idx+1,strand_idx+1,pos_idx+1
                    ).as_bytes()).await?;
                    let mut samples: Vec<OrderedFloat<f32>> = cont_arr.slice(
                        s![geno_idx, strand_idx, pos_idx, ..]
                    ).to_vec();
                    samples.sort();
                    let lower = f32::from(samples[lower_idx]);
                    let median = f32::from(samples[249]+samples[250])/2.0;
                    let upper = f32::from(samples[upper_idx]);
                    let mean = f32::from(samples.iter().sum::<OrderedFloat<f32>>()
                        / samples.len() as f32);
                    let positive_numbers: f32 = samples.iter().filter(|&&x| x > OrderedFloat(0.0)).count() as f32;
                    let negative_numbers: f32 = samples.len() as f32 - positive_numbers;
                    let K_gt = positive_numbers / negative_numbers;
                    let K_lt = negative_numbers / positive_numbers;
 
                    writer.write(format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\n",
                        lower,median,upper,mean,K_gt,K_lt
                    ).as_bytes()).await?;
                }
            }
        }
    }
    writer.flush().await?;
    Ok(())
}

#[tokio::main]
async fn main() {
    let args: Vec<String> = env::args().collect();
    let infname = &args[1];
    let summary_outfname = &args[2];
    let samp_num: usize = args[3].parse().unwrap();
    let var_name = &args[4];
    let contrast_arg = &args[5];
    let lut_fname = &args[6];
    let cont_type = &args[7];
    let contrasts: Vec<&str> = contrast_arg.split(",").collect();
    let geno_lut = read_genotype_lut(lut_fname).await.unwrap();
    println!("Geno LUT: {:?}: ", geno_lut);
    println!("Will perform the following contrasts for {}: {:?}", var_name, contrasts);
    println!("Reading {} from {}", var_name, infname);
    let samp_arr = parse_lines(infname, samp_num, var_name).await.unwrap();
    let results = do_contrasts(&contrasts, samp_arr, geno_lut, cont_type.to_string()).unwrap();
    write_summaries(
        summary_outfname,
        var_name,
        results,
        contrasts,
        cont_type.to_string(),
    ).await.unwrap();
}

