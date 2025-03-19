use parse_csv_to_blob::{
    config::{Args, ExptKey, Config, PositionMap},
    index::{Index, ParamID},
    blob::Blob,
    BoxedResult,
};
use std::fs::File;
use std::io::{BufReader, BufRead, Seek, SeekFrom, Result};
use itertools::Itertools;

use clap::Parser;

use plotters::prelude::*;

fn cubic_interpolation(x: &f32, a: f32, b: f32, c: f32) -> f32 {
    a*x + b*x.powf(2.0) + c*x.powf(3.0)
}

fn plot_relationships(
        path: &str,
        geno_name: &str,
        non_sub: Vec<f32>,
        beta: Vec<f32>,
        gamma: Vec<f32>,
        a: f32,
        b: f32,
        c: f32,
) -> BoxedResult<()> {

    let y_min = non_sub.iter().cloned().reduce(f32::min).unwrap() - 0.1;
    let x_min = beta.iter().cloned().reduce(f32::min).unwrap() - 0.1;
    let y_max = non_sub.iter().cloned().reduce(f32::max).unwrap() + 0.1;
    let x_max = beta.iter().cloned().reduce(f32::max).unwrap() + 0.1;

    // get 100 values between x_min and x_max
    let step = x_max / (100 - 1) as f32;
    let line_x: Vec<f32> = (0..100).map(|i| 0.0 + step * i as f32).collect();
    let line_y: Vec<f32> = line_x.iter().map(|x| cubic_interpolation(x, a, b, c)).collect();
    let mut line: Vec<(f32, f32)> = vec![(x_min, 0.0), (0.0,0.0)];
    let mut curve: Vec<(f32, f32)> = line_x.iter().zip(line_y).map(|(x,y)| (*x as f32,y as f32)).collect();
    line.append(&mut curve);
    let points: Vec<(f32,f32)> = beta.iter().zip(non_sub).map(|(x,y)| (*x,y)).collect();
    let gamma_points: Vec<(f32,f32)> = beta.iter().zip(gamma).map(|(x,y)| (*x,y)).collect();

    let fname = format!("{}/{}_non_sub_vs_beta.png", path, geno_name);
    let root = BitMapBackend::new(&fname, (640, 480)).into_drawing_area();
    root.fill(&WHITE);
    let root = root.margin(10, 10, 10, 10);
    // After this point, we should be able to construct a chart context
    let mut chart = ChartBuilder::on(&root)
        // Set the caption of the chart
        .caption(format!("{} extract vs control relationship", geno_name), ("sans-serif", 40).into_font())
        // Set the size of the label region
        .x_label_area_size(20)
        .y_label_area_size(40)
        // Finally attach a coordinate on the drawing area and make a chart context
        .build_cartesian_2d(x_min..x_max, y_min..y_max)?;

    // Then we can draw a mesh
    chart
        .configure_mesh()
        // We can customize the maximum number of labels allowed for each axis
        .x_labels(5)
        .y_labels(5)
        // We can also change the format of the label text
        .y_label_formatter(&|y| format!("{:.1}", y))
        .draw()?;

    // We can draw point series
    chart.draw_series(PointSeries::of_element(
        points,
        5,
        &BLACK,
        &|c, s, st| {
            return EmptyElement::at(c)    // We want to construct a composed element on-the-fly
            + Circle::new((0,0),s,st.filled()); // At this point, the new pixel coordinate is established
            //+ Text::new(format!("{:?}", c), (10, 0), ("sans-serif", 10).into_font());
        },
    ))?;

    // And we can draw something in the drawing area
    chart.draw_series(LineSeries::new(line, &RED,))?;
    root.present()?;
    println!("Plot {fname} written to disk");

    let fname = format!("{}/{}_gamma_vs_beta.png", path, geno_name);
    let gamma_root = BitMapBackend::new(&fname, (640, 480)).into_drawing_area();
    gamma_root.fill(&WHITE);
    let gamma_root = gamma_root.margin(10, 10, 10, 10);
    // After this point, we should be able to construct a chart context
    let mut chart = ChartBuilder::on(&gamma_root)
        // Set the caption of the chart
        .caption(format!("{} subtracted (Gamma) vs control relationship", geno_name), ("sans-serif", 40).into_font())
        // Set the size of the label region
        .x_label_area_size(20)
        .y_label_area_size(40)
        // Finally attach a coordinate on the drawing area and make a chart context
        .build_cartesian_2d(x_min..x_max, y_min..y_max)?;

    // Then we can draw a mesh
    chart
        .configure_mesh()
        // We can customize the maximum number of labels allowed for each axis
        .x_labels(5)
        .y_labels(5)
        // We can also change the format of the label text
        .y_label_formatter(&|y| format!("{:.1}", y))
        .draw()?;

    // We can draw point series
    chart.draw_series(PointSeries::of_element(
        gamma_points,
        5,
        &BLACK,
        &|c, s, st| {
            return EmptyElement::at(c)    // We want to construct a composed element on-the-fly
            + Circle::new((0,0),s,st.filled()); // At this point, the new pixel coordinate is established
            //+ Text::new(format!("{:?}", c), (10, 0), ("sans-serif", 10).into_font());
        },
    ))?;

    gamma_root.present()?;
    println!("Plot {fname} written to disk");

    Ok(())
}

fn main() -> BoxedResult<()> {

    let args = Args::parse();
    println!("Reading config file at {}", &args.cfg_file);
    let mut cfg = Config::from_file(&args.cfg_file)?;
    //println!("Config: {:?}", &cfg);
    let expt_key = ExptKey::read(&cfg.input.expt_key_file)?;
    println!("Reading position map at {}", &cfg.input.pos_file_name);
    let pos_file = File::open(&cfg.input.pos_file_name)?;
    let mut pos_map_reader = BufReader::new(pos_file);
    let mut pos_line = String::new();
    let mut pos_map_line : Vec<&str> = Vec::new();

    // will need to gather some info for plotting the fitted cubic relationship
    // between non_sub and beta. To do this will require gathering beta and non_sub, though.
    // So, throw a big warning if they're not present in config parse section.
    let cubic_params = [
        "a_cubic",
        "b_cubic",
        "c_cubic",
    ];
    // these vecs will store the values for a given genotype
    let mut plot_beta: Vec<f32> = Vec::new();
    let mut plot_nonsub: Vec<f32> = Vec::new();
    let mut plot_gamma: Vec<f32> = Vec::new();

    // these vecs will store all values to later be iterated over for plotting
    // here we can use the indices used by stan to determine which nonsub/gamma go with which beta
    // Later, we'll sort each inner Vec<f32> by the values in the associated beta vec and
    // take a fraction of the values.
    let mut geno_plot_betas: Vec<Vec<f32>> = Vec::new();
    let mut geno_plot_nonsubs: Vec<Vec<f32>> = Vec::new();
    let mut geno_plot_gammas: Vec<Vec<f32>> = Vec::new();

    let out_direc: &str = &cfg.output.out_direc;
    println!("Reading blob index at {}/index.json", &out_direc);
    let index = Index::read(&out_direc)?;
    println!("Successfully read blob index from {}/index.json", &out_direc);

    // Get buffer size and position number to enable reading correct number of
    // bytes from blob, starting at correct position for every parameter.
    let blob_fname = format!("{out_direc}/params.blob");
    let buffer_size: usize = expt_key.sample_num[0] * 4;
    let pos_num: usize = expt_key.pos_num[0];

    // what indices from sorted samples to take the upper and lower 90% credible limits?
    let samp_num: u64 = expt_key.sample_num[0].try_into()?;
    let cl_low_idx: usize = ((samp_num / 20) - 1).try_into()?;
    let cl_up_idx: usize = (samp_num - (samp_num / 20)).try_into()?;

    // Open the Blob to manage reading from binary file
    let mut blob_file = Blob::open(
        &blob_fname,
        buffer_size,
    )?;
    println!("Successfully opened connection to blob");

    // Ideal scenario would be to:
    // 1. Get the index parameters in the order their data appear in `blob_file`. 
    // 2. Gather the distinct parameters present in the parsers and contrasts,
    //      and associate each parser and contrast with the appropriate parameters.
    // 3. For each param in ordered index, modify its parsers and contrasts appropriately.
    //
    // Below is not ideal, since it re-reads many paramters. Better would be the steps above.
    let mut geno_names = Vec::new();
    let param_thresh = &cfg.parse.threshold;
    println!("Parsing parameters of interest to save bedgraph summaries");
    for param_parser in &mut cfg.parse.param_ids {
        println!("Parsing parameter {:?}", &param_parser.param_id.to_string());
        let offset_idx = index.param_ids.iter()
            .position(|s| *s == param_parser.param_id.to_string())
            .ok_or("Could not parse param_id. No such param_id. Exiting now.")?;
        let base_offset: u64 = index.offsets[offset_idx].try_into()?;
        for pos_idx in 0..pos_num {
            // get the &("ctg", start, end) for this position
            let pos_line_bytes = pos_map_reader.read_line(&mut pos_line);
            //println!("pos_line: {:?}", &pos_line);
            pos_map_line = pos_line.trim().split(' ').collect();
            //println!("pos_map_line: {:?}", &pos_map_line);
            let ctg = pos_map_line[0];
            let start: usize = pos_map_line[1].parse::<usize>()?;
            let end: usize = pos_map_line[2].parse::<usize>()?;

            let pos: u64 = pos_idx.try_into()?;
            let this_offset: u64 = base_offset
                + pos * samp_num * 4;
            // reads the bytes into `blob_file.buffer`,
            // converts to f32 values, returns the Vec<f32>
            let mut samples = blob_file.read_f32_from_position(this_offset)?;
            // sort the samples to get conf lims
            samples.sort_unstable_by(|x,y| x.total_cmp(&y));
            let lower = samples[cl_low_idx];
            let upper = samples[cl_up_idx];
            let sum: f32 = samples.iter().sum();
            let mean: f32 = sum / samp_num as f32;

            // get number of samples less than threshold
            let lt_count = samples.partition_point(|x| x < param_thresh);
            let ge_count = samp_num - lt_count as u64;

            // Calculate the ratio
            let Klt = lt_count as f32 / ge_count as f32;
            let Kge = ge_count as f32 / lt_count as f32;

            if param_parser.param_id.parameter == "Beta" {
                plot_beta.push(mean.clone());
            } else if param_parser.param_id.parameter == "Gamma" {
                plot_gamma.push(mean.clone());
            } else if param_parser.param_id.parameter == "non_sub" {
                plot_nonsub.push(mean.clone());
            }

            param_parser.write_vals_to_files(
                &ctg,
                &start,
                &end,
                mean,
                lower,
                upper,
                Kge,
                Klt,
            )?;
            pos_line.clear();

        }
        // drain appropriate vec to appropriate vec of vecs
        if param_parser.param_id.parameter == "Beta" {
            geno_plot_betas.push(plot_beta.drain(..).collect());
        } else if param_parser.param_id.parameter == "Gamma" {
            geno_plot_gammas.push(plot_gamma.drain(..).collect());
            geno_names.push(param_parser.param_id.geno.clone());
        } else if param_parser.param_id.parameter == "non_sub" {
            geno_plot_nonsubs.push(plot_nonsub.drain(..).collect());
        }

        // return pos_map_reader to beginning of file
        pos_map_reader.seek(SeekFrom::Start(0))?;
        param_parser.flush_buffers();
    }

    //println!("cfg.geno_contrasts: {:?}", &cfg.geno_contrasts);
    if let Some(geno_contrasts) = &mut cfg.geno_contrasts {
        println!("Preparing genotype contrasts of interest to save bedgraph summaries");
        let geno_contrast_thresh = geno_contrasts.threshold;
        for contrast in &mut geno_contrasts.contrasts {
            println!("Preparing genotype contrast {}", contrast.to_string());

            let mut param_ids: Vec<ParamID> = Vec::new();
            let mut add: Vec<bool> = Vec::new();
            for param_id in &contrast.numerators {
                param_ids.push(param_id.clone());
                add.push(true);
            }
            for param_id in &contrast.denominators {
                param_ids.push(param_id.clone());
                add.push(false);
            }

            for pos_idx in 0..pos_num {

                // allocate vec for contrasts
                let mut contrast_vals = Vec::<f32>::with_capacity(samp_num as usize);
                for _ in 0..samp_num {
                    contrast_vals.push(0.0);
                }

                // get the &("ctg", start, end) for this position
                let pos_line_bytes = pos_map_reader.read_line(&mut pos_line);
                pos_map_line = pos_line.trim().split(' ').collect();
                let ctg = pos_map_line[0];
                let start: usize = pos_map_line[1].parse::<usize>()?;
                let end: usize = pos_map_line[2].parse::<usize>()?;

                let pos: u64 = pos_idx.try_into()?;
                for (param_id,numerator) in param_ids.iter().zip(add.clone()) {
                    let offset_idx = index.param_ids.iter()
                        .position(|s| *s == param_id.to_string())
                        .ok_or("No such param_id")?;
                    let base_offset: u64 = index.offsets[offset_idx].try_into()?;
                    let this_offset: u64 = base_offset
                        + pos * samp_num * 4;
                    // reads the bytes into `blob_file.buffer`,
                    // converts to f32 values, returns the Vec<f32>
                    let samples = blob_file.read_f32_from_position(this_offset)?;
                    for (i,val) in contrast_vals.iter_mut().enumerate() {
                        if numerator {
                            *val += samples[i];
                        } else {
                            *val -= samples[i];
                        }
                    }
                    // sort the samples to get conf lims
                    contrast_vals.sort_unstable_by(|x,y| x.total_cmp(&y));
                    let lower = contrast_vals[cl_low_idx];
                    let upper = contrast_vals[cl_up_idx];
                    let sum: f32 = contrast_vals.iter().sum();
                    let mean: f32 = sum / samp_num as f32;

                    // get number of samples less than threshold
                    let lt_count = contrast_vals.partition_point(|x| x < param_thresh);
                    let ge_count = samp_num - lt_count as u64;

                    // Calculate the ratio
                    let Klt = lt_count as f32 / ge_count as f32;
                    let Kge = ge_count as f32 / lt_count as f32;
                    contrast.write_vals_to_files(
                        &ctg,
                        &start,
                        &end,
                        mean,
                        lower,
                        upper,
                        Kge,
                        Klt,
                    )?;
                }
                pos_line.clear();

            }
            // return pos_map_reader to beginning of file
            pos_map_reader.seek(SeekFrom::Start(0))?;
            contrast.flush_buffers();
        }
    }    

    if let Some(strand_contrasts) = &mut cfg.strand_contrasts {
        println!("Preparing strand contrasts of interest to save bedgraph summaries");
        let strand_contrast_thresh = strand_contrasts.threshold;
        for contrast in &mut strand_contrasts.contrasts {
            println!("Preparing strand contrast {}", contrast.to_string());

            let mut param_ids: Vec<ParamID> = Vec::new();
            let mut add: Vec<bool> = Vec::new();
            for param_id in &contrast.numerators {
                param_ids.push(param_id.clone());
                add.push(true);
            }
            for param_id in &contrast.denominators {
                param_ids.push(param_id.clone());
                add.push(false);
            }

            for pos_idx in 0..pos_num {

                // allocate vec for contrasts
                let mut contrast_vals = Vec::<f32>::with_capacity(samp_num as usize);
                for _ in 0..samp_num {
                    contrast_vals.push(0.0);
                }

                // get the &("ctg", start, end) for this position
                let pos_line_bytes = pos_map_reader.read_line(&mut pos_line);
                pos_map_line = pos_line.trim().split(' ').collect();
                let ctg = pos_map_line[0];
                let start: usize = pos_map_line[1].parse::<usize>()?;
                let end: usize = pos_map_line[2].parse::<usize>()?;

                let pos: u64 = pos_idx.try_into()?;
                for (param_id,numerator) in param_ids.iter().zip(add.clone()) {
                    let offset_idx = index.param_ids.iter()
                        .position(|s| *s == param_id.to_string())
                        .ok_or("No such param_id")?;
                    let base_offset: u64 = index.offsets[offset_idx].try_into()?;
                    let this_offset: u64 = base_offset
                        + pos * samp_num * 4;
                    // reads the bytes into `blob_file.buffer`,
                    // converts to f32 values, returns the Vec<f32>
                    let samples = blob_file.read_f32_from_position(this_offset)?;
                    for (i,val) in contrast_vals.iter_mut().enumerate() {
                        if numerator {
                            *val += samples[i];
                        } else {
                            *val -= samples[i];
                        }
                    }
                    // sort the samples to get conf lims
                    contrast_vals.sort_unstable_by(|x,y| x.total_cmp(&y));
                    let lower = contrast_vals[cl_low_idx];
                    let upper = contrast_vals[cl_up_idx];
                    let sum: f32 = contrast_vals.iter().sum();
                    let mean: f32 = sum / samp_num as f32;

                    // get number of samples less than threshold
                    let lt_count = contrast_vals.partition_point(|x| x < param_thresh);
                    let ge_count = samp_num - lt_count as u64;

                    // Calculate the ratio
                    let Klt = lt_count as f32 / ge_count as f32;
                    let Kge = ge_count as f32 / lt_count as f32;
                    contrast.write_vals_to_files(
                        &ctg,
                        &start,
                        &end,
                        mean,
                        lower,
                        upper,
                        Kge,
                        Klt,
                    )?;
                }
                pos_line.clear();
            }
            // return pos_map_reader to beginning of file
            pos_map_reader.seek(SeekFrom::Start(0))?;
            contrast.flush_buffers();
        }
    }    

    // get the order of elements in each betas vec
    let mut beta_orders: Vec<Vec<usize>> = Vec::new();
    for plot_beta_vec in &geno_plot_betas {
        // Create a vector of tuples (index, value)
        let mut indexed_values: Vec<(usize, f32)> = plot_beta_vec.iter().cloned().enumerate().collect();
        // Sort by the value
        indexed_values.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        // Extract the indices into a separate vector
        let order: Vec<usize> = indexed_values.into_iter().map(|(index, _)| index).collect();
        beta_orders.push(order);
    }

    let mut this_beta_vec = Vec::new();
    let mut this_gamma_vec = Vec::new();
    let mut this_nonsub_vec = Vec::new();
    // extract every element for which the beta_order mod 100 == 0
    for (i, (nonsub_vec, gamma_vec)) in geno_plot_nonsubs.iter().zip(geno_plot_gammas).enumerate() {
        let beta_idx = cfg.beta_sub_idx[i] - 1;
        // borrow mutable reference to the inner vec of interest
        let beta_vec = &geno_plot_betas[beta_idx];
        let beta_order = &beta_orders[beta_idx];

        // keep every 100th ordered value
        for (j,val) in beta_order.iter().enumerate() {
            if val % 100 == 0 {
                this_beta_vec.push(beta_vec[j]);
                this_gamma_vec.push(gamma_vec[j]);
                this_nonsub_vec.push(nonsub_vec[j]);
            }
        }
        
        let geno_name = &geno_names[i];

        // now read the proper values from blob and get their means
        let param = format!("a_cubic.{geno_name}");
        println!("Parsing {param}");
        let offset_idx = index.param_ids.iter()
            .position(|s| *s == param)
            .ok_or("Could not parse param_id. No such param_id. Exiting now.")?;
        let this_offset: u64 = index.offsets[offset_idx].try_into()?;
        let mut samples = blob_file.read_f32_from_position(this_offset)?;
        let sum: f32 = samples.iter().sum();
        let a: f32 = sum / samp_num as f32;

        let param = format!("a_cubic.{geno_name}");
        println!("Parsing {param}");
        let offset_idx = index.param_ids.iter()
            .position(|s| *s == param)
            .ok_or("Could not parse param_id. No such param_id. Exiting now.")?;
        let this_offset: u64 = index.offsets[offset_idx].try_into()?;
        let mut samples = blob_file.read_f32_from_position(this_offset)?;
        let sum: f32 = samples.iter().sum();
        let b: f32 = sum / samp_num as f32;

        let param = format!("a_cubic.{geno_name}");
        println!("Parsing {param}");
        let offset_idx = index.param_ids.iter()
            .position(|s| *s == param)
            .ok_or("Could not parse param_id. No such param_id. Exiting now.")?;
        let this_offset: u64 = index.offsets[offset_idx].try_into()?;
        let mut samples = blob_file.read_f32_from_position(this_offset)?;
        let sum: f32 = samples.iter().sum();
        let c: f32 = sum / samp_num as f32;

        let _ = plot_relationships(
            &out_direc,
            &geno_name,
            this_nonsub_vec.clone(),
            this_beta_vec.clone(),
            this_gamma_vec.clone(),
            a,
            b,
            c,
        )?;
        this_beta_vec.clear();
        this_gamma_vec.clear();
        this_nonsub_vec.clear();
    }

    Ok(())
}
