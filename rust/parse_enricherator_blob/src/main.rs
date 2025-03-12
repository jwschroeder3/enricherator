use parse_csv_to_blob::{
    config::{Args, ExptKey, Config, PositionMap},
    index::Index,
    blob::Blob,
    BoxedResult,
};
use clap::Parser;

fn main() -> BoxedResult<()> {

    let args = Args::parse();
    let expt_key = ExptKey::read(&args.expt_key)?;
    let mut cfg = Config::from_file(&args.cfg_file, &expt_key)?;
    let position_map = PositionMap::read(
        &cfg.input.pos_file_name,
        &cfg.output.out_direc,
    )?;

    let cubic_params = [
        "a_cubic",
        "b_cubic",
        "c_cubic",
    ];

    let out_direc: &str = &cfg.output.out_direc;
    let index = Index::read(&out_direc)?;
    let blob_fname = format!("{out_direc}/params.blob");
    let buffer_size: usize = expt_key.sample_num[0] * 4;
    let pos_num: usize = expt_key.pos_num[0];
    let samp_num: u64 = expt_key.sample_num[0].try_into()?;
    let cl_low_idx: usize = ((samp_num / 20) - 1).try_into()?;
    let cl_up_idx: usize = (samp_num - (samp_num / 20)).try_into()?;

    // instantiate the Blob to manage IO to binary file
    let mut blob_file = Blob::open(
        &blob_fname,
        buffer_size,
    )?;

    // Ideal scenario would be to:
    // 1. Get the index parameters in the order their data appear in `blob_file`. 
    // 2. Gather the distinct parameters present in the parsers and contrasts,
    //      and associate each parser and contrast with the appropriate parameters.
    // 3. For each param in ordered index, modify its parsers and contrasts appropriately.
    //
    // I need to get the position map in here. Could just be ctg_id, start
    //
    let param_thresh = &cfg.parse.threshold;
    for param_parser in &mut cfg.parse.param_ids {
        let offset_idx = index.param_ids.iter()
            .position(|s| *s == param_parser.param_id.to_string())
            .ok_or("No such param_id")?;
        let base_offset: u64 = index.offsets[offset_idx].try_into()?;
        for pos_idx in 0..pos_num {
            // get the &("ctg", start, end) for this position
            let (ctg,start,end) = &position_map.map[pos_idx];
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

            param_parser.write_vals_to_files(
                &ctg,
                start,
                end,
                mean,
                lower,
                upper,
                Kge,
                Klt,
            )?;
        }
    }

    Ok(())
}
