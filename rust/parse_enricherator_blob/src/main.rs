use parse_csv_to_blob::{
    config::{Args, ExptKey, ConfigFile},
    index::Index,
    blob::Blob,
    BoxedResult,
};

fn main() -> BoxedResult<()> {

    let args = Args::parse();
    let expt_key = ExptKey::read(&args.expt_key)?;
    let cfg = ConfigFile::read(&args.cfg_file)?;
    println!("Config: {:?}", &cfg);
    println!("expt_key: {:?}", &expt_key);

    let cubic_params = [
        "a_cubic",
        "b_cubic",
        "c_cubic",
    ];

    let out_direc: &str = &cfg.output.out_direc;
    let offset_idx = Index::read(&out_direc)?;
    let blob_fname = format!("{out_direc}/params.blob");
    let buffer_size: usize = expt_key.sample_num[0];

    // instantiate the Blob to manage IO to binary file
    let mut blob_file = Blob::open(
        &blob_fname,
        buffer_size,
    )?;


}
