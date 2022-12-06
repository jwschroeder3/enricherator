# Enricherator

A set of tools to use Bayesian inference to infer enrichment
of sequencing fragments in an extracted sample, i.e.,
Chip-seq, IPOD-HR, HBD-seq, etc., vs. input DNA fragments.

## Conda env setup

To set up the conda environment, run

```bash
conda env create -f conda_environment.yaml
```

This will create the `rstan` conda environment.

NOTE: The environment still will not contain `cmdstanr`,
and you will have to install that separately, from
within the environment.

To set up cmdstanr, run the following:

```bash
conda activate rstan
R
```

Then, from the R terminal, run:

```R
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# check whether cmdstanr is properly set up
library(cmdstanr)
check_cmdstan_toolchain()
```

The `check_cmdstan_toolchain` call has always returned that
the toolchain is properly set up, so I haven't had to
debug any issues that could arise here. However, if you
have issues, note that `check_cmdstan_toolchain` has
a `fix` argument that could be helpful.

## Running enrichment inference

Enter directory containing sample info file
(see `examples/nonstranded_sample_info.txt` and
`examples/stranded_sample_info.txt` for sample info file templates).
Run the following,
substituting the location of the `enricherator` source tree for
`<srcdir>`. You will also need to adjust several of the arguments
passed to the scripts below for your own usage. As a rule of thumb,
set `--ext_subsample_dist` to half the mean extracted data fragment length,
set `--ext_fragment_length` to the mean extracted data fragment length,
set `--input_subsample_dist` to half the mean input data fragment length,
and set `--input_fragment_length` to the mean input data fragment length.

If you want to ignore any contigs in your reference fasta file, you
can set those contigs as a comma-separated lest using the `--ignore_ctgs`
argument as shown below. To include all contigs, simply omit the
`--ignore_ctgs` argument.

The `--libsize_key` argument should be exactly as shown below in order
to use a trimmed mean coverage as the size factor for normalization.
Also, since we don't yet have spike-in normalization working as
we'd like, keep the `--norm_method` argument set to `libsize`.

Also note that whether your data and analysis are stranded or not,
you'll use `ab_smoothed_stranded_enrichment.stan` as your `--stanfile`
argument.

```bash
cd <top_direc>

mkdir enricherator_results
mkdir enricherator_results/draws

SRCDIR="<srcdir>"

Rscript $SRCDIR/R/cmdstan_fit_enrichment_model.R \
    --info stan_sample_info.txt \
    --stan_file "${SRCDIR}/stan/ab_smoothed_stranded_enrichment.stan" \
    --ignore_ctgs P2918_rnadna_spikein,NC_011916.1 \
    --norm_method libsize \
    --fit_file "enricherator_results/fit.RData" \
    --data_file "enricherator_results/data.RData" \
    --ext_subsample_dist 30 \
    --ext_fragment_length 60 \
    --input_subsample_dist 60 \
    --input_fragment_length 120 \
    --libsize_key "tmm_size_factors" \
    --draws_direc "enricherator_results/draws" \
    --cores 32 \
    > "enricherator_results/fit.log" \
    2> "enricherator_results/fit.err"
```

Once complete, the following files should be present in `enricherator_results`:

```bash
fit.log
fit.err
fit.RData
data.RData
draws/draws-1.csv
```

Note that if you fit the model to genome-wide data,
`draws-1.csv` will be several hundreds of Gigabytes,
since it is the full output for 500 samples of the
approximate posterior for tens- to hundreds-of-millions of
parameters. These parameters are filtered, and the 500 samples
are summarized in the next step.

## Extracting quantities of interest from the cmdstan output

After the fitting step completes, run the following
to get summaries for each parameter of interest.
We write a bedgraph file for each of the following summary statistics
for the 500 samples from the approximate posterior:
mean, median, lower 90% quantile interval, and upper 90% quantile interval.
Bedgraph files are written for each genotype and strand in your analysis.

```bash
cd <top_direc>
mkdir enricherator_results/out_files
SRCDIR="<srcdir>"

Rscript $SRCDIR/R/cmdstan_gather_estimates_from_stanfit.R \
        --fit_file "enricherator_results/fit.RData" \
        --data_file "enricherator_results/data.RData" \
        --out_direc "enricherator_results/out_files" \
        --params Alpha,Beta \
        > "enricherator_results/gather.log" \
        2> "enricherator_results/gather.err"
```

The bedgraph files will appear in `enricherator_results/out_files`.
Robust z-scores and rolling means can be calculated using `bgtools`.

At this point, if you have inspected the quantities output by the
above script, and if they seem reasonable, it is recommended that you
delete the `draws-1.csv` file.

