# Enricherator

A set of tools to use Bayesian inference to infer enrichment
of sequencing fragments in an extracted sample, i.e.,
Chip-seq, IPOD-HR, HBD-seq, etc., vs. input DNA fragments.

For the under-the-hood workings of Enricherator,
see supplemental material in the original Enricherator [paper](https://www.science.org/doi/10.1126/sciadv.adi5945),
which also happens to be the work you should cite if using Enricherator in your own research:

    Schroeder, et al. 2023. Science Advances 9 (30): eadi5945. 

## Conda env setup

To start using Enricherator, clone this repository and enter
the directory created.

Then, to set up the conda environment, run

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
`examples/stranded_sample_info.txt` for sample info file templates.
Note: the norm_factor column of these files is currently ignored, but
just has to be present.).
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
you'll use `sig_noise_alloc.stan` as your `--stanfile`
argument.

```bash
cd <top_direc>
conda activate rstan

mkdir enricherator_results
mkdir enricherator_results/draws

SRCDIR="<srcdir>"

Rscript $SRCDIR/R/cmdstan_fit_enrichment_model.R \
    --info stan_sample_info.txt \
    --stan_file "${SRCDIR}/stan/sig_noise_alloc.stan" \
    --ignore_ctgs P2918_rnadna_spikein,NC_011916.1 \
    --norm_method libsize \
    --fit_file "enricherator_results/fit.RData" \
    --data_file "enricherator_results/data.RData" \
    --ext_subsample_dist 30 \
    --ext_fragment_length 60 \
    --input_subsample_dist 60 \
    --input_fragment_length 120 \
    --libsize_key "tm_size_factors" \
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
to filter the large csv file, keeping the parameters of interest, and to
get summaries for each parameter of interest.
We write a bedgraph file for each of the following summary statistics for each
genotype and strand:
mean of the 500 samples, median of the 500 samples, lower 90% quantile of the 
500 samples, and the upper 90% quantile of the 500 samples.

```bash
cd <top_direc>
conda activate rstan
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

Several new files will be produced by `cmdstan_gather_estimates_from_stanfit.R`.

1. Bedgraph files containing scores on the log2 scale will appear in `enricherator_results/out_files`. The files with "Beta" in their name contain log2(enrichment) scores, and the files with "Alpha" in their name are a measure of the log2(normalized input abundance). Robust z-scores, rolling medians, and rolling means can be calculated using [`bgtools`](https://github.com/jwschroeder3/bgtools).
2. The file `enricherator_results/draws/summaries.txt` contains the information also found in the bedgraph files already mentioned, and is one output of the binary `$SRCDIR/bin/parse_cmdstan_csv`
3. The file `enricherator_results/draws/samples.csv` contains the samples for only the parameters of interest (indicated by the `--params` argument to the `$SRCDIR/R/cmdstan_gather_estimates_from_stanfit.R` script). It is also an output of the binary `$SRCDIR/bin/parse_cmdstan_csv`. We typically keep this file for potential later use in setting up new contrasts of interest. For instance, if your analysis contains multiple genotypes or is stranded, you could set up a contrast to compare enrichments across those genotypes or across strands.

At this point, if you have inspected the quantities output by the
above script, and if they seem reasonable, it is recommended that you
delete the `draws-1.csv` file.

## Running additional contrasts (currently only supported for genotype-level or strand-level contrasts)

To run contrasts using the samples from the approximate posterior, do the following:

```bash
cd <top_direc>
conda activate rstan
mkdir enricherator_results/contrasts
SRCDIR="<srcdir>"

Rscript $SRCDIR/R/get_contrasts.R \
    --type <contrast_type> \
    --data_file enricherator_results/data.RData \
    --samples_file enricherator_results/draws/samples.csv \
    --contrasts <contrast_arg> \
    --out_direc enricherator_results/contrasts \
    > enricherator_results/contrast.log \
    2> enricherator_results/contrast.err
```

Note that `<contrast_type>` must be replaced with the type of
contrast you're performing, i.e., "genotype" or "strand", and
`<contrast_arg>` must be replaced by the contrasts of interest.
For example, if you had three genotypes in your analysis,
you could supply as the `--type` argument `genotype`, and for the
`--contrasts` argument `genoB-genoA,genoC-genoA,genoC-genoB`.
This would perform pairwise comparisons of enrichments from genotype B
to genotype A, genotype C to genotype A, and genotype C to genotype A,
respectively. The list of contrasts can be arbitrarily long.
The contrasts must be comma-separated with no spaces.

Similarly, to perform a contrast of minus strand over plus strand,
you would set the `--type` argument to "strand" and the `--contrasts`
argument to `minus-plus`, assuming your original sample info file
assigned the minus strand as "minus" and the plus strand as "plus".
