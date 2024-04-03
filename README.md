# Enricherator

A set of tools to use Bayesian inference to infer enrichment
of sequencing fragments in an extracted sample, i.e.,
Chip-seq, IPOD-HR, HBD-seq, etc., vs. input DNA fragments.

We note that Enricherator takes the term "sequencing fragments"
seriously. Its statistical model works well for count data that quantify alignments of
paried-end sequencing fragments. If you have single-end data, you must
estimate the fragment size and extend the 3-prime ends of your
alignments by the appropriate ammount to count pseudo-fragments.


The current version of Enricherator is in press at the Journal of Molecular Biology
and will be published in the 2024 special issue "Computation Resources for Molecular Biology"


The original implementation of Enricherator was published in this [paper](https://www.science.org/doi/10.1126/sciadv.adi5945), but the current version is substantially improved.
Enricherator now allocates sequencing counts between "signal" and "noise" components.

## Enricherator workflow

Enricherator uses bedgraph files (in the near future we look forward to
supporting use of bigwig files)
containing fragment counts to infer enrichments.
To use enricherator, the user must have counted paired-end alignment counts
at their desired genome resolution (we often use 5-bp resolution), and have
estimated the mean fragment size for their extracted data and their input
sequencing data. Enricherator will use the fragment sizes to inform
how much local genomic space over which to smooth enrichment estimates.

Once the above information and files are in hand, an "info file" must
be prepared to indicate to Enircherator the salient information about
each sequencing library, i.e., the condition of interest, the biological repliate
id a given library represents, and the strand represented by the data if a
strand-specific analysis is desired.

Examples of bedgraph files and info files can be found in the "examples"
directory of this repository.

The steps of running enricherator include:

1. Fit the Enricherator model to data
    * Enricherator will fit the model 5 independent times, concurrently,
    and will choose to proceed with the fit that achieved the highest ELBO.
2. Gather enrichment summaries
    * Bedgraph files will be written containing enrichments
3. Calculate contrasts (optional)
    * Bedgraph files will be written containing estimates of the contrast
    of interest. Currently Enricherator provides the ability to contrast
    across conditions (drug vs control or knock-out vs wild type)
    OR between strands.

## Using containerized Enricherator

We provide a containerized version of Enricherator that should simplify its use. The container
is build using Apptainer, which will allow the container to run on any Linux operating system
with Apptainer installed. 

To pull the container from its public repository, install Apptainer by following
its [online instructions](https://apptainer.org/docs/user/main/quick_start.html).

Once Apptainer is installed, add cloud.sylabs.io as a remote endpoint by running the
following:

```bash
apptainer remote add SylabsCloud cloud.sycloud.io
apptainer remote use SylabsCloud
```

Now you may pull the Enricherator container from the repository by
running the following, replacing `<container/path>` with the
location in which you would like the container to reside:

```bash
cd </container/path>
apptainer pull enricherator.sif library://schroedj/appliances/enricherator:latest
```

## Building Enricherator (only recommended if you cannot use the Apptainer container)

If you are able to use our containerized application, skip to the [Running enrichment inference](#running-enrichment-inference)
section.

### Conda env setup

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
Note: the norm\_factor column of these files is currently ignored, but
just has to be present. In future versions it may be used to
enable spike-in normalization as opposed to Enricherator's current
library size normalization).
Run the code below; you will need to adjust several of the arguments
passed to the scripts below for your own usage. First,
set `--ext_subsample_dist` to half the mean extracted data fragment length (`<ext_frag_len/2>` below),
set `--ext_fragment_length` to the mean extracted data fragment length (`<ext_frag_len>` below),
set `--input_subsample_dist` to half the mean input data fragment length (`<inp_frag_len/2>` below),
and set `--input_fragment_length` to the mean input data fragment length (`<inp_frag_len>` below).

If you want to ignore any contigs in your reference fasta file, you
can set those contigs as a comma-separated lest using the `--ignore_ctgs`
argument as shown below. To include all contigs, simply omit the
`--ignore_ctgs` argument.

In future versions of Enricherator we plan to support spike-in normalization,
but for now, the `--libsize_key` argument should be exactly as
shown below in order
to use a trimmed mean coverage as the size factor for normalization.

If not using apptainer, substitute your location of the Enricherator
source code for "/src" in the code examples below and omit the
line beginning with "apptainer".

```bash
cd <top_direc>

SRCDIR="/src"

# this script will perform 5 independent fits of the enricherator model
# to the input data. It can transiently require a lot of storage space (100s of GB)
apptainer exec -B $(pwd) /path/to/enricherator.sif \
    Rscript $SRCDIR/R/enricherator_fit.R \
    --info stan_sample_info.txt \
    --compiled_model ${SRCDIR}/stan/sig_noise_alloc \
    --ignore_ctgs P2918_rnadna_spikein,NC_011916.1 \
    --norm_method libsize \
    --out_direc enricherator_results \
    --ext_subsample_dist <ext_frag_len/2> \
    --ext_fragment_length <ext_frag_len> \
    --input_subsample_dist <inp_frag_len/2> \
    --input_fragment_length <inp_frag_len> \
    --libsize_key tm_size_factors
```

Once complete, the following files should be present in `enricherator_results`,
where the `i` in `draws_i` will be replaced with a number from 1-5, depending
on which of the 5 fits had the highest ELBO:

```bash
fit.log
fit.err
fit.RData
data.RData
draws_i/draws-1.csv
```

Note that if you fit the model to genome-wide data,
`draws-1.csv` will be several hundreds of Gigabytes,
since it is the full output for 500 samples of the
approximate posterior for tens- to hundreds-of-millions of
parameters. These parameters are filtered, and the 500 samples
are summarized in the next step, after which we typically
delete `draws-1.csv`.

## Extracting quantities of interest from the cmdstan output

After the fitting step completes, run the following
to filter the large csv file, keeping the parameters of interest, and to
get summaries for each parameter of interest.
We write a bedgraph file for each of the following summary statistics for each
genotype and strand:
mean of the 500 samples, median of the 500 samples, lower 90% quantile of the 
500 samples, the upper 90% quantile of the 500 samples, the evidence ratio
(K) for the coefficient being greater than the value provided at the
`--threshold` argument,
and K for the coefficient being less than the `--threshold` value.
Note that `--threshold` is on the log2-scale, so setting `--threshold 1.0`
denotes a fold-enrichment of 2, and setting `--threshold 2.0` would
denote a fold-enrichment of 4.

```bash
cd <top_direc>
SRCDIR="/src"

apptainer exec -B $(pwd) /path/to/enricherator.sif \
    Rscript $SRCDIR/R/cmdstan_gather_estimates_from_stanfit.R \
    --fit_file enricherator_results/fit.RData \
    --data_file enricherator_results/data.RData \
    --out_direc enricherator_results/out_files \
    --params Alpha,Beta \
    --threshold 1.0 \
    > enricherator_results/gather.log \
    2> enricherator_results/gather.err
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
SRCDIR="/src"

apptainer exec -B $(pwd) /path/to/enricherator_<version>.sif \
    Rscript $SRCDIR/R/get_contrasts.R \
    --type <contrast_type> \
    --data_file enricherator_results/data.RData \
    --samples_file enricherator_results/draws/samples.csv \
    --contrasts <contrast_arg> \
    --out_direc enricherator_results/contrasts \
    --threshold <your_threshold> \
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

Consider your choice of threshold carefully above. If your hypothesis
is that knocking out a protein will decrease enrichment, you likely
will want to set `--threshold` to be a negative number to collect
evidence ratios testing whether the knockout does, indeed, come with
less occupancy.

## Running on provided example data

To run Enricherator on our provided example simulated datasets, assuming your
current working directory is your local clone of this repository, do the following
for a non-stranded analysis:

```bash
cd examples/nonstranded

SRCDIR="/src"
OUTDIR="example_results"

apptainer exec -B $(pwd) /path/to/enricherator.sif \
    Rscript $SRCDIR/R/enricherator_fit.R \
    --info example_info.csv \
    --compiled_model ${SRCDIR}/stan/sig_noise_alloc \
    --norm_method libsize \
    --out_direc ${OUTDIR} \
    --ext_subsample_dist 30 \
    --ext_fragment_length 60 \
    --input_subsample_dist 60 \
    --input_fragment_length 120 \
    --libsize_key tm_size_factors

apptainer exec -B $(pwd) /path/to/enricherator.sif \
    Rscript $SRCDIR/R/cmdstan_gather_estimates_from_stanfit.R \
    --fit_file ${OUTDIR}/fit.RData \
    --data_file ${OUTDIR}/data.RData \
    --out_direc ${OUTDIR}/out_files \
    --params Alpha,Beta \
    > ${OUTDIR}/gather.log \
    2> ${OUTDIR}/gather.err
```

or the following for a stranded analysis:

```bash
cd examples/stranded

SRCDIR="/src"
OUTDIR="example_results"

apptainer exec -B $(pwd) /path/to/enricherator.sif \
    Rscript $SRCDIR/R/enricherator_fit.R \
    --info example_info.csv \
    --compiled_model ${SRCDIR}/stan/sig_noise_alloc \
    --norm_method libsize \
    --out_direc ${OUTDIR} \
    --ext_subsample_dist 30 \
    --ext_fragment_length 60 \
    --input_subsample_dist 60 \
    --input_fragment_length 120 \
    --libsize_key tm_size_factors

apptainer exec -B $(pwd) /path/to/enricherator.sif \
    Rscript $SRCDIR/R/cmdstan_gather_estimates_from_stanfit.R \
    --fit_file ${OUTDIR}/fit.RData \
    --data_file ${OUTDIR}/data.RData \
    --out_direc ${OUTDIR}/out_files \
    --params Alpha,Beta \
    --threshold 1.0 \
    > ${OUTDIR}/gather.log \
    2> ${OUTDIR}/gather.err
```

## Fully Bayesian peak calling using evidence ratios

The evidence ratios provided by both the `cmdstan_gather_estimates_from_stanfit.R`
script and the `get_contrasts.R` script can be used to call peaks at any desired
level of evidence. We often use the evidence ratios pointed out by
Kass and Raftery[^1], with K $\geq$ 3 representing "positive" evidence,
K $\geq$ 20 representing "strong" evidence, and K $\geq$ 150 representing
"very strong" evidence in favor of the hypothesis.

While we do not provide peak calling utilities directly in this repository
or in the Enricherator apptainer container, we share our basic workflow for peak calling
below.

### Calling peaks

Let us say that we ran Enricherator on ChIP-seq data and want to call peaks
in enrichment, and that we define our threshold on the log2-scale to be 1.0.
To call peaks we would then need to identify contiguous regions of the genome
with an evidence ratio supporting the hypothesis for enrichment above 1.0 that
exceeds our evidence threshold.

For this example, we will require "very strong" evidence (K $\geq$ 150 ) for enrichment
to call peaks. We will also remove any enriched regions narrower than 50 base pairs,
simply as an example here to demonstrate how that could be done. We use 
[`bgtools`](https://github.com/jwschroeder3/bgtools)
in this example to merge contigous bedgraph files into bed formatted regions.

```bash
OUTDIR="enricherator_results/out_files"
K_thresh=150
min_width=50

# we will call peaks in every condition using the following loop
# The bedgraph files with K_gt in their names contain the evidence ratios for
# samples from the approximate posterior being greater than the threshold, so
# those are the files we'll use for peak calling
for bgfile in ${OUTDIR}/*_K_gt.bedgraph; do
    # strip the file extension and path to generate new file name
    base=$(basename $bgfile .bedgraph)
    # get the path only
    direc=$(dirname $bgfile)
    # now make the new file path
    outfile="${direc}/${base}_Kgeq_${K_thresh}.bed"
    echo "Creating ${outfile} by filtering ${bgfile}"
    # use awk to print only lines exceeding the evidence threshold,
    # pipe the result to bgtools to get just the contiguous regions as bed format
    # pipe the bed formatted regions to awk to filter by size.
    awk -v thresh=$K_thresh '$4 >= thresh {print}' $bgfile \
        | bgtools contiguous_regions -i - \
        | awk -v width=$min_width 'BEGIN{OFS=FS="\t"} $3-$2 >= width {print}' \
        > $outfile
done
```

Similarly, negative peaks can be called. We find this useful for assessing evidence
for a given genotype of condition _decreasing_ enrichment of the signal of interest.
To identify regions greater than or equal to 50 bp wide with "strong evidence"
of decreased enrichment in genotype A vs genotype B, the code would look something like
the following:

```bash
OUTDIR="enricherator_results/contrasts"
K_thresh=150
min_width=50

# The bedgraph files with K_lt in their names contain the evidence ratios for
# samples from the approximate posterior being less than the threshold, so
# those are the files we'll use to test for loss of enrichment
bgfile="${OUTDIR}/genoA-genoB_both_strand_Beta_K_lt.bedgraph"

# strip the file extension and path to generate new file name
base=$(basename $bgfile .bedgraph)

# now make the new file path
outfile="${OUTDIR}/${base}_Kgeq_${K_thresh}.bed"
echo "Creating ${outfile} by filtering ${bgfile}"

# use awk to print only lines exceeding the evidence threshold,
# pipe the result to bgtools to get just the contiguous regions as bed format
# pipe the bed formatted regions to awk to filter by size.
awk -v thresh=$K_thresh '$4 >= thresh {print}' $bgfile \
    | bgtools contiguous_regions -i - \
    | awk -v width=$min_width 'BEGIN{OFS=FS="\t"} $3-$2 >= width {print}' \
    > $outfile
```

[^1]: Kass, R. E., & Raftery, A. E. (1995). Bayes Factors. Journal of the American Statistical Association, 90(430), 773–795. https://doi.org/10.1080/01621459.1995.10476572

