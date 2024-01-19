#!/usr/bin/env Rscript

library(cmdstanr)
library(optparse)
library(tidyverse)

get_exec_file = function() {
    # copied from https://stackoverflow.com/questions/1815606/determine-path-of-the-executing-script
    cmdArgs = commandArgs(trailingOnly = FALSE)
    needle = "--file="
    match = grep(needle, cmdArgs)
    if (length(match) > 0) {
        # Rscript
        return(normalizePath(sub(needle, "", cmdArgs[match])))
    } else {
        # 'source'd via R console
        return(normalizePath(sys.frames()[[1]]$ofile))
    }
}

options(stringsAsFactors = FALSE)
this_file = get_exec_file()
this_direc = dirname(this_file)
source(file.path(this_direc,"stan_helpers.R"))

option_list = list(
    make_option(
        c("-i", "--info"), type="character",
        help="File containing all sample information",
    ),
    make_option(
        c("--ignore_ctgs"), type="character", default=NULL,
        help="Comma-separated list of contigs to remove from analysis."
    ),
    make_option(
        c("--norm_method"), type="character", default="libsize",
        help="The normalization method to use. Can be one of 'spikein' or 'libsize'"
    ),
    make_option(
        c("--spikein"), type="character", default=NULL,
        help="The name of the contig to use to calculate library size (or spike-in) normalization factors",
    ),
    make_option(
        c("--spikein_rel_abund"), type="numeric", default=NULL,
        help="The fraction of reads expected to align to spikein contig. Set this carefully, it's used to set prior on alpha",
    ),
    make_option(
        c("--stan_file"), type="character",
        help="The name of the file containing the stan model. You must include either this argument or the --compiled_model argument.",
    ),
    make_option(
        c("--compiled_model"), type="character",
        help="The name of the executable file containing the stan model. You must include either this argument or the --stan_file model.",
    ),
    make_option(
        c("--draws_direc"), type="character",
        help="The name of the directory into which draws.csv will be written",
    ),
    make_option(c("--frac_genome_enriched"), type="double", default=NULL,
        help="Prior expectation for the fraction of the genome enriched (positive and negative) for the signal of interest. Sets the degree of shrinkage applied to enrichment estimates.",
    ),
    make_option(
        c("--seed"), type="integer",
        help="Sets the seed to the (P)RNG for variational inference. Not usually set, but helpful for testing reproducibility.",
    ),
    make_option(
        c("--fit_file"), type="character",
        help="The name of the RData file containing stanfit object",
    ),
    make_option(
        c("--data_file"), type="character",
        help="The name of the RData file containing stan data list",
    ),
    make_option(
        c("--ext_fragment_length"), type="integer", default=60,
        help="Sets the width of the smoothing kernel applied to sampled betas."
    ),
    make_option(
        c("--ext_subsample_dist"), type="integer", default=30,
        help="Distance (in bp) between sampled betas, which will later be mixed using overlapping exponential kernels defined by --ext_fragment_length."
    ),
    make_option(
        c("--input_fragment_length"), type="integer", default=120,
        help="Sets the width of the smoothing kernel applied to sampled betas."
    ),
    make_option(
        c("--input_subsample_dist"), type="integer", default=60,
        help="Distance (in bp) between sampled alphas, which will later be mixed using overlapping exponential kernels defined by --input_fragment_length."
    ),
    make_option(
        c("--cores"), type="integer", default=4,
        help="The number of cores to use for parallel chains when doing MCMC."
    ),
    make_option(
        c("--libsize_key"), type="character", default="tm_size_factors",
        help="The key in stan_list to use as the libsize for size factor exposure."
    ),
    make_option(
        c("--debug"), action="store_true", default=FALSE,
        help="include if debugging. Will use small version of the data."
    ),
    make_option(
        c("--log_lik"), action="store_true", default=FALSE,
        help="include if you want the point-wise log-likelihood calculated and returned from the stanfit object."
    ),
    make_option(
        c("--no_beta"), action="store_true", default=FALSE,
        help="include if you only want to use input data and fit only alphas."
    ),
    make_option(
        c("--load_data_file"), action="store_true", default=FALSE,
        help="include if you want to load data already prepared from --data_file"
    ),
    make_option(
        c("--genome_data_file"), type="character", default=NULL,
        help="RData file containing initial genome data"
    ),
    make_option(
        c("--grad_samps"), type="integer", default=1,
        help="Sets the value to grad_samples for variational inference. Do not adjust from its default of 1 unless you have a good reason to do so."
    )
)
print("Reading command line options")
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

options(mc.cores = opt$cores)

debug = opt$debug

model_file = opt$stan_file
bin_file = opt$compiled_model
if (is.null(model_file)) {
    if (is.null(bin_file)) {
        stop("Either the --stan_file argument or the --compiled_model argument is required, but neither was present. Exiting now.")
    }
}
if (is.null(bin_file)) {
    if (is.null(model_file)) {
        stop("Either the --stan_file argument or the --compiled_model argument is required, but neither was present. Exiting now.")
    }
}
if (!is.null(bin_file)) {
    if (!is.null(model_file)) {
        stop("You included both the --stan_file and the --compiled_model arguments to enricherator. Re-run, including only one of those arguments. Exiting now.")
    }
}
if (is.na(opt$draws_direc)) {
    stop("--draws_direc argument is required, but was absent. Exiting now.")
}
if (is.na(opt$fit_file)) {
    stop("--fit_file argument is required, but was absent. Exiting now.")
}
if (is.na(opt$data_file)) {
    stop("--data_file argument is required, but was absent. Exiting now.")
}

if (is.null(bin_file)) {
    print(paste0("Compiling model in ", model_file))
    sm = cmdstan_model(model_file, cpp_options = list(stan_threads = TRUE))
} else {
    print(paste0("Using precombiled model in ", bin_file))
    sm = cmdstan_model(exe_file=bin_file)
}

spikein = opt$spikein
if (!is.null(opt$ignore_ctgs)) {
    ignores = str_split(opt$ignore_ctgs, ",", simplify=TRUE)[1,]
} else {
    ignores = ""
}
load = opt$load_data_file

if (load) {
    load(opt$data_file)
} else {
    if (is.null(opt$genome_data_file)) {
        experiment_info = read_csv(opt$info) %>%
            mutate(
                sample_x = ifelse(sample=="input", 0, ifelse(sample=="inp", 0, 1)), # set hbd to 1
                strand_x = ifelse(strand=="both", 1, ifelse(strand=="plus", 1, 2)) #plus strand is idx 1, minus is idx2
            )

        if (opt$norm_method == "spikein") {
            if (is.null(opt$spikein_rel_abund)) {
                stop("You set --norm_method to 'spikein', but you did not include the --spikein_rel_abund argument at the command line. To do spikein normalization you must include both. Enricherator was not run. Exiting now.")
            }
            if (is.null(experiment_info$norm_factor)) {
                stop("You chose to do spike-in normalization, but no column named 'norm_factor' was found in the info file. Spike-in normalization factors should be a column in the info file.")
            }
        }

        genotype_factor = factor(experiment_info$genotype)
        experiment_info$geno_x = as.integer(genotype_factor)

        replicate_factor = factor(experiment_info$rep)
        experiment_info$rep_id = as.integer(replicate_factor)

        if (opt$no_beta) {
            experiment_info = experiment_info %>%
                filter(sample == "input")
        }

        experiment_info = experiment_info %>%
            mutate(
                sample_id = as.character(interaction(genotype,sample,rep)),
                data = purrr::map(file, read_file, ignores=ignores, debug=debug)
            )
        save(experiment_info, file="intermediate.RData")
    } else {
        load(opt$genome_data_file)
    }

    #save(experiment_info, file="debug_stan.RData")
    #load("debug_stan.RData")
    stan_list = prep_stan_data(
        experiment_info,
        opt$norm_method,
        spikein,
        opt$spikein_rel_abund,
        opt$input_subsample_dist,
        opt$input_fragment_length,
        opt$ext_subsample_dist,
        opt$ext_fragment_length,
        opt$log_lik,
        opt$frac_genome_enriched
    )
    save(stan_list, file=opt$data_file)
}

if (opt$norm_method == "libsize") {
    stan_list[["libsize"]] = stan_list[[opt$libsize_key]] 
} else if (opt$norm_method == "spikein") {
    stan_list[["libsize"]] = stan_list[["spikein_norm_factors"]]
}

print("Fitting model using variational inference")

newlist = list()
include_vars = c(
    "L","S","B","A","G","Q","alpha_prior","geno_x","sample_x",
    "Y", "libsize", "hs_df", "hs_df_global", "hs_df_slab",
    "hs_scale_global", "hs_scale_slab", "a_sub_L", "b_sub_L",
    "b_num_non_zero", "b_weights_vals", "b_col_accessor",
    "b_row_non_zero_number", "a_num_non_zero", "a_weights_vals",
    "a_col_accessor", "a_row_non_zero_number", "gather_log_lik"
)
for (var in include_vars) {
    newlist[[var]] = stan_list[[var]]
}
if (!dir.exists(opt$draws_direc)) {
    dir.create(opt$draws_direc, recursive=TRUE)
}
grad_samps = opt$grad_samps

fit = sm$variational(
    data = newlist,
    seed = opt$seed,
    threads = opt$cores,
    output_samples = 500,
    output_dir = opt$draws_direc,
    output_basename = "draws",
    algorithm = "meanfield",
    grad_samples = grad_samps
)

save(fit, file=opt$fit_file)

