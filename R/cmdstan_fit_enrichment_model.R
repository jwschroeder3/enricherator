#!/usr/bin/env Rscript

library(cmdstanr)
library(optparse)
library(tidyverse)
library(stringr)

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
        c("--stan_file"), type="character",
        help="The name of the file containing the stan model",
    ),
    make_option(
        c("--compiled_model"), type="character", default=NULL,
        help="The name of the binary file containing the compiled stan model",
    ),
    make_option(
        c("--draws_direc"), type="character",
        help="The name of the directory into which draws.csv will be written",
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
        c("--grad_samps"), type="integer", default=1,
        help="Sets the value to grad_samples for variational inference. Do not adjust from its default of 1 unless you have a good reason to do so."
    ),
    make_option(
        c("--cores"), type="integer", default=4,
        help="The number of cores to use for parallel chains when doing MCMC."
    ),
    make_option(
        c("--debug"), action="store_true", default=FALSE,
        help="include if debugging. Will use small version of the data."
    ),
    make_option(
        c("--debug_outdir"), type="character", default=NULL,
        help="Directory into which files for debugging will be written"
    )
)

print("Reading command line options")
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

options(mc.cores = opt$cores)

debug = opt$debug
debug_outdir = opt$debug_outdir

if (debug) {
    if (is.null(debug_outdir)) {
        stop("--debug was set without --debug_outdir. To debug, you must set --debug_outdir. Exiting now.")
    }
}

if (!is.null(opt$compiled_model)) {
    model_file = opt$compiled_model
    print(paste0("Reading pre-built model in ", model_file))
    sm = cmdstan_model(exe_file=model_file)
} else {
    model_file = opt$stan_file
    print(paste0("Compiling model in ", model_file))
    sm = cmdstan_model(model_file, cpp_options = list(stan_threads = TRUE))
}

print(paste0("Reading stan_list from ", opt$data_file))
load(opt$data_file)
print("Fitting model using variational inference")

newlist = list()
include_vars = c(
    "L","S","B","A","G","Q","C","alpha_prior","sample_x", "param_idx_arr", "beta_sub_idx",
    "Y", "libsize", "hs_df_beta", "hs_df_global_beta", "hs_df_slab_beta",
    "hs_df_gamma", "hs_df_global_gamma", "hs_df_slab_gamma",
    "hs_scale_global_gamma", "hs_scale_slab_gamma", 
    "hs_scale_global_beta", "hs_scale_slab_beta", "a_sub_L", "b_sub_L",
    "g_sub_L",
    "g_num_non_zero", "g_weights_vals", "g_col_accessor",
    "g_row_non_zero_number", 
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

