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
        c("-f", "--fit_file"), type="character",
        help="Name of RData file containing stanfit object",
    ),
    make_option(
        c("-s", "--data_file"), type="character",
        help="Name of RData file containing stan data list",
    ),
    make_option(
        c("-q", "--quantile_interval"), type="numeric", default=90,
        help="Quantile intervale to return.",
    ),
    make_option(
        c("-p", "--params"), type="character",
        help="Comma-separated list of parameter names to gather",
    ),
    make_option(
        c("-o", "--out_direc"), type="character",
        help="Name of directory to save output bedgraph files",
    )
)
 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.na(opt$params)) {
    stop("--params argument was not provided, but is mandatory. Exiting now.")
}

params = str_split(opt$params, ",", simplify=TRUE)[1,]
fit_file = opt$fit_file
data_file = opt$data_file
out_direc = opt$out_direc

if (is.na(fit_file)) {
    stop("--fit_file argument was not provided, but is mandatory. Exiting now.")
}
if (is.na(data_file)) {
    stop("--data_file argument was not provided, but is mandatory. Exiting now.")
}
if (is.na(out_direc)) {
    stop("--out_direc argument was not provided, but is mandatory. Exiting now.")
}

print("Reading data used by stan")
load(data_file)
print("Reading cmdstanVBfit object")
load(fit_file)

gather_vb_estimates(
    fit,
    stan_list,
    direc=opt$out_direc,
    interval=opt$quantile_interval,
    cmdstan=TRUE,
    params=params
)
#save(param_df, file=opt$out_file)
#precision = list_of_draws[["prec"]]
#precision_q = apply(precision, FUN=quant, MAR=2)
#save(precision_q, file="full_vb_precision.RData")
