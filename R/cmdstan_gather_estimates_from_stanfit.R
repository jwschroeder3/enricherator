#!/usr/bin/env Rscript

library(cmdstanr)
library(optparse)
library(tidyverse)

options(stringsAsFactors = FALSE)
source("stan_helpers.R")

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
    #make_option(
    #    c("--method"), type="character",
    #    help="The name of the fitting method (can be either 'EM' or 'VB')",
    #)
)
 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

params = str_split(opt$params, ",", simplify=TRUE)[1,]
fit_file = opt$fit_file
data_file = opt$data_file

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
