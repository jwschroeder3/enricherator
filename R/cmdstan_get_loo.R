#!/usr/bin/env Rscript

library(tidyverse)
library(optparse)
library(cmdstanr)
library(loo)

options(stringsAsFactors = FALSE)
#source("stan_helpers.R")

option_list = list(
    make_option(
        c("-f", "--fit_file"), type="character",
        help="Name of RData file containing stanfit object",
    ),
    make_option(
        c("-o", "--out_file"), type="character",
        help="Name of file to save loo result",
    ),
    make_option(
        c("-c", "--cores"), type="integer",
        help="The number of cores to use",
    )
)
 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

fit_file = opt$fit_file

print("Reading cmdstanfit object")
load(fit_file)

print("Gathering point-wise log-likelihood")
log_lik = fit$draws("log_lik")
loo_res = loo::loo(log_lik, cores=opt$cores)
print(loo_res)
save(loo_res, file=opt$out_file)

