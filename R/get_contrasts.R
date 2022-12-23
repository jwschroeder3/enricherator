#!/usr/bin/env Rscript

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
        c("-d", "--data_file"), type="character",
        help="Name of RData file containing stan data list"
    ),
    make_option(
        c("-s", "--samples_file"), type="character",
        help="Name of csv file containing samples"
    ),
    make_option(
        c("-c", "--contrasts"), type="character",
        help="Comma-separated list of contrasts to perform"
    ),
    make_option(
        c("--out_direc"), type="character",
        help="Directory into which to write output bedgraph files."
    )
)
 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

contrasts = opt$contrasts
data_file = opt$data_file
samples = opt$samples_file
out_dir = opt$out_direc

if (is.na(out_dir)) {
    stop("--out_direc argument was not provided, but is mandatory. Exiting now.")
}
if (is.na(data_file)) {
    stop("--data_file argument was not provided, but is mandatory. Exiting now.")
}
if (is.na(samples)) {
    stop("--samples_file argument was not provided, but is mandatory. Exiting now.")
}
if (is.na(contrasts)) {
    stop("--contrasts argument was not provided, but is mandatory. Exiting now.")
}

print("Reading data used by stan")
load(data_file)

draws_direc = dirname(samples)
contrast_fname = file.path(draws_direc, "contrast_summaries.txt")

data_direc = dirname(data_file)
geno_fname = file.path(data_direc, "geno_lut.csv")

info = stan_list[["info"]]
geno_lut = info[!duplicated(info$geno_x), c("genotype","geno_x")]
write.csv(geno_lut, geno_fname, row.names=FALSE, quote=FALSE)

exec_f_path = get_exec_file()
exec_direc = dirname(exec_f_path)
bin_path = file.path(exec_direc, "../bin/get_contrasts")

args_vec = c(
    samples,
    contrast_fname,
    "501",
    "Beta",
    contrasts,
    geno_fname
)

exec = c(bin_path, args_vec)
print(paste0("Running the following: ", paste(exec, collapse=" ")))

res = system2(
    bin_path,
    args_vec,
    env="RUST_BACKTRACE=1",
    stderr="contrast.err",
    stdout="contrast.log"
)
if (res != 0) {
    stop("Error running get_contrasts. Check contrast.err and contrast.log. Exiting now.")
}
param_summaries = read_tsv(contrast_fname)
write_cmdstan_summaries(contrast_fname, stan_list, out_dir, "Beta", contrasts)

