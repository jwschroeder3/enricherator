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
        help="Name of RData file containing stan data list",
        default=NA
    ),
    make_option(
        c("-s", "--samples_file"), type="character",
        help="Name of csv file containing samples",
        default=NA
    ),
    make_option(
        c("-c", "--contrasts"), type="character",
        help="Comma-separated list of contrasts to perform",
        default=NA
    ),
    make_option(
        c("-t", "--type"), type="character",
        help="Set whether to do genotype or strand contrasts. Allowed values are 'genotype' and 'strand'",
        default=NA
    ),
    make_option(
        c("--out_direc"), type="character",
        help="Directory into which to write output bedgraph files.",
        default=NA
    ),
    make_option(
        c("--parameter"), type="character", default="Beta",
        help="Name of the parameter to perform contrasts for. Defaults to Beta.",
    )
)
 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

contrasts = opt$contrasts
con_type = opt$type
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
if (is.na(con_type)) {
    stop("--type argument was not provided, but is mandatory. Exiting now.")
}
# if con_type is neither of the possible choices, exit.
if (!(con_type == "genotype") & (!(con_type == "strand"))) {
    stop(paste0(
        "--type argument must be either 'genotype' or 'strand', but you provided, '",
        con_type,
        "'. Exiting now."
    ))
}

if (!dir.exists(opt$out_dir)) {
    dir.create(opt$out_dir, recursive=TRUE)
}

print("Reading data used by stan")
load(data_file)

draws_direc = dirname(samples)
contrast_fname = file.path(draws_direc, "contrast_summaries.txt")

data_direc = dirname(data_file)
geno_fname = file.path(data_direc, "geno_lut.csv")
strand_fname = file.path(data_direc, "strand_lut.csv")

info = stan_list[["info"]]
geno_lut = info[!duplicated(info$geno_x), c("genotype","geno_x")]
strand_lut = info[!duplicated(info$strand_x), c("strand","strand_x")]
write.csv(geno_lut, geno_fname, row.names=FALSE, quote=FALSE)
write.csv(strand_lut, strand_fname, row.names=FALSE, quote=FALSE)

exec_f_path = get_exec_file()
exec_direc = dirname(exec_f_path)
bin_path = file.path(exec_direc, "../bin/get_contrasts")
if (con_type == "genotype") {
    lut_fname = geno_fname
}
if (con_type == "strand") {
    lut_fname = strand_fname
}

args_vec = c(
    samples,
    contrast_fname,
    "501",
    param,
    contrasts,
    lut_fname,
    con_type
)

exec = c(bin_path, args_vec)
print(paste0("Running the following: ", paste(exec, collapse=" ")))

res = system2(
    bin_path,
    args_vec,
    env="RUST_BACKTRACE=1",
    stderr="",
    stdout=""
)
if (res != 0) {
    stop("Error getting contrasts. Check err and log files. Exiting now.")
}
param_summaries = read_tsv(contrast_fname)
#print(head(param_summaries))
write_cmdstan_summaries(param_summaries, stan_list, out_dir, param, contrasts, con_type)

