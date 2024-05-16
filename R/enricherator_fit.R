#!/usr/bin/env Rscript

library(optparse)
library(stringr)
library(parallel)

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

fit_path = file.path(this_direc,"cmdstan_fit_enrichment_model.R")

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
        c("--out_direc"), type="character",
        help="The name of the directory into which all output will be written",
    ),
    make_option(c("--frac_genome_enriched"), type="double", default=NULL,
        help="Prior expectation for the fraction of the genome enriched (positive and negative) for the signal of interest. Sets the degree of shrinkage applied to enrichment estimates.",
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
        c("--cores"), type="integer", default=5,
        help="The number of cores to use for parallel chains when fitting model."
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
    ),
    make_option(
        c("--shared_input"), action="store_true", default=FALSE,
        help="Include at command line if input replicates are to be shared across ChIP-seq genotypes/conditions"
    )
)
print("Reading command line options")
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

options(mc.cores = opt$cores)


fit_enricherator = function(.i, .opt) {
    #print(.i)

    args = NA
    args_length=0

    for (arg_name in names(.opt)) {
        new_arg_name = arg_name
        arg = opt[[arg_name]]
        if (arg_name == "out_direc") {
            new_arg_name = "draws_direc"
            arg = paste0(opt[["out_direc"]], "/draws_", .i)
        }
        args_length = args_length + 1
        args[args_length] = paste0("--",new_arg_name)
        args_length = args_length + 1
        #print(arg)
        args[args_length] = arg
    }
    args_length = args_length + 1
    args[args_length] = "--fit_file"
    args_length = args_length + 1
    args[args_length] = paste0(.opt[["out_direc"]], "/fit_", .i, ".RData")
    args_length = args_length + 1
    args[args_length] = "--data_file"
    args_length = args_length + 1
    args[args_length] = paste0(.opt[["out_direc"]], "/data_", .i, ".RData")

    outfile = paste0(.opt[["out_direc"]], "/fit_", .i, ".log")
    errfile = paste0(.opt[["out_direc"]], "/fit_", .i, ".err")

    res = system2(
        "Rscript",
        c(
            fit_path,
            args
        ),
        stderr=errfile,
        stdout=outfile
    )
    if (res != 0) stop(paste0("Error in fit ", .i, ". Do not use results of this fit. Check logs carefully."))
}

choose_fit = function(inputs, .opt) {
    elbos = NA
    for (i in inputs) {
        outfile = paste0(.opt[["out_direc"]], "/fit_", i, ".log")
        fcon = file(outfile, "r")
        while(TRUE) {
            line = readLines(fcon, n=1)
            if (grepl("ELBO CONVERGED", line)) {
                break
            }
        }
        close(fcon)
        elbo = str_split(line, "\\s+", simplify=TRUE)[1,3]
        elbos[i] = as.numeric(elbo)
    }
    ordered = rev(order(elbos))
    return(ordered[1])
}

prune_fits = function(retained_fit, inputs, .opt) {

    final_fit_file = paste0(.opt[["out_direc"]], "/fit.RData")
    final_fit_log = paste0(.opt[["out_direc"]], "/fit.log")
    final_fit_err = paste0(.opt[["out_direc"]], "/fit.err")
    final_data_file = paste0(.opt[["out_direc"]], "/data.RData")

    for (i in inputs) {

        fit_file = paste0(.opt[["out_direc"]], "/fit_", i, ".RData")
        fit_log = paste0(.opt[["out_direc"]], "/fit_", i, ".log")
        fit_err = paste0(.opt[["out_direc"]], "/fit_", i, ".err")
        data_file = paste0(.opt[["out_direc"]], "/data_", i, ".RData")
        draws_direc = paste0(opt[["out_direc"]], "/draws_", i)

        if (i == retained_fit) {
            unlink(final_fit_file)
            unlink(final_fit_log)
            unlink(final_fit_err)
            unlink(final_data_file)
            file.rename(from=fit_file, to=final_fit_file)
            file.rename(from=fit_log, to=final_fit_log)
            file.rename(from=fit_err, to=final_fit_err)
            file.rename(from=data_file, to=final_data_file)
        } else {
            unlink(fit_file)
            unlink(fit_log)
            unlink(fit_err)
            unlink(data_file)
            unlink(draws_direc, recursive=TRUE)
        }
    }
}

inputs = list()
for (i in seq(5)) {
    inputs[i] = i
}

results = mclapply(inputs, fit_enricherator, opt, mc.cores = opt$cores)
retained_i = choose_fit(inputs, opt)
prune_fits(retained_i, inputs, opt)

