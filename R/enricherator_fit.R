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
        c("--data_direc"), type="character",
        help="The name of the directory containing data.RData file",
    ),
    make_option(
        c("--stan_file"), type="character",
        help="The name of the file containing the stan model",
    ),
    make_option(
        c("--seed"), type="integer",
        help="Sets the seed to the (P)RNG for variational inference. Not usually set, but helpful for testing reproducibility.",
    ),
    make_option(
        c("--grad_samps"), type="integer", default=1,
        help="Sets the value to grad_samples for variational inference. Do not adjust from its default of 1 unless you have a good reason to do so."
    ),
    make_option(
        c("--cores"), type="integer", default=7,
        help="The number of cores to use for parallel fits with best chosen to move forward."
    )
)
print("Reading command line options")
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# set out_direc to directory above data_direc
opt[["out_direc"]] = dirname(opt[["data_direc"]])
options(mc.cores = opt$cores)

fit_enricherator = function(.i, .opt) {
    #print(.i)

    args = NA
    args_length=0

    for (arg_name in names(.opt)) {
        new_arg_name = arg_name
        arg = .opt[[arg_name]]
        if (arg_name == "out_direc") {
            new_arg_name = "draws_direc"
            arg = paste0(.opt[["out_direc"]], "/draws_", .i)
            dir.create(arg, recursive=TRUE)
        }
        if (arg_name == "data_direc") {
            new_arg_name = "data_file"
            arg = paste0(.opt[["out_direc"]], "/data/data.RData")
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
        fit_converged = FALSE
        while(TRUE) {
            line = readLines(fcon, n=1)
            if (length(line) == 0) {
                break
            }
            if (grepl("ELBO CONVERGED", line)) {
                fit_converged = TRUE
                break
            }
        }
        close(fcon)
        if (fit_converged) {
            elbo = str_split(line, "\\s+", simplify=TRUE)[1,3]
        } else {
            elbo=NA
        }
        elbos[i] = as.numeric(elbo)
    }
    ordered = rev(order(elbos, na.last=FALSE))
    return(ordered[1])
}

prune_fits = function(retained_fit, inputs, .opt) {

    final_fit_file = paste0(.opt[["out_direc"]], "/fit.RData")
    final_fit_log = paste0(.opt[["out_direc"]], "/fit.log")
    final_fit_err = paste0(.opt[["out_direc"]], "/fit.err")

    for (i in inputs) {

        fit_file = paste0(.opt[["out_direc"]], "/fit_", i, ".RData")
        fit_log = paste0(.opt[["out_direc"]], "/fit_", i, ".log")
        fit_err = paste0(.opt[["out_direc"]], "/fit_", i, ".err")
        draws_direc = paste0(opt[["out_direc"]], "/draws_", i)

        if (i == retained_fit) {
            unlink(final_fit_file)
            unlink(final_fit_log)
            unlink(final_fit_err)
            file.rename(from=fit_file, to=final_fit_file)
            file.rename(from=fit_log, to=final_fit_log)
            file.rename(from=fit_err, to=final_fit_err)
        } else {
            unlink(fit_file)
            unlink(fit_log)
            unlink(fit_err)
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

