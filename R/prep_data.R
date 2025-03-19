#!/usr/bin/env Rscript

library(cmdstanr)
library(optparse)
library(tidyverse)
library(stringr)
library(blogdown)

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
        c("--draws_direc"), type="character",
        help="The name of the directory into which draws.csv will be written",
    ),
    make_option(c("--frac_genome_enriched"), type="double", default=NULL,
        help="Prior expectation for the fraction of the genome enriched (positive and negative) for the signal of interest. Sets the degree of shrinkage applied to enrichment estimates.",
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
        c("--ipod_fragment_length"), type="integer", default=90,
        help="Sets the width of the smoothing kernel applied to ipod fits"
    ),
    make_option(
        c("--ipod_subsample_dist"), type="integer", default=45,
        help="Distance (in bp) between sampled ipod factors."
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
        c("--debug_outdir"), type="character", default=NULL,
        help="Directory into which files for debugging will be written"
    ),
    make_option(
        c("--log_lik"), action="store_true", default=FALSE,
        help="include if you want the point-wise log-likelihood calculated and returned from the stanfit object."
    ),
    make_option(
        c("--load_data_file"), action="store_true", default=FALSE,
        help="include if you want to load data already prepared from --data_file"
    ),
    make_option(
        c("--genome_data_file"), type="character", default=NULL,
        help="RData file containing initial genome data"
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

if (is.null(opt$draws_direc)) {
    stop("--draws_direc argument is required, but was absent. Exiting now.")
}

# Check if the directory exists
if (!file.exists(opt$draws_direc)) {
      # Create the directory
      dir.create(opt$draws_direc)
  cat("Directory created:", opt$draws_direc, "\n")
} else {
      cat("Directory already exists:", opt$draws_direc, "\n")
}

if (is.null(opt$data_file)) {
    stop("--data_file argument is required, but was absent. Exiting now.")
}

spikein = opt$spikein
if (!is.null(opt$ignore_ctgs)) {
    ignores = str_split(opt$ignore_ctgs, ",", simplify=TRUE)[1,]
} else {
    ignores = ""
}
if (is.null(opt$genome_data_file)) {
    experiment_info = read_csv(opt$info) %>%
        mutate(
            sample_x = ifelse(sample=="input", 0, ifelse(sample=="inp", 0, ifelse(sample=="extract", 2, 1))),
            strand_x = ifelse(
                strand=="both", 1,
                ifelse(
                    strand=="plus", 1,
                    ifelse(
                        strand=="minus", 2,
                        stop(paste0("For a non-stranded analysis, every value in the 'strand' column of your info file should be 'both'. For a stranded analysis, you must have either 'plus' or 'minus' in the strand column. Instead, ", strand, " was found. Exiting without fitting a model. Update your info file and re-run.")))))
        )

    # Check for disallowed characters in the genotype column
    if (any(str_detect(experiment_info$genotype, "[-+]"))) {
        print()
        stop("Error: '-' or '+' characters are not allowed in genotype names. Replace all instances of those characters and try again. Exiting now.")
    }

    if (opt$norm_method == "spikein") {
        if (is.null(opt$spikein_rel_abund)) {
            stop("You set --norm_method to 'spikein', but you did not include the --spikein_rel_abund argument at the command line. To do spikein normalization you must include both. Enricherator was not run. Exiting now.")
        }
        if (is.null(experiment_info$norm_factor)) {
            stop("You chose to do spike-in normalization, but no column named 'norm_factor' was found in the info file. Spike-in normalization factors should be a column in the info file.")
        }
    }

    replicate_factor = factor(experiment_info$rep)
    experiment_info$rep_id = as.integer(replicate_factor)

    experiment_info = experiment_info %>%
        mutate(
            sample_id = as.character(interaction(genotype,sample,rep,input_group,control_group)),
            data = purrr::map(file, read_file, ignores=ignores, debug=debug)
        )
    if (debug) {
        save(experiment_info, file=file.path(debug_outdir, "intermediate.RData"))
    }
} else {
    load(opt$genome_data_file)
}

print(opt)
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
    opt$ipod_subsample_dist,
    opt$ipod_fragment_length,
    opt$log_lik,
    opt$frac_genome_enriched,
    opt$debug,
    opt$debug_outdir
)
stan_list[["expt_key_file"]] = paste0(opt$draws_direc, "/expt_key.json")
stan_list[["pos_map_file"]] = paste0(opt$draws_direc, "/position_map.txt")
stan_list[["config_file"]] = paste0(opt$draws_direc, "/config.toml")

print(paste0("Writing covariate key to ", stan_list[["expt_key_file"]]))
jsonlite::toJSON(stan_list[["expt_key"]]) %>%
    jsonlite::prettify() %>%
    write(stan_list[["expt_key_file"]])

print(paste0("Writing position map to ", stan_list[["pos_map_file"]]))
write.table(
    stan_list[["position_mapper"]],
    file=stan_list[["pos_map_file"]],
    col.names=FALSE,
    row.names=FALSE,
    quote=FALSE
)

if (opt$norm_method == "libsize") {
    stan_list[["libsize"]] = stan_list[[opt$libsize_key]] 
} else if (opt$norm_method == "spikein") {
    stan_list[["libsize"]] = stan_list[["spikein_norm_factors"]]
}

print(paste0("Writing config file to ", stan_list[["config_file"]]))
cfg_list = list(
    "input" = list(
        "expt_key_file" = stan_list[["expt_key_file"]],
        "draws_file" = paste0(opt$draws_direc, "/draws.csv"),
        "pos_file_name" = stan_list[["pos_map_file"]]
    ),
    "output" = list(
        "out_direc" = opt$draws_direc
    ),
    "beta_sub_idx" = stan_list[["beta_sub_idx"]]
)

write_toml(cfg_list, output = stan_list[["config_file"]])
print(paste0("Writing stan_list to ", opt$data_file))
save(stan_list, file=opt$data_file)

