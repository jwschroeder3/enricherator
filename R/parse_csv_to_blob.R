library(optparse)

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
        c("-c", "--cfg_file"), type="character",
        help="config.toml file with information on inputs, output, parsers, and contrasts of interest",
    )
)

print("Reading command line options")
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

parse_csv_to_blob(opt$cfg_file)

