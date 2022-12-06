#thisFile = function() {
#    cmdArgs = commandArgs(trailingOnly = FALSE)
#    needle = "--file="
#    match = grep(needle, cmdArgs)
#    if (length(match) > 0) {
#            # Rscript
#            return(normalizePath(sub(needle, "", cmdArgs[match])))
#    } else {
#            # 'source'd via R console
#            return(normalizePath(sys.frames()[[1]]$ofile))
#    }
#}

source("stan_helpers.R")

print(get_exec_file())
