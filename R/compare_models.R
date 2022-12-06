#!/usr/bin/env Rscript

library(loo)
library(tidyverse)
library(optparse)

read_loo_result = function(.fname) {
    return(get(load(.fname)))
}

options(stringsAsFactors = FALSE)
#source("stan_helpers.R")

option_list = list(
    make_option(
        c("-d", "--in_direc"), type="character",
        help="Name of directory containing all RData files will loo results",
    )
)
 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

in_direc = opt$in_direc

print("Reading loo objects")
fnames = list.files(in_direc, "loo\\_result\\_extdist\\_.*", full.names=TRUE)
print("Found the following files:")
print(fnames)

fname_split_mat = str_extract_all(
    fnames,
    "(?<=extdist\\_)\\d+|(?<=extflen\\_)\\d+|(?<=inpflen\\_)\\d+|(?<=inpdist\\_)\\d+",
    simplify=TRUE
)
ext_distances = as.integer(fname_split_mat[,1])
ext_f_lengths = as.integer(fname_split_mat[,2])
inp_distances = as.integer(fname_split_mat[,3])
inp_f_lengths = as.integer(fname_split_mat[,4])

loo_tib = tibble(
    fname=fnames,
    ext_distance=ext_distances,
    inp_distance=inp_distances,
    ext_frag_len=ext_f_lengths,
    inp_frag_len=inp_f_lengths
) %>% mutate(loo_res = purrr::map(fname, read_loo_result))

comparison = loo_compare(loo_tib$loo_res)

row_model_numbers = as.integer(
    str_split(rownames(comparison), "model", simplify=TRUE)[,2]
)
row_sort_order = order(row_model_numbers)
effect_sizes = comparison[,1] / comparison[,2]
effect_sizes[comparison[,1] == 0] = 0
res_tib = loo_tib %>%
    mutate(
        loo_effect_size = effect_sizes[row_sort_order],
        loo_val=comparison[row_sort_order,1]
    ) %>% select(ext_distance, inp_distance, ext_frag_len, inp_frag_len, loo_effect_size, loo_val)
print(res_tib %>% arrange(desc(loo_effect_size)) %>% head)
print(res_tib %>% arrange(desc(loo_effect_size)) %>% tail)
write_csv(res_tib, file=paste(in_direc, "loo_comparison.csv", sep="/"))

min_eff_size = -5
#min_eff_size = min(res_tib$loo_effect_size)
grad_lims = c(min_eff_size, 0)

for (ext_len in unique(res_tib$ext_frag_len)) {
    for (distance in unique(res_tib$ext_distance)) {
        res_tib %>% filter(
            ext_frag_len == ext_len,
            ext_distance == distance
        ) %>% mutate(
            loo_effect_size = ifelse(loo_effect_size < min_eff_size, -5, loo_effect_size)
        ) %>% ggplot(aes(x=inp_frag_len, y=inp_distance, fill=loo_effect_size)) + 
            geom_tile() +
            scale_fill_gradient(
                low="white",
                high="red",
                name="elpd_loo effect size (higher is better)",
                limits=grad_lims
            ) +
            theme_classic() +
            labs(
                y="Distance between input DNA parameters (bp)",
                x="Assumed mean input DNA fragment length (bp)"
            )
        plot_file = paste0("loo_comparison_ext_dist_", distance, "_ext_len_", ext_len, ".png")
        ggsave(filename=paste(in_direc, plot_file, sep="/"))
    }
}
