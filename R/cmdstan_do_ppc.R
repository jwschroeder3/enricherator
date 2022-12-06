#!/usr/bin/env Rscript

library(tidyverse)
library(tidybayes)
library(optparse)
library(cmdstanr)
library(bayesplot)

options(stringsAsFactors = FALSE)
#source("stan_helpers.R")

option_list = list(
    make_option(
        c("-f", "--fit_file"), type="character",
        help="Name of RData file containing stanfit object",
    ),
    make_option(
        c("-d", "--data_file"), type="character",
        help="Name of RData file containing stanfit data",
    ),
    make_option(
        c("-p", "--ppc_file"), type="character",
        help="Name of file to save posterior predictive check result",
    )
)
 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

fit_file = opt$fit_file

print("Reading stanfit object")
load(fit_file)
load(opt$data_file)
Y = stan_list[["Y"]]
samp_ids = stan_list[["sample_x"]]
post_pred_df = gather_draws(fit, post_pred[sample,strand,position])
ndraws = max(post_pred_df$.draw)
#Y_hat = exp(extract(fit, pars="Y_hat")$Y_hat)
y_dims = dim(Y)
y_vec_len = y_dims[1] * y_dims[2] * y_dims[3]

pred_mat = matrix(0, nrow=50, ncol=y_vec_len)
Y_vec = c()
pred_vec = c()
group_vec = c()
for (i in 1:y_dims[1]) {
    genidx = stan_list[["geno_x"]][i]
    genoname = unique(stan_list[["info"]][stan_list[["info"]]$geno_x == genidx,]$genotype) 
    samptype = stan_list[["sample_x"]][i]
    if (samptype == 0) {
        sampname = "input"
    } else {
        sampname = "hbd"
    }
    for (j in 1:y_dims[2]) {
        Y_vec = c(Y_vec, Y[i,j,])
        group_vec = c(group_vec, rep(paste0(i, "; ", genoname, "; ", sampname), y_dims[3]))
        fity_count = 1
        start = 1 + 2*y_dims[3]*(i-1) + y_dims[3]*(j-1)
        end = start + y_dims[3] - 1
        for (l in as.integer(seq(1,ndraws,length.out=50))) {
            this_vec = post_pred_df %>%
                filter(.draw==l, sample==i, strand==j) %>% 
                .$.value
            pred_mat[fity_count,start:end] = this_vec
            fity_count = fity_count + 1
        }
    }
}


p = ppc_dens_overlay_grouped(Y_vec, pred_mat, group_vec) +
    coord_cartesian(xlim=c(0,250))
ggsave(filename=opt$ppc_file, plot=p)

