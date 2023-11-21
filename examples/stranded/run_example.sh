#!/bin/bash

singularity exec -B $(pwd) /corexfs/schroedj/appliances/Enricherator/dev/enricherator_0.0.1.sif \
    Rscript /src/R/cmdstan_fit_enrichment_model.R \
        --info stan_sample_info.txt \
        --stan_file /src/stan/ab_smoothed_stranded_enrichment.stan \
        --ignore_ctgs P2918_rnadna_spikein,NC_011916.1 \
        --norm_method libsize \
        --fit_file "example_results/fit.RData" \
        --data_file "example_results/data.RData" \
        --ext_subsample_dist 30 \
        --ext_fragment_length 60 \
        --input_subsample_dist 60 \
        --input_fragment_length 120 \
        --libsize_key "tmm_size_factors" \
        --draws_direc "example_results/draws" \
        --cores 8 \
        > "example_results/fit.log" \
        2> "example_results/fit.err"
