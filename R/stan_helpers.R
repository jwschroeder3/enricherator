library(tidyverse)
#library(rstan)

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

read_file = function(fname, ignores, debug=FALSE) {
    print(paste0("Reading data from ", fname))
    tmp_tib = read_delim(
        file=fname,
        col_names=c("seqname","start","end","score","strand"),
        col_types="ciidc",
        delim="\t"
    ) %>% filter(
        !(seqname %in% ignores)
    ) %>% mutate(
        score = as.integer(score)
    )

    if (debug) {
        tmp_tib = tmp_tib %>%
            filter(
                ((seqname == "CP006881.1" & start > 1180500 & end < 1185500)
                | seqname == "NC_011916.1")
            )
    }

    print(tmp_tib %>% head)
    if (all(is.na(tmp_tib$strand))) {
        tmp_tib$strand = "*"
    }
    return(tmp_tib)
}

set_up_mapper = function(vec) {
    mapper = list()
    for (i in 1:length(vec)) {
        mapper[[vec[i]]] = i
    }
    return(mapper)
}

place_val_in_list = function(.val_list, .row, .col, .value) {
    # if the named position does not exist as a key in val_list, add it
    next_idx = .val_list[[.row]][["len"]] + 1
    .val_list[[.row]][["values"]][next_idx] = .value
    .val_list[[.row]][["cols"]][next_idx] = .col
    .val_list[[.row]][["len"]] = next_idx
    if (any(is.na(.val_list[[.row]][["values"]]))) {
        print(.row)
        print(.col)
        print(.value)
        print(.val_list[[.row]][["values"]])
        print(next_idx)
        stop("na in values in place_val_in_list")
    }
    if (any(is.na(.val_list[[.row]][["cols"]]))) {
        print(.row)
        print(.col)
        print(.val_list[[.row]][["cols"]])
        stop("na in cols in place_val_in_list")
    }
    return(.val_list)
}

initialize_val_list = function(L) {
    val_list = vector(mode="list", L)
    for (i in 1:length(val_list)) {
        val_list[[i]] = list("values"=c(), "cols"=c(), "len"=0)
    }
    return(val_list)
}

make_weights_array = function(L, C, gauss_dist) {
    print(paste0("Constructing sparse weights array for correlation length of ", C, " positions"))
    #sub_L = as.integer(L/C)
    sub_L = L

    x = -5000:5000
    dist = dnorm(x, 0, C/3)
    dist2 = dist[dist>=0.01]
    gauss_kern = dist2 / sum(dist2)
    gauss_dist = length(gauss_kern)

    val_list = initialize_val_list(L)

    report = 50
    bottom_idx = gauss_dist
    for (sl in 1:sub_L) {
        if (sl %% report == 0) {
            print(paste0("Placing data for subset Beta ", sl, " of ", sub_L, " into val_list"))
        }
        top_idx = 1;
        top = sl - C;
        bottom = sl + C;

        #print("gauss_dist")
        #print(gauss_dist)
        #print("top_idx")
        #print(top_idx)
        #print("center")
        #print(sl)
        #print("top")
        #print(top)
        #print("bottom")
        #print(bottom)
        #print("total positions")
        #print(L)
        
        if (top < 1) {
            top_idx = top * -1 + 2;
            top = 1;
        }

        if (bottom > L) {
            bottom_clip = bottom - L;
            bottom_idx = gauss_dist - bottom_clip
            bottom = L;
        }

        #print("adjusted_top_idx")
        #print(top_idx)
        #print("adjusted_bottom_idx")
        #print(bottom_idx)

        rows = top:bottom
        values = gauss_kern[top_idx:bottom_idx]
        print(rows)
        if (any(is.na(values))) {
            print(sl)
            print(values)
            stop("stopped after grabbing values from gauss_kern")
        }

        for (i in 1:length(rows)) {
            if (rows[i] == 7) {
                print(values)
            }
            val_list = place_val_in_list(val_list, rows[i], sl, values[i])
        }
    }

    # fill in final spots with 1.0
    #sl = sub_L
    #remaining = L-bottom;
    #if (remaining > 0) {
    #    rows = (bottom+1):L
    #    for (l in rows) {
    #        val_list = place_val_in_list(val_list, l, sl, 1.0)
    #    }
    #}

    # normalize each position values to sum to 1.0
    print("Normalizing each position's weights to sum to 1.0.")
    for (l in 1:L) {
        these_vals = val_list[[l]][["values"]]
        if (any(is.na(these_vals))) {
            print(l)
            print(these_vals)
            stop("stopped in normalizing vals")
        }
        val_list[[l]][["values"]] = these_vals / sum(these_vals)
    }

    #return(val_list)

    ######################################################################
    ## w = values ########################################################
    ## v = column accessor ###############################################
    ## u = number non-zero ###############################################
    ######################################################################

    print("Creating sparse matrix values from compiled list")
    sparse_parts = get_sparse_vals(val_list)
    #save(L, C, sub_L, val_list, remaining, bottom, sparse_parts, file="debug_weights.RData")
    #print(sparse_parts[["w"]][which(is.na(sparse_parts[["w"]]))])
    return(sparse_parts)
}

get_sparse_vals = function(.val_list) {
    w = c()
    v = c()
    u = c(1)
    len = 1
    for (i in 1:length(.val_list)) {
        l_info = .val_list[[i]]
        w = c(w, l_info[["values"]])
        if (any(is.na(w))) {
            print(i)
            print(w)
            stop()
        }
        v = c(v, l_info[["cols"]])
        len = len + l_info[["len"]]
        u = c(u, len)
    }
    return(list("w"=w, "v"=v, "u"=u))
}

geom_mean = function(x) {
    return(exp(mean(log(x))))
}

get_size_factors = function(mat, summary_fn=get_trimmed_mean) {
    referenced = matrix(0, dim(mat)[1], dim(mat)[2])
    for (i in 1:dim(mat)[1]) {
        gm_i = geom_mean(mat[i,])
        for (j in 1:dim(mat)[2]) {
            referenced[i,j] = mat[i,j]/gm_i
        }
    }
    size_factors = NA
    for (j in 1:dim(mat)[2]) {
        size_factors[j] = summary_fn(referenced[,j])
    }
    return(size_factors)
}

get_trimmed_mean = function(x) {
    x = sort(x)
    len = length(x)
    start = as.integer(len/40)
    end = as.integer(len-start)
    return(mean(x[start:end]))
}

insert_data = function(data, spike_name, spikein_rel_abund) {
    data_list = list()
    data_list[["alpha_prior"]] = log(1/spikein_rel_abund)
    pos_num = nrow(data$data[[1]] %>% filter(seqname != spike_name))
    data_list[["position_mapper"]] = tibble(
        seqnames=data$data[[1]] %>% filter(seqname != spike_name) %>% .$seqname,
        start=data$data[[1]] %>% filter(seqname != spike_name) %>% .$start,
        end=data$data[[1]] %>% filter(seqname != spike_name) %>% .$end,
        position_row=1:pos_num
    )

    #print("norm_df line 252")
    #print(norm_df)
    if (spike_name != "") {
        norm_df = data %>% 
            select(-c(strand,file,norm_factor)) %>%
            unnest(data)
        norm_factors = data$norm_factor[!duplicated(data$sample_id)]
        norm_df = norm_df %>% filter(seqname == spike_name)
    } else {
        norm_df = data %>% 
            select(-c(strand,file)) %>%
            unnest(data)
        norm_factors = 1
    }
    #print("norm_df line 260")
    #print(norm_df)
    # works, whether stranded or not
    norm_df = norm_df %>%
        select(sample_id, seqname, start, strand, score) %>%
        group_by(sample_id, seqname, start) %>%
        summarize(score = sum(score))
        #spread(strand, score) 
    #print("norm_df line 268")
    #print(norm_df)
    #norm_df = norm_df %>%
    #    group_by(sample_id,locus) %>%
    #    mutate(score=sum(`-`,`+`))
        
    norm_df = norm_df %>%
        select(sample_id,seqname,start,score) %>%
        ungroup() %>%
        spread(sample_id,score) %>%
        select(-c(seqname,start))
    #print("norm_df line 279")
    #print(norm_df)
    #sample_ids = names(spikein_df)
    norm_mat = as.matrix(norm_df)
    #tmm_deseq_size_factors = get_size_factors(spike_mat, summary_fn=get_trimmed_mean)
    #med_deseq_size_factors = get_size_factors(spike_mat, summary_fn=median)
    #mean_size_factors = apply(spike_mat, FUN=mean, MAR=2)
    tmm_size_factors = apply(norm_mat, FUN=get_trimmed_mean, MAR=2)
    #print("size factor shape")
    #print(dim(tmm_size_factors))
    #print("norm factors")
    #print(norm_factors)
    #stop()
    
    #stan_list[["mean_size_factors"]] = mean_size_factors
    data_list[["tmm_size_factors"]] = tmm_size_factors
    data_list[["spikein_norm_factors"]] = norm_factors
    #stan_list[["tmm_deseq_size_factors"]] = tmm_deseq_size_factors
    #stan_list[["med_deseq_size_factors"]] = med_deseq_size_factors

    return(list(pos_num, data_list))
}


prep_stan_data = function(data_df, norm_method, spikein, spikein_rel_abund=0.05, input_subsample_dist_bp=50, input_frag_len_bp=125, ext_subsample_dist_bp=50, ext_frag_len_bp=50, log_lik=FALSE) {
    # this sort step is VERY USEFUL for arranging data later on to keep size factors
    # associated with correct samples
    data_df = data_df %>% arrange(sample_id)

    if (norm_method == "libsize") {
        spikein=""
        res = insert_data(
            data_df,
            spike_name=spikein,
            spikein_rel_abund=1
        )
    } else if (norm_method == "spikein") {
        res = insert_data(
            data_df,
            spike_name=spikein,
            spikein_rel_abund=spikein_rel_abund
        )
    } else {
        stop(paste0("--norm_method must be either 'libsize' or 'spikein', but it is currently set to ", norm_method, ". Exiting now."))
    }

    pos_num = res[[1]]
    stan_list = res[[2]]
    num_geno = length(unique(data_df$genotype))
    levels = unique(data_df$sample_id)
    strand_names = unique(data_df$strand)
    samp_num = length(levels)
    strand_num = length(strand_names)
    info_df = data_df %>% select(-data) %>% filter(!duplicated(sample_id))
    stan_list[["L"]] = pos_num
    stan_list[["S"]] = samp_num
    stan_list[["B"]] = num_geno
    stan_list[["A"]] = num_geno
    stan_list[["G"]] = num_geno
    stan_list[["Q"]] = strand_num
    stan_list[["geno_x"]] = info_df$geno_x
    stan_list[["sample_x"]] = info_df$sample_x
    stan_list[["strand_x"]] = info_df$strand_x
    stan_list[["info"]] = data_df %>% select(-data)
    stan_list[["gather_log_lik"]] = log_lik
    data_arr = base::array(0, dim=c(samp_num,strand_num,pos_num))
    for (s in 1:samp_num) {
        level = levels[s]
        for (q in 1:strand_num) {
            strand_name = strand_names[q]
            which_row = which(data_df$sample_id == level & data_df$strand == strand_name)
            sq_data  = unlist(data_df$data[[which_row]] %>%
                filter(seqname != spikein) %>%
                select(score))
            data_arr[s,q,] = sq_data
        }
    }
    stan_list[["Y"]] = data_arr

    resolution = stan_list[["position_mapper"]]$start[2] - stan_list[["position_mapper"]]$start[1]

    # round to nearest multiple of "resolution"
    input_subsample_dist_bp = round(input_subsample_dist_bp/resolution)*resolution
    input_frag_len_bp = round(input_frag_len_bp/resolution)*resolution
    ext_subsample_dist_bp = round(ext_subsample_dist_bp/resolution)*resolution
    ext_frag_len_bp = round(ext_frag_len_bp/resolution)*resolution

    a_C = floor(input_subsample_dist_bp / resolution)
    b_C = floor(ext_subsample_dist_bp / resolution)
    b_K = ext_frag_len_bp / resolution
    a_K = input_frag_len_bp / resolution
    # enforce that C must be odd
    if (a_C %% 2 == 0) {
        a_C = a_C+1
    }
    if (b_C %% 2 == 0) {
        b_C = b_C+1
    }
    stan_list[["a_C"]] = a_C
    a_sub_L = as.integer(pos_num/a_C)
    stan_list[["a_sub_L"]] = a_sub_L
    stan_list[["b_C"]] = b_C
    b_sub_L = as.integer(pos_num/b_C)
    stan_list[["b_sub_L"]] = b_sub_L
    #stan_list[["sub_L"]] = pos_num

    by_ctg_df = data_df$data[[1]] %>% filter(seqname != spikein) %>%
        group_by(seqname) %>%
        summarize(ctg_posnum = n())
    ctg_end_vec = by_ctg_df$ctg_posnum
    #print("ctg_end_vec")
    #print(ctg_end_vec)
    # start with zero-indexted end
    accum_ctg_ends = ctg_end_vec[1] - 1
    for (i in seq_along(ctg_end_vec)) {
        if (i == 1) {
            next
        }
        # accumulate lengths for each ctg_end
        accum_ctg_ends[i] = accum_ctg_ends[i-1] + ctg_end_vec[i]
    }
    # make ends into comma-separated list as required by make_sparse_matrix
    ctg_ends = paste(accum_ctg_ends, collapse=",")

    w_out = paste("w", b_K, b_sub_L, "vals.txt", sep="_")
    v_out = paste("v", b_K, b_sub_L, "vals.txt", sep="_")
    u_out = paste("u", b_K, b_sub_L, "vals.txt", sep="_")

    print("Building sparse matrix of weights for beta smoothing")
    print(paste0("Using beta subsampling distance of ", ext_subsample_dist_bp, " to perform smoothing."))
    print(paste0("Using ", ctg_ends, " as list of 0-indexed contig end positions, at resolution ", resolution))

    exec_f_path = get_exec_file()
    exec_direc = dirname(exec_f_path)
    bin_path = file.path(exec_direc, "../bin/make_sparse_matrix")
    
    #if (!file.exists(w_out)) {
    print("Running command:")
    print(paste("make_sparse_matrix", pos_num, b_sub_L, b_C, b_K, ctg_ends, w_out, v_out, u_out, sep=" "))
    res = system2(
        bin_path,
        c(
            as.character(pos_num),
            as.character(b_sub_L),
            as.character(b_C),
            as.character(b_K),
            ctg_ends,
            w_out,
            v_out,
            u_out
        ),
        env="RUST_BACKTRACE=1",
        stderr="",
        stdout=""
    )
    if (res != 0) stop("Error in sparse matrix creation. Check err and log files.")
    #} else {
    #    print(paste("File", w_out, "already exists. Reading it.", sep=" "))
    #}
    w = unlist(read_table(w_out, col_names=F), use.names=F)
    v = unlist(read_table(v_out, col_names=F), use.names=F)
    u = unlist(read_table(u_out, col_names=F), use.names=F)
    #stan_list[["gauss_dist"]] = C*2+1
    #sparse_weights_info = make_weights_array(pos_num, C, stan_list[["gauss_dist"]])
    #save(sparse_weights_info, file="debug_weights.RData")
    stan_list[["b_weights_vals"]] = w
    stan_list[["b_num_non_zero"]] = length(stan_list[["b_weights_vals"]])
    stan_list[["b_col_accessor"]] = v
    stan_list[["b_row_non_zero_number"]] = u

    w_out = paste("w", a_K, a_sub_L, "vals.txt", sep="_")
    v_out = paste("v", a_K, a_sub_L, "vals.txt", sep="_")
    u_out = paste("u", a_K, a_sub_L, "vals.txt", sep="_")

    print("Building sparse matrix of weights for alpha smoothing")
    print(paste0("Using alpha subsampling distance of ", input_subsample_dist_bp, " to perform smoothing."))
    print(paste0("Using ", ctg_ends, " as list of 0-indexed contig end positions, at resolution ", resolution))
    #if (!file.exists(w_out)) {

    print("Running command:")
    print(paste("make_sparse_matrix ", pos_num, " ", a_sub_L, " ", a_C, " ", a_K, " ", ctg_ends, w_out, v_out, u_out, sep=" "))

    res = system2(
        bin_path,
        c(
            as.character(pos_num),
            as.character(a_sub_L),
            as.character(a_C),
            as.character(a_K),
            ctg_ends,
            w_out,
            v_out,
            u_out
        ),
        env="RUST_BACKTRACE=1",
        stderr="",
        stdout=""
    )
    if (res != 0) {
        stop("Error in sparse matrix creation. Check err and log files.")
    }
    #} else {
    #    print(paste("File", w_out, "already exists. Reading it.", sep=" "))
    #}
    w = unlist(read_table(w_out, col_names=F), use.names=F)
    v = unlist(read_table(v_out, col_names=F), use.names=F)
    u = unlist(read_table(u_out, col_names=F), use.names=F)

    stan_list[["a_weights_vals"]] = w
    stan_list[["a_num_non_zero"]] = length(stan_list[["a_weights_vals"]])
    stan_list[["a_col_accessor"]] = v
    stan_list[["a_row_non_zero_number"]] = u

    stan_list[["hs_df"]] = 1
    stan_list[["hs_scale_global"]] = 1
    stan_list[["hs_df_global"]] = 1
    stan_list[["hs_scale_slab"]] = 2
    stan_list[["hs_df_slab"]] = 4
    return(stan_list)
}

quant = function(vec, interval) {
    quantile(vec, probs=quant_vec)
}

summarize_rstan_inference = function(param_mat, stan_data, quant_vec, out_direc) {

    print("Calculating quantiles of interest from samples")
    quantiles_arr = apply(
        param_mat,
        FUN=quantile,
        MAR=2:4,
        probs=quant_vec
    )
    mean_arr = apply(param_mat, FUN=mean, MAR=2:4)

    print("Writing quantiles of interest to files")
    info_df = stan_data[["info"]]
    genotype_ids = sort(unique(stan_data[["geno_x"]]))
    for (i in 1:length(quant_vec)) {
        quantile = quant_vec[i]
        if (i == 1) {
            quant = "lower"
        } else if (i == 2) {
            quant = "median"
        } else if (i == 3) {
            quant = "upper"
        } else {
            stop("Only three quantiles can be returned.")
        }
        values = quantiles_arr[i,,,]/log(2)
        for (geno_id in genotype_ids) {
            geno_vals = values[geno_id,,]
            genotype = info_df$genotype[info_df$geno_x == geno_id][1]
            for (k in 1:stan_data[["Q"]]) {
                strand_vals = geno_vals[k,]
                strand = info_df$strand[info_df$strand_x == k][1]
                this_df = tibble(
                    seqnames=stan_data[["position_mapper"]]$seqnames,
                    start=stan_data[["position_mapper"]]$start,
                    end=stan_data[["position_mapper"]]$end,
                    score=strand_vals,
                    strand = ifelse(strand == "plus", "+", ifelse(strand == "minus", "-", ifelse(strand=="both", "*", strand)))
                )
                out_fname = paste0(out_direc, "/", genotype, "_", strand, "_strand_beta_", quant, ".bedgraph")
                print(paste0("Writing ", quant, " beta estimates to ", out_fname))
                readr::write_tsv(
                    this_df,
                    path=out_fname,
                    col_names=FALSE
                )
            }
        }
    }

    #m_dims = dim(mean_arr)
    #mean_df = tibble(
    #    geno_x=rep(1:m_dims[1],m_dims[3]*m_dims[2]),
    #    strand_x=rep(1:m_dims[2], each=m_dims[1], m_dims[3]),
    #    position_row=rep(1:m_dims[3], each=m_dims[1]*m_dims[2]),
    #    value=as.vector(mean_arr),
    #    quant="mean"
    #)

    #quantiles = as_tibble(quantiles_arr)
    #cols = ncol(quantiles)
    #quantiles$quant = as.character(quant_vec * 100)

    #param_df = gather(quantiles, key="param", value="value", 1:cols) %>%
    #    separate(
    #        param,
    #        into=c("geno_x","strand_x","position_row"),
    #        convert=TRUE
    #    ) %>%
    #    bind_rows(mean_df) %>%
    #    # convert ln(x) to log2(x)
    #    mutate(value = value/log(2)) %>%
    #    left_join(stan_data[["info"]] %>%
    #        filter(sample != "input", rep_id==1) %>%
    #        select(-c(file,sample,rep,sample_x,sample_id)),
    #        by=c("geno_x","strand_x")
    #    ) %>%
    #    select(-c(geno_x, strand_x)) #%>%
    #    #left_join(stan_data[["position_mapper"]])
    #return(param_df)
}

prep_param_df_for_file = function(df, stan_data, .genotype, .strand, .quant) {
    this_df = stan_data[["position_mapper"]] %>%
        select(-position_row) %>%
        bind_cols(
            df %>% filter(
                genotype==.genotype,
                strand==.strand,
                quant==.quant
            ) %>% select(value)
        ) %>%
        mutate(
            strand = ifelse(.strand == "plus", "+", ifelse(.strand == "minus", "-", ifelse(.strand=="both", "*", .strand)))
        )
    return(this_df)
}


write_cmdstan_summaries = function(summary_df, stan_data, out_direc, params, contrasts=NULL, contrast_type=NULL) {

    info_df = stan_data[["info"]]
    if (is.null(contrasts)) {
        genotype_ids = sort(unique(stan_data[["geno_x"]]))
        strand_ids = sort(unique(stan_data[["strand_x"]]))
    } else {
        contrast_vec = str_split(contrasts, ",", simplify=TRUE)[1,]
        if (contrast_type == "genotype") {
            genotype_ids = 1:length(contrast_vec)
            strand_ids = sort(unique(stan_data[["strand_x"]]))
        }
        if (contrast_type == "strand") {
            genotype_ids = sort(unique(stan_data[["geno_x"]]))
            strand_ids = 1:length(contrast_vec)
        }
    }

    print(strand_ids)
    print(genotype_ids)
    print("Writing quantiles of interest to files")
    for (param in params) {
        print(paste0("Param: ", param))
        for (geno_id in genotype_ids) {
            if (is.null(contrasts)) {
                genotype = info_df$genotype[info_df$geno_x == geno_id][1]
            } else {
                if (contrast_type == "genotype") {
                    genotype = contrast_vec[geno_id]
                } else {
                    genotype = info_df$genotype[info_df$geno_x == geno_id][1]
                }
            }
            for (k in strand_ids) {

                .strand = info_df$strand[info_df$strand_x == k][1]

                if (!is.null(contrasts)) {
                    if (contrast_type == "genotype") {
                        .strand = info_df$strand[info_df$strand_x == k][1]
                    } else { if (contrast_type == "strand") {
                        .strand = contrast_vec[k]
                    }}
                }

                print(paste0("Geno_id: ", geno_id))
                this_df = summary_df %>%
                    filter(genotype == geno_id, strand == k, var_name == param) %>%
                    mutate(
                        seqnames=stan_data[["position_mapper"]]$seqnames,
                        start=stan_data[["position_mapper"]]$start,
                        end=stan_data[["position_mapper"]]$end,
                        strand = ifelse(.strand == "plus", "+",
                            ifelse(.strand == "minus", "-",
                            ifelse(.strand=="both", "*", .strand)))
                    )
                for (quant in c("lower", "median", "mean", "upper")) {
                    out_fname = paste0(out_direc, "/", genotype, "_",
                        .strand, "_strand_", param, "_", quant, ".bedgraph")
                    print(paste0("Writing ", quant,
                            " ", param, " estimates to ", out_fname))
                    readr::write_tsv(
                        this_df[,c("seqnames","start","end",quant,"strand")],
                        file=out_fname,
                        col_names=FALSE
                    )
                }
            }
        }
    }
}

gather_vb_estimates = function(fit, stan_data, direc, interval=90, cmdstan=FALSE, params=NULL) {
    tail = (100-interval)/2
    quant_vec = c(tail/100, 0.5, (100-tail)/100)

    if (cmdstan) {
        print("Running command:")
        exec_f_path = get_exec_file()
        exec_direc = dirname(exec_f_path)
        bin_path = file.path(exec_direc, "../bin/parse_cmdstan_csv")
 
        in_file = fit$output_files()
        summary_out_file = file.path(dirname(in_file), "summaries.txt")
        sample_out_file = file.path(dirname(in_file), "samples.csv")
        parse_cmd = paste(
            bin_path, in_file, summary_out_file, sample_out_file, "500", paste(params, collapse=","),
            sep=" "
        )
        print(parse_cmd)
        res = system2(
            bin_path,
            c(in_file, summary_out_file, sample_out_file, "500", paste(params, collapse=",")),
            env="RUST_BACKTRACE=1",
            stderr="",
            stdout=""
        )
        if (res != 0) {
            stop("Error in summarizing cmdstan csv. Check err and log files. Exiting now.")
        }
        param_summaries = read_tsv(summary_out_file)
        write_cmdstan_summaries(param_summaries, stan_data, direc, params)
    } else {
        list_of_draws = rstan::extract(fit, pars=c("Beta","prec"))
        betas = list_of_draws[["Beta"]]
        summarize_rstan_inference(betas, stan_data, quant_vec, direc)
    }
}

gather_em_estimates = function(stanfit, stan_data) {
    par_vals = stanfit$par
    par_names = names(par_vals)
    
    alpha_mask = grepl("Alpha\\[", par_names)
    beta_mask = grepl("Beta\\[", par_names)
    alpha_beta_mask = (alpha_mask | beta_mask)
    betas = par_vals[beta_mask]
    beta_names = names(betas)

    par_names = str_remove(beta_names, "\\]")
    par_names = str_split(par_names, "\\[", simplify=TRUE)
    par_name = par_names[,1]
    geno_loc = str_split(par_names[,2], ",", simplify=TRUE)
    geno_idx = as.integer(geno_loc[,1])
    loc_idx = as.integer(geno_loc[,2])

    geno_names = NA
    ctg_names = NA
    starts = NA
    for (i in 1:length(geno_idx)) {
        i_geno = geno_idx[i]
        i_loc = loc_idx[i]
        genotype = stan_data[["reverse_geno_mapper"]][[i_geno]]
        locus = stan_data[["pos_mapper"]]$position[stan_data[["pos_mapper"]]$pos_idx == i_loc]
        loc = str_split(locus, "_", simplify=TRUE)
        ctg_name = paste(loc[1,1:(ncol(loc)-1)], collapse="_")
        start = as.integer(loc[1,ncol(loc)])
        geno_names[i] = genotype
        ctg_names[i] = ctg_name
        starts[i] = start
    }
    param_df = tibble(param=par_name, quant="estimate", value=betas, genotype=geno_names, seqnames=ctg_names, start=starts)
    return(param_df)
}

