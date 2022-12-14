functions {
    // Taken straight from stan code generated by brms
    /* Efficient computation of the horseshoe prior
        * see Appendix C.1 in https://projecteuclid.org/euclid.ejs/1513306866
        * Args:
        *   z: standardized population-level coefficients
        *   lambda: local shrinkage parameters
        *   tau: global shrinkage parameter
        *   c2: slab regularization parameter
        * Returns:
        *   population-level coefficients following the horseshoe prior
    */
    vector horseshoe(vector z, vector lambda, real tau, real c2) {
        int K = rows(z);
        vector[K] lambda2 = square(lambda);
        vector[K] lambda_tilde = sqrt(c2 * lambda2 ./ (c2 + tau^2 * lambda2));
        return z .* lambda_tilde * tau;
    }
}

data {
    int<lower=1> L; // number of positions
    //int<lower=1> a_C; // correlation distance (in positions, not basepairs)
    //int<lower=1> b_C; // correlation distance (in positions, not basepairs)
    int<lower=1> S; // number of samples
    int<lower=1> B; // number of betas
    int<lower=1> A; // number of alphas
    int<lower=1> G; // number of genotypes
    int<lower=1, upper=2> Q; // number of strands
    real alpha_prior;
    array[S] int geno_x; // genotype for each sample
    array[S] int sample_x; // sample_id for each datum
    array[S,Q,L] int Y; // counts
    array[S] real<lower=0> libsize; // exposure term for each sample
    real<lower=0> hs_df; // local df
    real<lower=0> hs_df_global; // global df
    real<lower=0> hs_df_slab; // slab df
    real<lower=0> hs_scale_global; // global prior scale
    real<lower=0> hs_scale_slab; // slab prior scale
    //int<lower=1> gauss_dist;
    int<lower=1> a_sub_L;
    int<lower=1> b_sub_L;
    int<lower=1> b_num_non_zero;
    vector<lower=0,upper=1>[b_num_non_zero] b_weights_vals;
    array[b_num_non_zero] int b_col_accessor;
    array[L+1] int b_row_non_zero_number;
    int<lower=1> a_num_non_zero;
    vector<lower=0,upper=1>[a_num_non_zero] a_weights_vals;
    array[a_num_non_zero] int a_col_accessor;
    array[L+1] int a_row_non_zero_number;
    int<lower=0, upper=1> gather_log_lik;
}

transformed data {
    array[S] real cent_loglibsize;
    int N = L*Q*S;
    int num_hs = b_sub_L*B*Q;
    {
        array[S] real log_libsize;
        real mean_loglib;
        for (s in 1:S) {
            log_libsize[s] = log(libsize[s]);
        }
        mean_loglib = mean(log_libsize);
        for (s in 1:S) {
            cent_loglibsize[s] = log_libsize[s] - mean_loglib;
        }
    }
}

parameters {
    array[2] real log_prec; // stratify global precision inference by extracted vs input
    array[G,Q] vector[a_sub_L] sub_Alpha; // one intercept for each genotype/sub-position
    vector[B] Gamma; // an intercept to offset hbd samples by

    // local parameters for horseshoe prior
    array[B,Q] vector[b_sub_L] zbeta;
    array[B,Q] vector<lower=0>[b_sub_L] hs_local;
    // horseshoe shrinkage parameters 
    real<lower=0> hs_global; 
    real<lower=0> hs_slab; 
    real<lower=0> shape; 
}

transformed parameters {
    array[2] real<lower=0> prec = exp(log_prec);
    real lprior = 0; // prior contributions to log posterior
    array[B,Q] vector[L] tmp_Beta;
    array[B,Q] vector[L] Beta;
    array[G,Q] vector[L] Alpha; // one intercept for each genotype/position combination
    array[S,Q] vector[L] Y_hat;

    lprior += student_t_lpdf(Gamma | 3, -5, 5);

    {
        vector[b_sub_L] sub_Beta;
        for (b in 1:B) {
            for (q in 1:Q) {
                sub_Beta = horseshoe(
                    zbeta[b,q],
                    hs_local[b,q],
                    hs_global,
                    hs_scale_slab^2 * hs_slab
                );
                //print(sub_Beta[1:5])
                tmp_Beta[b,q] = csr_matrix_times_vector(
                    L,
                    b_sub_L,
                    b_weights_vals,
                    b_col_accessor,
                    b_row_non_zero_number,
                    sub_Beta
                );
                Beta[b,q] = tmp_Beta[b,q] + Gamma[b];
            }
        }
    }

    for (g in 1:G) {
        for (q in 1:Q) {
            lprior += normal_lpdf(sub_Alpha[g,q] | alpha_prior, 4);
            Alpha[g,q] = csr_matrix_times_vector(
                L,
                a_sub_L,
                a_weights_vals,
                a_col_accessor,
                a_row_non_zero_number,
                sub_Alpha[g,q]
            );
        }
    }

    lprior += student_t_lpdf(hs_global | hs_df_global, 0, hs_scale_global)
        - 1 * log(0.5);
    lprior += inv_gamma_lpdf(hs_slab | 0.5 * hs_df_slab, 0.5 * hs_df_slab);
    lprior += gamma_lpdf(shape | 0.01, 0.01);
    lprior += normal_lpdf(log_prec | 2, 1);

    {
        int genotype;
        int sample_type;
        for (s in 1:S) {
            genotype = geno_x[s];
            sample_type = sample_x[s];
            for (q in 1:Q) {
                Y_hat[s,q] = Alpha[genotype,q]
                    + sample_type * Beta[genotype,q]
                    + cent_loglibsize[s];
            }
        }
    }
}

model {
    int grainsize = 1;
    int sample_type;

    target += lprior;
    for (b in 1:B) {
        for (q in 1:Q) {
            target += std_normal_lpdf(zbeta[b,q]);
            target += student_t_lpdf(hs_local[b,q] | hs_df, 0, 1)
              - num_hs * log(0.5);
        }
    }

    for (s in 1:S) {
        sample_type = sample_x[s];
        for (q in 1:Q) {
            target += neg_binomial_2_log_lpmf(Y[s,q] | Y_hat[s,q], prec[sample_type+1]);
        }
    }
}

generated quantities {
    vector[N] log_lik;
    array[S,Q] vector[L] post_pred;
    {
        int n = 1;
        int sample_type;
        if (gather_log_lik) {
            for (s in 1:S) {
                sample_type = sample_x[s];
                for (q in 1:Q) {
                    for (l in 1:L) {
                        post_pred[s,q][l] = neg_binomial_2_rng(exp(Y_hat[s,q][l]), prec[sample_type+1]);
                        log_lik[n] = neg_binomial_2_log_lpmf(Y[s,q,l] | Y_hat[s,q][l], prec[sample_type+1]);
                        n += 1;
                    }
                }
            }
        }
    }
}

