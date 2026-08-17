// HiCBZIP-N(GS) spatial/GM model without covariates.
// This is the cleaned release copy of BHZIP_GM_nocov_250909.stan.
data {
    int<lower=1> N;                // Number of locus pairs
    int<lower=1> K;                // Number of cells
    vector[K] lambda;              // Cell-specific coverage scaling
    vector[N] b_i;                 // Prior mean for logit(pi)
    array[N, K] int<lower=0> Y;    // Observed counts
    vector[N] theta;               // Prior mean for log(mu)
    real<lower=0> tau;             // Prior SD for log(mu)
}

parameters {
    real<lower=1e-6> sigma_b;
    real<lower=1e-6> sigma_a;
    vector[N] a_i;
    matrix[N, K] logit_pi;
    matrix[N, K] log_mu;
}

transformed parameters {
    matrix[N, K] pi;
    matrix[N, K] mu;

    for (i in 1:N) {
        for (k in 1:K) {
            pi[i, k] = inv_logit(logit_pi[i, k]);
            mu[i, k] = exp(log_mu[i, k]);
        }
    }
}

model {
    a_i ~ normal(theta, tau);
    sigma_b ~ normal(0, 2);
    sigma_a ~ normal(0, 2);

    for (i in 1:N) {
        for (k in 1:K) {
            logit_pi[i, k] ~ normal(b_i[i], sigma_b);
            log_mu[i, k] ~ normal(a_i[i], sigma_a);
        }
    }

    for (i in 1:N) {
        for (k in 1:K) {
            if (Y[i, k] == 0) {
                target += log_sum_exp(
                    bernoulli_lpmf(1 | pi[i, k]),
                    bernoulli_lpmf(0 | pi[i, k]) + poisson_lpmf(0 | lambda[k] * mu[i, k])
                );
            } else {
                target += bernoulli_lpmf(0 | pi[i, k]) +
                          poisson_lpmf(Y[i, k] | lambda[k] * mu[i, k]);
            }
        }
    }
}

generated quantities {
    array[N, K] int S;
    matrix[N, K] muS;
    for (i in 1:N) {
        for (k in 1:K) {
            real posterior_prob_S1;
            if (Y[i, k] == 0) {
                posterior_prob_S1 = pi[i, k] /
                    (pi[i, k] + (1 - pi[i, k]) * exp(-lambda[k] * mu[i, k]));
            } else {
                posterior_prob_S1 = 0;
            }
            S[i, k] = bernoulli_rng(posterior_prob_S1);
            muS[i, k] = mu[i, k] * (1 - S[i, k]);
        }
    }
}
