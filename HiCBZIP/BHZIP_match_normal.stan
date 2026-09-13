// correct 11/08: lambda now is a vector
data {
  int<lower=1> N;         // number of observations
  array[N] int<lower=0> Y;      // observed counts
  vector<lower=0>[N] lambda;            // per-observation lambda
  real a_norm;            // mean for log(mu)
  real<lower=0> sigma_mu; // sd for log(mu)
  real b_norm;            // mean for logit(pi)
  real<lower=0> sigma_pi; // sd for logit(pi)
}
parameters {
  vector[N] logmu;       // latent log(mu)
  vector[N] logitpi;     // latent logit(pi)
}
transformed parameters {
  vector<lower=0>[N] mu = exp(logmu);
  vector<lower=0,upper=1>[N] pi = inv_logit(logitpi);
}
model {
  // Priors on the latent parameters:
  logmu ~ normal(a_norm, sigma_mu);
  logitpi ~ normal(b_norm, sigma_pi);
  
  // Likelihood: marginalize out S
  for (n in 1:N) {
    real rate = lambda[n] * mu[n];
    if (Y[n] == 0) {
      target += log( pi[n] + (1 - pi[n]) * exp(-rate) );
    } else {
      target += log1m(pi[n]) + poisson_lpmf(Y[n] | rate);
    }
  }
}
generated quantities {
  // Back out S and compute imputed signal mu*(1-S)
  array[N] int<lower=0,upper=1> S;
  vector[N] mu_tilde;
  for (n in 1:N) {
    real rate = lambda[n] * mu[n];
    if (Y[n] > 0)
      S[n] = 0;
    else {
      real pS1 = pi[n] / (pi[n] + (1 - pi[n]) * exp(-rate));
      S[n] = bernoulli_rng(pS1);
    }
    mu_tilde[n] = mu[n] * (1 - S[n]);
  }
}
