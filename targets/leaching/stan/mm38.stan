functions {

  #include /m70.stan
  #include /m71.stan

}

data {

  // controls
  int<lower = 0, upper = 1> has_leaching;
  int<lower = 0, upper = 1> has_alpha;
  int<lower = 0, upper = 1> prior_only;

  // data
  int<lower=0> N;
  int<lower = 0> N_has_precision;
  vector[N] m;
  vector[N] t;
  vector<lower = 0>[N_has_precision] precision;

  // indices
  int<lower = 0, upper = N> index_measurement_has_precision[N_has_precision];
  int<lower = 0> N_samples;
  int<lower = 0, upper = N_samples> index_samples_to_measurements[N];
  int<lower = 0> N_species;
  int<lower = 0, upper = N_species> index_species_to_samples[N_samples];
  int<lower = 0> N_speciesxstudies;
  int<lower = 0, upper = N_speciesxstudies> index_speciesxstudies_to_samples[N_samples];

  // parameters
  real k_2_p1_p1;
  //real<lower = 0> k_2_p1_p2;
  real phi_2_p2_p1_p1;
  //real<lower = 0> phi_2_p2_p1_p2;
  real k_2_p2_p1;
  //real<lower = 0> k_2_p2_p2;
  real phi_2_p2_p2_p1;
  //real<lower = 0> phi_2_p2_p2_p2;
  real k_2_p3_p1;
  //real<lower = 0> k_2_p3_p2;
  real phi_2_p2_p3_p1;
  //real<lower = 0> phi_2_p2_p3_p2;
  real k_2_p4_p1;
  //real<lower = 0> k_2_p4_p2;
  real phi_2_p2_p4_p1;
  //real<lower = 0> phi_2_p2_p4_p2;

  real<lower = 0> phi_2_p3;

  real<lower = 0> phi_2_p1_p1;
  real<lower = 0> phi_2_p1_p2;
  real<lower = 0> phi_2_p1_p3;

  // leaching
  real l_2_p1_p1;
  //real<lower = 0> l_2_p1_p2;
  real l_2_p2_p1;
  //real<lower = 0> l_2_p2_p2;
  real l_2_p3_p1;
  //real<lower = 0> l_2_p3_p2;
  real l_2_p4_p1;
  //real<lower = 0> l_2_p4_p2;

  // alpha
  real alpha_2_p1_p1;
  real<lower = 0> alpha_2_p1_p2;
  real alpha_2_p2_p1;
  real<lower = 0> alpha_2_p2_p2;
  real alpha_2_p3_p1;
  real<lower = 0> alpha_2_p3_p2;
  real alpha_2_p4_p1;
  real<lower = 0> alpha_2_p4_p2;

  // priors for level standard deviations
  real<lower = 0> k_2_p1_p2_p1;
  real<lower = 0> k_2_p2_p2_p1;
  real<lower = 0> k_2_p3_p2_p1;
  real<lower = 0> k_2_p4_p2_p1;
  real<lower = 0> phi_2_p2_p1_p2_p1;
  real<lower = 0> phi_2_p2_p2_p2_p1;
  real<lower = 0> phi_2_p2_p3_p2_p1;
  real<lower = 0> phi_2_p2_p4_p2_p1;
  real<lower = 0> l_2_p1_p2_p1;
  real<lower = 0> l_2_p2_p2_p1;
  real<lower = 0> l_2_p3_p2_p1;
  real<lower = 0> l_2_p4_p2_p1;
  //real<lower = 0> alpha_2_p2_p2_p1;
  //real<lower = 0> alpha_2_p3_p2_p1;

  // constants
  real<lower = 0> s;

}

parameters {

  real k_2_p1;
  real phi_2_p2_p1;
  vector[N_species] k_2_p2;
  vector[N_species] phi_2_p2_p2;
  vector[N_speciesxstudies] k_2_p3;
  vector[N_speciesxstudies] phi_2_p2_p3;
  vector[N_samples] k_2_p4;
  vector[N_samples] phi_2_p2_p4;
  real<lower = 0> phi_2_p1;
  vector<lower = 0>[N_samples] phi_2;

  // leaching
  real l_2_p1[has_leaching];
  vector[N_species * has_leaching] l_2_p2;
  vector[N_speciesxstudies * has_leaching] l_2_p3;
  vector[N_samples * has_leaching] l_2_p4;

  // alpha
  real alpha_2_p1[has_alpha];
  vector[N_species * has_alpha] alpha_2_p2;
  vector[N_speciesxstudies * has_alpha] alpha_2_p3;
  vector[N_samples * has_alpha] alpha_2_p4;

  // group level standard deviations
  real<lower = 0> phi_2_p2_p1_p2;
  vector<lower = 0>[N_species] phi_2_p2_p2_p2;
  vector<lower = 0>[N_speciesxstudies] phi_2_p2_p3_p2;
  vector<lower = 0>[N_samples] phi_2_p2_p4_p2;
  real<lower = 0> k_2_p1_p2;
  vector<lower = 0>[N_species] k_2_p2_p2;
  vector<lower = 0>[N_speciesxstudies] k_2_p3_p2;
  vector<lower = 0>[N_samples] k_2_p4_p2;
  real<lower = 0> l_2_p1_p2[has_leaching];
  vector<lower = 0>[N_species * has_leaching] l_2_p2_p2;
  vector<lower = 0>[N_speciesxstudies * has_leaching] l_2_p3_p2;
  vector<lower = 0>[N_samples * has_leaching] l_2_p4_p2;
  //real<lower = 0> alpha_2_p2_p2[has_alpha];
  //real<lower = 0> alpha_2_p3_p2[has_alpha];

}

transformed parameters {

  real lprior;
  vector<lower = 0, upper = 1>[N] mu;
  vector<lower = 0>[N] phi;
  vector<lower = 0>[N_samples] k_2;
  vector<lower = 0>[N_samples] phi_2_p2;

  // leaching
  vector<lower = 0, upper = 1>[N_samples * has_leaching] l_2;

  // alpha
  vector<lower = 0>[N_samples * has_alpha] alpha_2;

  // priors

  lprior = normal_lpdf(phi_2_p2_p1 | 0, 1);
  //if(has_leaching) {
  //
  //}
  if(has_alpha) {
    lprior += normal_lpdf(alpha_2_p1 | 0, 1);
  }
  lprior += gamma_lpdf(phi_2_p1 | phi_2_p1_p1, phi_2_p1_p2 * phi_2_p1_p3);

  // group-level standard deviation
  lprior += normal_lpdf(k_2_p1_p2 | 0.0, 1.0);
  lprior += normal_lpdf(k_2_p2_p2 | 0.0, 1.0);
  lprior += normal_lpdf(k_2_p3_p2 | 0.0, 1.0);
  lprior += normal_lpdf(k_2_p4_p2 | 0.0, 1.0);
  lprior += normal_lpdf(phi_2_p2_p1_p2 | 0.0, 1.0);
  lprior += normal_lpdf(phi_2_p2_p2_p2 | 0.0, 1.0);
  lprior += normal_lpdf(phi_2_p2_p3_p2 | 0.0, 1.0);
  lprior += normal_lpdf(phi_2_p2_p4_p2 | 0.0, 1.0);
  if(has_leaching) {
    lprior += normal_lpdf(l_2_p1_p2 | 0.0, 1.0);
    lprior += normal_lpdf(l_2_p2_p2 | 0.0, 1.0);
    lprior += normal_lpdf(l_2_p3_p2 | 0.0, 1.0);
    lprior += normal_lpdf(l_2_p4_p2 | 0.0, 1.0);
  }
  //if(has_alpha) {
  //  lprior += normal_lpdf(alpha_2_p2_p2 | 0.0, 1.0);
  //  lprior += normal_lpdf(alpha_2_p3_p2 | 0.0, 1.0);
  //}

  {

    real k_2_p1_s = k_2_p1 * k_2_p1_p2 * k_2_p1_p2_p1 + k_2_p1_p1;
    real phi_2_p2_p1_s = phi_2_p2_p1 * phi_2_p2_p1_p2 * phi_2_p2_p1_p2_p1 + phi_2_p2_p1_p1;
    vector[N_species] k_2_p2_s = k_2_p2 .* k_2_p2_p2 * k_2_p2_p2_p1 + k_2_p2_p1;
    vector[N_species] phi_2_p2_p2_s = phi_2_p2_p2 .* phi_2_p2_p2_p2 * phi_2_p2_p2_p2_p1 + phi_2_p2_p2_p1;
    vector[N_speciesxstudies] k_2_p3_s = k_2_p3 .* k_2_p3_p2 * k_2_p3_p2_p1 + k_2_p3_p1;
    vector[N_speciesxstudies] phi_2_p2_p3_s = phi_2_p2_p3 .* phi_2_p2_p3_p2 * phi_2_p2_p3_p2_p1 + phi_2_p2_p3_p1;
    vector[N_samples] k_2_p4_s = k_2_p4 .* k_2_p4_p2 * k_2_p4_p2_p1 + k_2_p4_p1;
    vector[N_samples] phi_2_p2_p4_s = phi_2_p2_p4 .* phi_2_p2_p4_p2 * phi_2_p2_p4_p2_p1 + phi_2_p2_p4_p1;

    k_2 = exp(k_2_p1_s + k_2_p2_s[index_species_to_samples] + k_2_p3_s[index_speciesxstudies_to_samples] + k_2_p4_s);
    phi_2_p2 = exp(phi_2_p2_p1_s + phi_2_p2_p2_s[index_species_to_samples] + phi_2_p2_p3_s[index_speciesxstudies_to_samples] + phi_2_p2_p4_s);
    vector[N] k_1 = k_2[index_samples_to_measurements];

    // alpha
    vector[N] alpha_1;
    if(has_alpha) {
      real alpha_2_p1_s = alpha_2_p1[1] * alpha_2_p1_p2 + alpha_2_p1_p1;
      vector[N_species] alpha_2_p2_s = alpha_2_p2 * alpha_2_p2_p2 + alpha_2_p2_p1;
      vector[N_speciesxstudies] alpha_2_p3_s = alpha_2_p3 * alpha_2_p3_p2 + alpha_2_p3_p1;
      vector[N_samples] alpha_2_p4_s = alpha_2_p4 * alpha_2_p4_p2 + alpha_2_p4_p1;
      alpha_2 = 1.0 + exp(alpha_2_p1_s + alpha_2_p2_s[index_species_to_samples] + alpha_2_p3_s[index_speciesxstudies_to_samples] + alpha_2_p4_s);
      alpha_1 = alpha_2[index_samples_to_measurements];
    }

    // leaching
    vector[N] l_1 = rep_vector(0.0, N);
    if(has_leaching) {
      real l_2_p1_s = l_2_p1[1] * l_2_p1_p2[1] * l_2_p1_p2_p1 + l_2_p1_p1;
      vector[N_species] l_2_p2_s = l_2_p2 .* l_2_p2_p2 * l_2_p2_p2_p1 + l_2_p2_p1;
      vector[N_speciesxstudies] l_2_p3_s = l_2_p3 .* l_2_p3_p2 * l_2_p3_p2_p1 + l_2_p3_p1;
      vector[N_samples] l_2_p4_s = l_2_p4 .* l_2_p4_p2 * l_2_p4_p2_p1 + l_2_p4_p1;
      l_2 = inv_logit(l_2_p1_s + l_2_p2_s[index_species_to_samples] + l_2_p3_s[index_speciesxstudies_to_samples] + l_2_p4_s);
      l_1 = l_2[index_samples_to_measurements];
    }

    if(! has_alpha) {
      mu = m71(t, k_1, l_1, s);
    } else {
      mu = m70(t, k_1, l_1, alpha_1, s);
    }

    phi = phi_2[index_samples_to_measurements] * phi_2_p3;

  }

}

model {

  // priors
  target += lprior;

  // herarchical component
  if(has_leaching) {
    target += normal_lpdf(l_2_p1 | 0, 1);
    target += normal_lpdf(l_2_p2 | 0, 1);
    target += normal_lpdf(l_2_p3 | 0, 1);
    target += normal_lpdf(l_2_p4 | 0, 1);
  }
  if(has_alpha) {
    target += normal_lpdf(alpha_2_p2 | 0, 1);
    target += normal_lpdf(alpha_2_p3 | 0, 1);
    target += normal_lpdf(alpha_2_p4 | 0, 1);
  }

  target += normal_lpdf(k_2_p1 | 0, 1);
  target += normal_lpdf(k_2_p2 | 0, 1);
  target += normal_lpdf(k_2_p3 | 0, 1);
  target += normal_lpdf(k_2_p4 | 0, 1);
  target += normal_lpdf(phi_2_p2_p1 | 0, 1);
  target += normal_lpdf(phi_2_p2_p2 | 0, 1);
  target += normal_lpdf(phi_2_p2_p3 | 0, 1);
  target += normal_lpdf(phi_2_p2_p4 | 0, 1);

  target += gamma_lpdf(phi_2 | phi_2_p1 * phi_2_p1_p3, (phi_2_p1 * phi_2_p1_p3) ./ phi_2_p2 * phi_2_p3);

  // likelihood
  if(! prior_only) {
    target += gamma_lpdf(precision / phi_2_p3 | phi_2_p1 * phi_2_p1_p3, (phi_2_p1 * phi_2_p1_p3) ./ phi_2_p2[index_samples_to_measurements][index_measurement_has_precision] * phi_2_p3);
    target += beta_lpdf(m | mu .* phi, (1 - mu) .* phi);
  }

}

generated quantities {

  vector<lower = 0, upper = 1>[N] m_rep;
  vector[N] log_lik;

  for(n in 1:N) {
    m_rep[n] = beta_rng(mu[n] * phi[n], (1 - mu[n]) * phi[n]);
    log_lik[n] = beta_lpdf(m[n] | mu[n] .* phi[n], (1 - mu[n]) .* phi[n]);
  }

}
