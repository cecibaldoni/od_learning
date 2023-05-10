data {
  int K;              // num behaviors
  int N;              // num observations in dataset
  int J;              // num individuals
  int tech[N];        // techique chosen
  int bout[N];        // processing bout per individual
  int id[N];          // individual id
  int N_effects;      // number of learning parameters to estimate
  int Sn[N];          // season of observation
}

parameters {
  vector[J] theta;     
  vector[3] season_effect; // effect of season on learning rates
}

model {
  vector[K] AC;  // Attraction scores
  real logPrA;  // log Asocial learning Probability for observation i
  // Priors
  theta ~ normal(0, 1); // Prior on individual learning rates
  season_effect ~ normal(0, 1); // Prior on season effects
  // Likelihood loop
  for (i in 1:N) {
    // Update attractions
    for (j in 1:K) {
      AC[j] = 1 + season_effect[Sn[i]];
    }
    logPrA = AC[tech[i]] - log_sum_exp(AC);
    target += logPrA; //log-posterior density of the model
  }
}

generated quantities {
  vector[K] AC;  // Attraction scores
  real logPrA;  // Individual learning temp
  matrix[N, K] PrPreds;
  for (i in 1:N) {
    // Update attractions
    for (j in 1:K) {
      AC[j] = 1 + season_effect[Sn[i]];;
    }
    // Add season effect to individual learning rates
    logPrA = AC[tech[i]] - log_sum_exp(AC);
    for (j in 1:K) {
      PrPreds[i, j] = exp(AC[j] - log_sum_exp(AC));
    }
  }
}
