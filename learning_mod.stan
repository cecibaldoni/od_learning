data {
    int K;              // num behaviors
    int N;              // num observations in dataset
    int J;              // num individuals
    int tech[N];        // techique chosen
    int bout[N];        // processing bout per individual
    int id[N];          // individual id
    int N_effects;      // number of learning parameters to estimate
}

parameters {
  vector[J] theta;      // theta = vector of lenght J that represents the Individual learning rate
  //Each theta[j] represents the learning rate of the j-th individual in the population
} 

model {
  vector[K] AC;  // Attraction scores
  real logPrA;  // log Asocial learning Probability for observation i
   // Priors
  theta ~ normal(0, 1); // Prior on individual learning rates
  // Likelihood loop
  for (i in 1:N) {
    // Update attractions
    for (j in 1:K) {
      AC[j] = 1;
    }
    logPrA = AC[tech[i]] - log_sum_exp(AC);
    target += logPrA; //log-posterior density of the model
  }
} //goal of the model is to estimate the individual learning rates (theta) 
//that best explain the observed choices, given the attraction scores and the asocial learning process

generated quantities {
  vector[K] AC;  // Attraction scores
  real logPrA;  // Individual learning temp
  matrix[N, K] PrPreds;
  for (i in 1:N) {
    // Update attractions
    for (j in 1:K) {
      AC[j] = 1;
    }
    logPrA = AC[tech[i]] - log_sum_exp(AC);
    for (j in 1:K) {
      PrPreds[i, j] = exp(AC[j] - log_sum_exp(AC));
    }
  }
}
