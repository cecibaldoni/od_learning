library(dplyr)
library(tidyr)
library(rstan)
library(ggplot2)
library(bayesplot)

# Set random seed for reproducibility
set.seed(123)

# Define number of individuals, behaviors, and observations
N <- 1000 # number of observations
ID <- rep(1:20, each = 10) # 20 individuals, each with 10 observations
trial <- rep(1:10, times = 20) # 10 trials for each individual
season <- sample(c("winter", "spring", "summer"), size = N, replace = TRUE) # 3 seasons
Sn <- as.numeric(factor(season, levels = c("winter", "spring", "summer"))) # convert to integers, 1 to 3
success <- rbinom(N, 1, prob = ifelse(S == "summer", 0.8, 0.2)) # probability of success is higher in summer

# Save data to a CSV file
sim_data <- data.frame(ID, trial, season, Sn, success)
write.csv(sim_data, file = "sim_data.csv", row.names = FALSE)

# Load the simulated data from the CSV file
data <- read.csv("sim_data.csv")

# Assign values to variables
N <- nrow(data)  # Number of observations
K <- 2  # Number of behaviors (success or failure)
J <- length(unique(data$ID))  # Number of individuals
#S <- as.factor(data$season)
tech <- data$success + 1  # Technique chosen (0 or 1 -> 1 or 2)
bout <- data$trial  # Processing bout per individual
id <- data$ID  # Individual ID
N_effects <- 1
# Rename the columns to match the variable names in the model code
names(data) <- c("ID", "trial", "season", "success")

# Define the model in stan, and load the .stan file
model <- stan_model("learning_mod.stan")
print(model)
# Fit the model to the data, while defining the data for the model
fit <- stan(
  model_code = model@model_code,
  data = list(
    N = N,
    K = K,
    J = J,
    tech = tech,
    bout = bout,
    id = id,
    N_effects = N_effects
  ), 
  chains = 4, iter = 2000, warmup = 1000
)

list_draws <- extract(fit)
print(names(list_draws))

# Print a summary of the fitted model
#lp__ (the log of posterior density function up to an additive constant)
print(fit, pars=c("theta"))

#The default plot shows posterior uncertainty intervals (by default 80% (inner) and 95% (outer)) and the posterior median for all the parameters
plot(fit, pars=c("theta", "lp__"))

traceplot(fit, pars=c("theta"), inc_warmup= TRUE, nrow=2)

# Extracting the posterior draws
theta_draws = extract(fit)$theta

# Calculating posterior mean (estimator)
mean(theta_draws)
## [1] 0.2715866
# Calculating posterior intervals
quantile(theta_draws, probs=c(0.10, 0.90))
##       10%       90% 
## 0.1569165 0.3934832
theta_draws_df = data.frame(list(theta=theta_draws))
plotpostre = ggplot(theta_draws_df, aes(x=theta.1)) +
  geom_histogram(bins=20, color="gray")
plotpostre


######################
#Newmodel with season#
######################

# Define the model in stan, and load the .stan file
model2 <- stan_model("learing_season.stan")
# Fit the model to the data, while defining the data for the model
fit1 <- stan(
  model_code = model2@model_code,
  data = list(
    N = N,
    K = K,
    J = J,
    Sn = Sn,
    tech = tech,
    bout = bout,
    id = id,
    N_effects = N_effects
  ), 
  chains = 4, iter = 2000, warmup = 1000
)

list_draws1 <- extract(fit1)
print(names(list_draws1))

# Print a summary of the fitted model
#lp__ (the log of posterior density function up to an additive constant)
print(fit1, pars=c("season_effect"))

#The default plot shows posterior uncertainty intervals (by default 80% (inner) and 95% (outer)) and the posterior median for all the parameters
plot(fit1, pars=c("theta", "season_effect"))

traceplot(fit, pars=c("theta"), inc_warmup= TRUE, nrow=2)

# Extracting the posterior draws
theta_draws = extract(fit)$eta

# Calculating posterior mean (estimator)
mean(theta_draws)
## [1] 0.2715866
# Calculating posterior intervals
quantile(theta_draws, probs=c(0.10, 0.90))
##       10%       90% 
## 0.1569165 0.3934832
theta_draws_df = data.frame(list(theta=theta_draws))
plotpostre = ggplot(theta_draws_df, aes(x=theta)) +
  geom_histogram(bins=20, color="gray")
plotpostre