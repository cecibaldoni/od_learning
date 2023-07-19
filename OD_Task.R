library(dplyr) # for data manipulation and plots
library(haven) #for reading sav data
library(sjstats) #for calculating intra-class correlation (ICC)
library(ROCR)#for calculating area under the curve (AUC) statistics
library(pROC)
library(brms) #for Bayesian (multilevel) generalised linear modelling
library(modelr) #for data manipulation
library(tidybayes) #for analysis of posterior draws of a Bayesian model
library(ggplot2)
# library(rstan)
library(tidyr)

od_data = read.csv(file = "~/data/od_learning/OD_Task.csv",
                   header = TRUE, sep = ",", dec = ".", na.strings = "NA")
#od_data = read.csv(file = "OD_Task.csv",
#                   header = TRUE, sep = ",", dec = ".", na.strings = "NA")

# make winter reference category
od_data$seasonF <- relevel(as.factor(od_data$season), ref = "winter")
od_data$test_ID <- factor(od_data$test_ID, ordered = TRUE)
od_data$ID <- as.factor(od_data$ID)


subset <- od_data %>%
  select(ID, season, seasonF, test_ID, success)
# str(subset)

###############
#what is this again
contrasts(od_data$test_ID, how.many=3) <- contr.poly(10)
contrasts(od_data$test_ID)
#if I want second order I just have linear+quadratic -> how.many=2
#####
#I can try to give it more space by how.many=4

## which priors you can set
get_prior(success ~ seasonF + (1|ID), 
          data = subset, family = bernoulli(link="logit"))

prior1 <- c(prior(normal(-1, 1), class = Intercept),
            prior(normal(1, 0.5), class = b, coef = seasonFsummer),
            prior(normal(0, 0.5), class = b, coef = seasonFspring),
            prior(normal(0, 1), class = sd, group = ID))
#we have used a normal distribution with a mean of -1 and a standard deviation of 1 for the intercept, 
#reflecting our belief that success rates are lower in winter compared to other seasons. 
#We have also used a normal distribution with a mean of 1 and a standard deviation of 0.5 for the coefficient of "season[summer]", 
#reflecting our belief that success rates are higher in summer compared to winter, and a normal distribution with a mean of 0 and a standard deviation of 0.5 
#for the coefficient of "season[spring]", reflecting our belief that success rates are intermediate in spring. 
#Finally, we have used a normal distribution with a mean of 0 and a standard deviation of 1 for the random intercept, 
#reflecting our belief that there is some variation in success rates among individuals that is not explained by the "season" variable.

## get posterior predictions
Bayes_Model_Binary1 <- brm(formula = success ~ seasonF + (1|ID), 
                          data = subset, family = bernoulli(link="logit"), 
                          warmup = 500, iter = 20000, chains = 4, thin = 10, 
                          cores = 4, #backend = "cmdstanr", 
                          prior = prior1)

prior_summary(Bayes_Model_Binary1)
plot(Bayes_Model_Binary1)
mcmc_plot(Bayes_Model_Binary1,
          type = "trace")
mcmc_plot(Bayes_Model_Binary1, type = "hist") #show histograms of the posterior distributions
mcmc_plot(Bayes_Model_Binary1) #plot posterior intervals
mcmc_plot(Bayes_Model_Binary1,
          type = "acf_bar")

summary(Bayes_Model_Binary1)
#it cannot estimate the variance from ID, probably there's no effect of ID on predicting success rate
mcmc_plot(Bayes_Model_Binary1,
          type = "areas",
          prob = 0.95)
#No significative effect

get_variables(Bayes_Model_Binary1)


#######################
###new changes, 4/11###
#######################
get_prior(success ~ seasonF*test_ID + (1+test_ID|ID), 
          #(1+test_ID|ID): individual have different trajectories, but same trajectories between seasons
          data = od_data, family = bernoulli(link="logit"))
#lkj(1)

bm_prior_new <- c(prior(normal(0,1), class = b, coef = seasonFsummer),
                  prior(normal(0,1), class = b, coef = seasonFspring),
                  prior(normal(0,1), class = b, coef = test_ID.C),
                  prior(normal(0,1), class = b, coef = test_ID.L),
                  prior(normal(0,1), class = b, coef = test_ID.Q),
                  prior(normal(0,1), class = b, coef = seasonFspring:test_ID.C),
                  prior(normal(0,1), class = b, coef = seasonFspring:test_ID.Q),
                  prior(normal(0,1), class = b, coef = seasonFspring:test_ID.L),
                  prior(normal(0,1), class = b, coef = seasonFsummer:test_ID.C),
                  prior(normal(0,1), class = b, coef = seasonFsummer:test_ID.L),
                  prior(normal(0,1), class = b, coef = seasonFsummer:test_ID.Q),
                  prior(exponential(1), class = sd))
#i could change the class = sd to the default (by not stating in the prior); the default is assuming that all individuals can behave VERY differently from each other
#exponential 1 says that they cant be so different from each other
#if the model behaves well I can then try to remove the sd and use the default, if the model fails with the default, then i put back the exponential(1)


# prior1 <- get_prior(success ~ seasonF*test_ID + (1+test_ID|ID), 
#                    data = od_data, family = bernoulli(link="logit"))

make_stancode(success ~ seasonF*test_ID + (1+test_ID|ID), 
              data = od_data, family = bernoulli(link="logit"), prior = bm_prior_new)
#(1+seasonF*test_ID|ID) this is not completely legitimate, because the number of trials per shrew are unbalanced between season (some shrew only have 10 (1 season) and very few 30( 3 season))
#we say theat ther is a unique effect of season per ID, but we cannot be sure about it with the data being so unbalanced 
Bayes_Model_Binary_new <- brm(formula = success ~ seasonF*test_ID + (1+seasonF*test_ID|ID), 
                              data = od_data, family = bernoulli(link="logit"), 
                              warmup = 500, iter = 2000, chains = 4, #thin = 10, 
                              cores = 4, 
                              prior = bm_prior_new)
Bayes_Model_Binary_new
#new number with the new priors, but result is consistent, they are more succesfull in summer, however the curve looks like, it is shifted upwards in summer compared to the intercept (winter)
#seasonFsummer               0.74      0.36     0.05     1.47 1.00     7253     4239
#upper and lower levels are on one side of the 0, so positive effect of season to success

#Test_ID.L ID.Q ID.C are of winter (intercept)
#everything else underneath are the adjustments to the winter effects (e.g. seasonFspring:testid 1.70 = -0.93+1.70)

#(1+test_ID|ID) individual have different trajectories, but same trajectories between seasons
#(1+seasonF*test_ID|ID) individuals have different trajectories, and i expect different trajectories per ID per season.

#everything else in the results share the explanatory effort, nothing stands out (except summer!)
#fitted model is numerically stable

#######################################
###New model, coping with unbalanced data
#######################################
bm_prior <- c(prior(normal(0,1), class = b, coef = seasonFsummer),
                  prior(normal(0,1), class = b, coef = seasonFspring),
                  prior(normal(0,1), class = b, coef = test_ID.C),
                  prior(normal(0,1), class = b, coef = test_ID.L),
                  prior(normal(0,1), class = b, coef = test_ID.Q),
                  prior(normal(0,1), class = b, coef = seasonFspring:test_ID.C),
                  prior(normal(0,1), class = b, coef = seasonFspring:test_ID.Q),
                  prior(normal(0,1), class = b, coef = seasonFspring:test_ID.L),
                  prior(normal(0,1), class = b, coef = seasonFsummer:test_ID.C),
                  prior(normal(0,1), class = b, coef = seasonFsummer:test_ID.L),
                  prior(normal(0,1), class = b, coef = seasonFsummer:test_ID.Q),
                  prior(exponential(1), class = sd))
#not as rich of a model, but i cannot learn about the unique season effect per ID yet
Bayes_Model_Binary <- brm(formula = success ~ seasonF*test_ID + (1+test_ID|ID), 
                              data = od_data, family = bernoulli(link="logit"), 
                              warmup = 500, iter = 2000, chains = 4, #thin = 10, 
                              cores = 4, 
                              prior = bm_prior_new)
summary(Bayes_Model_Binary)
#less complex, but still summer is supported.


#make plot with average trend and each shrew trend per season
library(here)
library(magrittr)
library(dplyr)
library(purrr)
library(forcats)
library(tidyr)
library(modelr)
library(ggdist)
library(tidybayes)
library(ggplot2)
library(cowplot)
library(rstan)
library(brms)
library(ggrepel)
library(RColorBrewer)
library(gganimate)
library(posterior)
theme_set(theme_tidybayes() + panel_border())

#se the get_variables() function to get a list of raw model variable names 
#so that we know what variables we can extract from the model
get_variables(Bayes_Model_Binary_new)
get_variables(Bayes_Model_Binary)
#https://cran.r-project.org/web/packages/tidybayes/vignettes/tidy-brms.html
#https://github.com/mvuorre/brmstools

#with this specification (r_condition[condition,term]) spread_draws() splits the variable indices by commas and spaces (you can split by other characters by changing the sep argument). 
#It lets you assign columns to the resulting indices in order. So r_ID[whatever,Intercept] has indices ID and intercept, and spread_draws() lets us extract these indices as columns in the resulting tidy data frame of draws from r_condition
Bayes_Model_Binary_new %>% 
  spread_draws(r_ID[20210802-1,Intercept]) %>% 
  head(10)
#In this particular model, there is only one term (Intercept), 
#thus we could omit that index altogether to just get each condition and the value of r_condition for that condition
Bayes_Model_Binary_new %>% 
  spread_draws(r_ID[20210802-1,]) %>% 
  head(10)

Bayes_Model_Binary_new %>%
  spread_draws(b_Intercept) %>%
  median_qi()


Bayes_Model_Binary_new %>%
  spread_draws(b_Intercept, r_ID[ID,], ) %>%
  median_qi(condition_mean = b_Intercept + r_ID, .width = c(.95, .66)) %>%
  ggplot(aes(y = ID, x = ID_mean, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() 
#If you would rather have a long-format list of intervals, use gather_draws() instead:
Bayes_Model_Binary_new %>%
  gather_draws(b_Intercept) %>%
  median_qi()

Bayes_Model_Binary_new %>%
  spread_draws(r_ID[ID,]) %>%
  median_qi()

Bayes_Model_Binary_new %>%
  spread_draws(b_Intercept, r_ID[ID,]) %>%
  summarise_draws()
#Within each draw, b_Intercept is repeated as necessary to correspond to every index of r_condition. 
#Thus, the mutate function from dplyr can be used to find their sum, condition_mean (which is the mean for each condition)
#we can simplify by moving the calculation of condition_mean from mutate into median_qi()
#median_qi() and its sister functions can produce an arbitrary number of probability intervals by setting the .width = argument
Bayes_Model_Binary_new %>%
  spread_draws(b_Intercept, r_ID[ID,]) %>%
  median_qi(condition_mean = b_Intercept + r_ID, .width = c(.95, .8, .5))
#The results are in a tidy format: one row per group and uncertainty interval width (.width). This facilitates plotting

Bayes_Model_Binary_new %>%
  spread_draws(b_Intercept, r_ID[ID,]) %>%
  median_qi(condition_mean = b_Intercept + r_ID, .width = c(.95, .66)) %>%
  ggplot(aes(y = ID, x = condition_mean, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() 

Bayes_Model_Binary1 %>%
  spread_draws(b_Intercept, r_ID[ID,]) %>%
  median_qi(condition_mean = b_Intercept + r_ID, .width = c(.95, .66)) %>%
  ggplot(aes(y = ID, x = condition_mean, xmin = .lower, xmax = .upper)) +
  stat_eye() +
  geom_vline(xintercept = c(-.5, -.1), linetype = "dashed") +
  scale_fill_manual(values = c("gray80", "skyblue"))
  #instead of stat_halfeye() I can try
  #geom_pointinterval() 




spaghetti_plot <- plot(
  conditional_effects(
    Bayes_Model_Binary_new,
    effects = "test_ID",
    prob = 0.89,
    spaghetti = TRUE
  ),
  theme = theme_bw()
)

#sample for sleepstudy is 18 subjects, 10 trial per subject (fully balanced)

#WTF. Before it worked and now it doesn't?
conditional_effects(Bayes_Model_Binary_new, conditions = distinct(success, ID))

out_f <- conditional_effects(Bayes_Model_Binary, spaghetti = TRUE)
                             #effects= "seasonF", conditions = distinct(od_data, test_ID))[[1]]
str(out_f)
#combination of fixed effects (season with trial_ID)

out_r <- conditional_effects(Bayes_Model_Binary, 
                             #effects = "seasonF", 
                             conditions = distinct(od_data, ID))
                             #re_formula = NULL)[[1]]
out_r
out_f
plot(out_f)
table(od_data$ID)

out_r %>% 
  ggplot(aes(seasonF*test_ID, estimate__)) +
  geom_ribbon(
    data = out_f,
    aes(ymin = lower__, ymax = upper__),
    alpha = .33
  ) +
  geom_line(data = out_f, size=2) +
  geom_line(aes(group = ID))



out_f <- marginal_effects(Bayes_Model_Binary_new)[[1]]
out_r <- marginal_effects(Bayes_Model_Binary_new, conditions = distinct(success, ID), re_formula = NULL)[[1]]
out_r %>% 
  ggplot(aes(test_ID, estimate__)) +
  geom_ribbon(
    data = out_f,
    aes(ymin = lower__, ymax = upper__),
    alpha = .33
  ) +
  geom_line(data= out_f, size=2) +
  geom_line(group = ID)
######################################
##
######################################
fitted_values <- fitted(Bayes_Model_Binary)
fitted_values

dat <- as.data.frame(cbind(Y = standata(Bayes_Model_Binary)$Y, fitted_values))
ggplot(dat) + geom_point(aes(x = Estimate, y = Y))

# hypothesis testing 
h1 <- hypothesis(Bayes_Model_Binary_new, c("seasonFsummer > 0", "seasonFsummer <0" ))
print(h1,digits = 3)
cond_eff <- conditional_effects(Bayes_Model_Binary_new)

  
df <- as.data.frame(cond_eff$seasonF)
cond_plot <- ggplot(df, aes(x=seasonF, y = estimate__)) +
  scale_color_manual(values = c("blue", "green", "red")) +
  geom_point(aes(color=seasonF), size = 4) +
  geom_linerange(aes(ymin = lower__, ymax = upper__), color = "darkgrey") +
  labs(y = "Success Rate", x = "Season") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__)) +
  theme(legend.position="none")

cond_plot
df$upper__
df$lower__

Season <- c("2", "3", "1")
new<- df %>%
  mutate(seasonF, levels = Season) %>%
  arrange(Season)
new$seasonF <- relevel(new$seasonF, ref = "summer")
cond_2 <- ggplot(new, aes(x=seasonF, y = estimate__)) +
  scale_color_manual(values = c("red", "blue", "green")) +
  geom_point(aes(color=seasonF), size = 4) +
  geom_linerange(aes(ymin = lower__, ymax = upper__), color = "darkgrey") +
  labs(y = "Success Rate", x = "Season") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__)) +
  theme(legend.position="none")
cond_2

Prob <- predict(Bayes_Model_Binary, type="response")
Prob <- Prob[,1]
Pred <- prediction(Prob, as.vector(pull(subset, success)))
AUC <- performance(Pred, measure = "auc")
AUC <- AUC@y.values[[1]]
AUC
#With an AUC score of close to 0.60, the model does not discriminate well.


subset$test_ID <- as.ordered(subset$test_ID)
subset$success <- as.numeric(subset$success)
subset$season <- as.factor(subset$season)


plot2 <- ggplot(od_data, aes(x = as.numeric(test_ID), y = success, color = season)) +
  geom_jitter(width =  0.05, height = 0.05) + 
  geom_smooth(mapping = aes(x = as.numeric(test_ID), y = success), ) + 
  labs(x = "Trial", y = "Success rate") +
  labs(col="Season") +
  scale_color_manual(values = c("green", "#D16103","#4E84C4")) 
plot2


df <- subset %>%  group_by(test_ID,season) %>%  summarise_at(vars(success), list(mean_success= mean,sd_success=sd))
df$se_success <-df$sd_success/sqrt(3)  #CONTROLLA SE HO CALCOLATO L'ERRORE STANDARD CORRETTAMENTE
plot3 <- ggplot(df, aes(x = test_ID, y = mean_success, color = season)) +
  geom_point(position = position_dodge(width = .6))+geom_errorbar(aes(ymin=mean_success-se_success, ymax=mean_success+se_success), width=.2,position = position_dodge(width = .6))
plot3


subset1 <- od_data %>%
  summarise(success = sum(success),
            total = n())
subset1

Bayes_Model_Prop <- brm(success | trials(250) ~ season + test_ID,
                        data = subset,
                        family = binomial(link = "logit"),
                        warmup = 500,
                        iter = 2000,
                        chains = 4,
                        #thin = 10,
                        cores = 4)

mcmc_plot(Bayes_Model_Prop,
          type = "trace")
mcmc_plot(Bayes_Model_Prop,
          type = "acf_bar")
summary(Bayes_Model_Prop)
fixef(Bayes_Model_Prop)
model_easy <- brm(success ~ season + test_ID,
                        data = subset,
                        family = binomial(link = "logit"),
                        warmup = 500,
                        iter = 2000,
                        chains = 4,
                        #thin = 10,
                        cores = 4)
summary(model_easy)
get_variables(model_easy)
subset %>% 
  data_grid(test_ID) %>% 
  add_epred_draws(model_easy, dpar = TRUE, category = "season") %>% 
  ggplot(aes(x = test_ID, y = .epred, color = season)) +
  stat_pointinterval(position = position_dodge(width = .4)) +
  scale_size_continuous(guide = "none") +
  scale_color_manual(values = brewer.pal(6, "Blues")[-c(1,2)])

subset %>% 
  add_epred_draws(model_easy, category = season) %>% 
  ggplot(aes(x = test_ID, y = .epred, color = season)) +
  stat_pointinterval(position = position_dodge(width = .4)) +
  scale_size_continuous(guide = "none")


subset %>% 
  select(test_ID) %>% 
  add_predicted_draws(model_easy, seed = 1234) %>% 
  ggplot(aes(x = test_ID, y = success)) +
  geom_count(color = "gray75") +
  geom_point(aes(fill = season), data = subset, shape = 21, size = 2) +
  scale_fill_brewer(palette = "Dark2")


#are there any other parameters that can explain the data?
library(lme4)
library(lmerTest)

subset2 <- select(od_data, ID, seasonF, test_ID, success, latency, test_time, c.of.d.)
subset2$test_ID <- factor(subset2$test_ID, ordered = TRUE)

str(subset2)
model1 <- lmer(formula = success ~ c.of.d. + scale(latency) + 
                 scale(test_time) + seasonF + (1|ID), data = subset2)
summary(model1)
#success is higher in summer? Higher than the Intercept (spring, latency zero, test time zero)

model2 <- lmer(formula = success ~ c.of.d.*seasonF + (1|ID), data = subset2)
summary(model2)

model3 <- lmer(formula = success ~ c.of.d.*test_ID  + (1|ID), data = subset2)
summary(model3)
#i want to see if change of directions (c.of.d.) increase/decrease 
#with the increasing/decreasing of the test_ID number, 
model4 <- lmer(formula = c.of.d. ~ test_ID  + (1|ID), data = subset2)
summary(model4)


### add ordered test number
subset$test_ID <- factor(subset$test_ID, ordered = TRUE)

Bayes_Model_Binary_2 <- brm(formula = success ~  (1|ID) + test_ID*seasonF, 
                          data = subset, family = bernoulli(link="logit"), 
                          warmup = 500, iter = 20000, chains = 4, thin = 10, cores = 4,
                          prior = bm_prior2)

mcmc_plot(Bayes_Model_Binary_2,
          type = "trace")
mcmc_plot(Bayes_Model_Binary_2, type = "hist") #show histograms of the posterior distributions
mcmc_plot(Bayes_Model_Binary_2) #plot posterior intervals
mcmc_plot(Bayes_Model_Binary_2,
          type = "acf_bar")
plot(Bayes_Model_Binary_2)
summary(Bayes_Model_Binary_2)

mcmc_plot(Bayes_Model_Binary_2,
          type = "areas",
          prob = 0.95)

# hypothesis testing 
h2 <- hypothesis(Bayes_Model_Binary_2, c("seasonFsummer > 0", "seasonFsummer <0" ))
print(h2,digits = 3)

plot(conditional_effects(Bayes_Model_Binary_2))

# complete model
Bayes_Model_Binary_3 <- brm(formula = success ~  (1|ID) + test_ID*seasonF 
                            + latency + test_time + c.of.d., data = od_data, 
                            family = bernoulli(link="logit"), warmup = 500, 
                            iter = 20000, chains = 4, thin = 10, cores = 4, 
                            prior = bm_prior2)

summary(Bayes_Model_Binary_3)

od_data$latency <- as.numeric(od_data$latency)
mcmc_plot(Bayes_Model_Binary_3,
          type = "areas",
          prob = 0.95)
mcmc_plot(Bayes_Model_Binary_3, variable = c("b_c.of.d."),
          type = "areas",
          prob = 0.95)
plot(conditional_effects(Bayes_Model_Binary_3))

?make_conditions
condition <- make_conditions(Bayes_Model_Binary_3, "c.of.d.")

cond <- conditional_effects(Bayes_Model_Binary_3, "c.of.d.", conditions = condition)
plot(cond)

cond2 <- conditional_effects(Bayes_Model_Binary_3, "latency", conditions = condition)
plot(cond2)

h3 <- hypothesis(Bayes_Model_Binary_3, c("c.of.d.>0", "c.of.d.<0" ))
print(h3,digits = 3)
h4 <- hypothesis(Bayes_Model_Binary_3, c("test_time>0", "test_time<0"))
print(h4, digits = 4)
Bayes_Model_Binary_3$data

full_dir_bm <- conditional_effects(Bayes_Model_Binary_3, effects = "c.of.d.", prob = 0.95)
full_season_bm <- conditional_effects(Bayes_Model_Binary_3, effects = "seasonF", prob = 0.95)
season_bm <- conditional_effects(Bayes_Model_Binary_2, effects = "seasonF", prob = 0.95)
print(season_bm)

