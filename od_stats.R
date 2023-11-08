library(tidyverse)
library(brms)
library(cmdstanr)
library(gridExtra)
library(sjstats) #for calculating intra-class correlation (ICC)
library(ROCR) #for calculating area under the curve (AUC) statistics
library(pROC)
library(modelr) #for data manipulation
library(tidybayes) #for analysis of posterior draws of a Bayesian model
library(mgcv)
library(BayesFactor)


od_data = read.csv(file = "~/data/od_learning/OD.csv",
                   header = TRUE, sep = ",", dec = ".", na.strings = "NA")
# od_data = read.csv(file = "~/data/od_learning/Od_sumwin_alltrial.csv",
#                    header = TRUE, sep = ",", dec = ".", na.strings = "NA")

od_data <- od_data %>% 
  rename(trial = test_ID) %>% 
  mutate(season = as.factor(season),
         ID = as.factor(ID),
         status = ifelse(status== "semicaptive", "captive", status),
         status = as.factor(status),
         sex = as.factor(sex)) %>% 
  filter(season != "spring") %>% 
  droplevels()
#new dataset with all the habituations, the following code remove all habituation trials:
od_data <- od_data %>%
  filter(!grepl("^H_", trial)) #The ^ character in the pattern specifies the start of the string

# od_data <- od_data %>% 
#   filter(!grepl("^spring", season))

#od_data$season <- as.integer(ifelse(od_data$season == "summer", 1, 2))
od_data$trial <- factor(as.numeric(od_data$trial), ordered = TRUE)
od_data$trial_n <- as.numeric(od_data$trial)

subset <- od_data[, c("ID", "season", "success", "trial", "trial_n", "status")]

ggplot(subset(od_data, status != "captive"), aes(x = trial_n, y = success, color = season, linetype = status)) +
  geom_jitter(width = 0.05, height = 0.05) + 
  geom_smooth(mapping = aes(x = trial_n, y = success), se = TRUE, alpha= 0.12) +
  labs(x = "Trial", y = "Success rate", col="Season") +
  scale_color_manual(values = c("#D16103","#4E84C4")) +
  scale_x_discrete(breaks = subset(od_data, status != "captive")$trial_n, labels = subset(od_data, status != "captive")$trial_n)

ggplot(od_data, aes(x = trial_n, y = success, color = season, linetype = status)) +
  geom_jitter(width =  0.05, height = 0.05) + 
  geom_smooth(mapping = aes(x = trial_n, y = success), se = TRUE, alpha= 0.12) +
  labs(x = "Trial", y = "Success rate") +
  labs(col="Season") +
  scale_color_manual(values = c("#D16103","#4E84C4")) +
  scale_x_discrete(breaks = od_data$trial_n, labels = od_data$trial_n)

####WILD DATA####

wild_data <- od_data %>%
  filter(status == "wild") %>%
  droplevels()
length(unique(wild_data$ID[wild_data$season == "summer"]))
#create base plot for later use
p_wild <- ggplot(wild_data, aes(x = trial, y = success, color = season)) +
  geom_jitter(width =  0.05, height = 0.05) + 
  geom_smooth(mapping = aes(x = trial_n, y = success), se = FALSE) +
  labs(x = "Trial", y = "Success rate") +
  labs(col="Season") +
  scale_color_manual(values = c("#D16103","#4E84C4")) +
  scale_x_discrete(breaks = wild_data$trial, labels = wild_data$trial)

fit1 <- brm(data = wild_data,
             family = bernoulli,
             success ~ 0 + s(trial_n, by = season, bs = "bs") + season + (0 + trial_n | ID),
             prior = c(prior(normal(0, 1), class = b),
                       prior(exponential(.65), class = sd)),
             chains = 4, iter = 2000,
             save_pars = save_pars(all = TRUE),
             backend = "cmdstanr",
             threads = threading(2),
             control=list(adapt_delta=0.99, max_treedepth = 10))
p1 = pp_check(fit1, type = 'stat', stat = 'min') +
  ggtitle('Prior predictive distribution of minimum values') + 
  #geom_vline(xintercept = minvalue) +
  theme_classic()
p2 = pp_check(fit1, type = 'stat', stat = 'max') +
  ggtitle('Prior predictive distribution of maximum values')  +
  #geom_vline(xintercept = maxvalue) +
  theme_classic()
p3 = pp_check(fit1, type = 'stat', stat = 'mean') +
  ggtitle('Prior predictive distribution of average values')   +
  #geom_vline(xintercept = meanvalue) +
  theme_classic()
p4 = pp_check(fit1, ndraws = 100)
ggpubr::ggarrange(p1,p2,p3, p4)

fit2 <- brm(data = wild_data,
            family = bernoulli,
            success ~ 0 + trial + season + (0  + trial | ID),
            prior = c(prior(normal(0, 1), class = b),
                      prior(exponential(.65), class = sd),
                      prior(lkj(1), class = cor)),
            chains = 4, iter = 2000,
            save_pars = save_pars(all = TRUE),
            backend = "cmdstanr",
            threads = threading(2),
            control=list(adapt_delta=0.8, max_treedepth = 10))
p1 = pp_check(fit2, type = 'stat', stat = 'min') +
  ggtitle('Prior predictive distribution of minimum values') + 
  #geom_vline(xintercept = minvalue) +
  theme_classic()
p2 = pp_check(fit2, type = 'stat', stat = 'max') +
  ggtitle('Prior predictive distribution of maximum values')  +
  #geom_vline(xintercept = maxvalue) +
  theme_classic()
p3 = pp_check(fit2, type = 'stat', stat = 'mean') +
  ggtitle('Prior predictive distribution of average values')   +
  #geom_vline(xintercept = meanvalue) +
  theme_classic()
p4 = pp_check(fit2, ndraws = 100)
ggpubr::ggarrange(p1,p2,p3, p4)

get_prior(success ~ 0 + s(trial_n, by = season, bs = "bs") + season, data = wild_data)
fit3 <- brm(data = wild_data,
             family = bernoulli,
             success ~ 0 + s(trial_n, by = season, bs = "bs") + season,
             prior = c(prior(normal(0, 1), class = "b", coef = "seasonsummer"),
                       prior(normal(0,1), class= "b", coef = "strial_n:seasonsummer_1"),
                       prior(exponential(.65), class = "sds")),
             chains = 4, iter = 2000,
             save_pars = save_pars(all = TRUE),
             backend = "cmdstanr",
             threads = threading(2),
             control=list(adapt_delta=0.99, max_treedepth = 10))

p1 = pp_check(fit3, type = 'stat', stat = 'min') +
  ggtitle('Prior predictive distribution of minimum values') + 
  #geom_vline(xintercept = minvalue) +
  theme_classic()
p2 = pp_check(fit3, type = 'stat', stat = 'max') +
  ggtitle('Prior predictive distribution of maximum values')  +
  #geom_vline(xintercept = maxvalue) +
  theme_classic()
p3 = pp_check(fit3, type = 'stat', stat = 'mean') +
  ggtitle('Prior predictive distribution of average values')   +
  #geom_vline(xintercept = meanvalue) +
  theme_classic()
p4 = pp_check(fit3, ndraws = 100)
ggpubr::ggarrange(p1,p2,p3, p4)

fit1 = add_criterion(fit1, c("loo", "loo_R2"), 
                      control=list(adapt_delta=0.8, max_treedepth = 10), 
                      backend = "cmdstanr",
                      reloo = TRUE)
fit2 = add_criterion(fit2, c("loo", "loo_R2"), 
                      control=list(adapt_delta=0.8, max_treedepth = 10), 
                      backend = "cmdstanr",
                      reloo = TRUE)
fit3 = add_criterion(fit3, c("loo", "loo_R2"), 
                      control=list(adapt_delta=0.8, max_treedepth = 10), 
                      backend = "cmdstanr",
                      reloo = TRUE)
loo_results <- loo_compare(fit1, fit2, fit3)
loo_results


results=summary(fit3)
results=data.frame(results$fixed)
results$covariate=rownames(results)
posteriors <- ggplot(results,aes(x=Estimate,y=covariate,color=I("blue")))+geom_point()+geom_linerange(aes(xmin=l.95..CI,xmax=u.95..CI,color=I("blue")))+
  geom_vline(xintercept = 0)+theme_classic()
posteriors

wd <-  wild_data %>% 
  distinct(ID, status, trial, trial_n, season, success)

fitted <- fitted(fit3, newdata = wd) %>% 
  data.frame() %>% 
  bind_cols(wd) %>% 
  mutate(#ID     = str_c("ID[", ID, "]"),
    season = factor(season),
    trial = factor(trial))

fitted_avg <- fitted %>%
  group_by(trial, season) %>%
  summarize(AvgEstimate = mean(Estimate),
            Lower = mean(Q2.5),
            Upper = mean(Q97.5))


fitted_avg_prob <- fitted_avg %>%
  mutate(AvgProbability = logistic_function(AvgEstimate),
         LowerProbability = logistic_function(Lower),
         UpperProbability = logistic_function(Upper))

plot_with_fitted <- ggplot(wild_data, aes(x = trial, y = success, color = season)) +
  geom_jitter(width =  0.05, height = 0.05) + 
  #geom_smooth(mapping = aes(x = trial_n, y = success), se = FALSE) +
  labs(x = "Trial", y = "Success rate") +

  labs(color="Season Model", fill="CI 95%", linetype="Season Smooth") +
  
  scale_color_manual(values = c("#D16103","#0A33A9"), labels = c("Summer", "Winter")) +
  scale_fill_manual(values = c("#D16103","#0A33A9"), labels = c("Summer", "Winter")) +
  scale_linetype_manual(values = c("solid", "dashed"), labels = c("Summer Smooth", "Winter Smooth")) +
  
  scale_x_discrete(breaks = wild_data$trial, labels = wild_data$trial) +
  geom_line(data = fitted_avg, aes(x = trial, y = AvgEstimate, group = season)) +
  geom_ribbon(data = fitted_avg, aes(x = trial, ymin = Lower, ymax = Upper, fill = season, group = season), alpha= 0.2, inherit.aes = FALSE)

print(plot_with_fitted)

h2e <- hypothesis(fit3, c("b_seasonsummer>0","b_seasonsummer<0",
                           "b_seasonwinter>0","b_seasonwinter<0"),class="" )
print(h2e, digits = 6)
#You have strong evidence that the effect of summer on success is positive.
#The effect of winter is less certain; there's weak to moderate evidence suggesting it's positive, but the evidence isn't conclusive.

library(emmeans)

EMM=emmeans(fit3, specs = pairwise ~ trial_n)
EMM
plot(EMM)

EMM2=emmeans(fit3, specs = pairwise ~ season)


### Quantify difference inestimates between trial 1 and 6
# Extracting values for key trials and the maximum trial
key_trials <- c(1, 6, 8, max(wild_data$trial_n))

key_fitted_avg <- fitted_avg %>% 
  filter(trial %in% key_trials) %>% 
  mutate(Odds = exp(AvgEstimate),
         LowerOdds = exp(Lower),
         UpperOdds = exp(Upper))

predicted_values <- posterior_predict(fit3, newdata = new_data)

# Generate predicted logits for the key trials for each season
new_data <- expand.grid(trial_n = key_trials,
                        season = unique(wd$season))

predicted_logits <- predict(fit3, newdata = new_data, summary = FALSE)

# Compute the odds for each key trial and season
predicted_odds <- exp(predicted_logits)
odds_season1 <- predicted_odds[,1]
odds_season2 <- predicted_odds[,2]
# Compute the odds ratio for the changes across the intervals for each season
odds_ratio_1_to_6_season1 <- odds_season1[2] / odds_season1[1]
odds_ratio_6_to_7_season1 <- odds_season1[3] / odds_season1[2]
odds_ratio_7_to_max_season1 <- odds_season1[4] / odds_season1[3]

odds_ratio_1_to_6_season2 <- odds_season2[2] / odds_season2[1]
odds_ratio_6_to_7_season2 <- odds_season2[3] / odds_season2[2]
odds_ratio_7_to_max_season2 <- odds_season2[4] / odds_season2[3]

# Store the odds ratios in a data frame for easier visualization and reporting
odds_ratios_df <- data.frame(
  Interval = rep(c("1 to 6", "6 to 7", "7 to max"), 2),
  Season = rep(c("Season1", "Season2"), each=3),
  OddsRatio = c(odds_ratio_1_to_6_season1, odds_ratio_6_to_7_season1, odds_ratio_7_to_max_season1,
                odds_ratio_1_to_6_season2, odds_ratio_6_to_7_season2, odds_ratio_7_to_max_season2)
)

print(odds_ratios_df)

####WINTER DATA####

winter_data <- subset %>%
  filter(season == "winter") %>%
  droplevels()
length(unique(winter_data$ID[winter_data$status == "wild"]))

#create base plot for later use
p_winter <- ggplot(winter_data, aes(x = trial, y = success, color = status)) +
  geom_jitter(width =  0.05, height = 0.05) + 
  geom_smooth(mapping = aes(x = trial_n, y = success), se = FALSE) +
  labs(x = "Trial", y = "Success rate") +
  labs(col="Status") +
  scale_color_manual(values = c("#0AA9A9","#0A33A9")) +
  scale_x_discrete(breaks = winter_data$trial, labels = winter_data$trial)
p_winter

fit1_winter <- brm(data = winter_data,
                    family = bernoulli,
                    success ~ 0 + s(trial_n) + status + (0 + trial_n | ID),
                    prior = c(prior(normal(0, 1.5), class = b),
                              prior(exponential(.65), class = sd)),
                    chains = 4, iter = 2000,
                    save_pars = save_pars(all = TRUE),
                    backend = "cmdstanr",
                    threads = threading(2),
                    control=list(adapt_delta=0.999, max_treedepth = 10))

get_prior(success ~ 0 + s(trial_n, by = status, bs = "bs") + status, data = winter_data)
fit2_winter <- brm(data = winter_data,
                    family = bernoulli,
                    success ~ 0 + s(trial_n, by = status, bs = "bs") + status,
                    prior = c(prior(normal(0, 1), class = "b", coef = "statuswild"),
                              prior(normal(0,1), class= "b", coef = "strial_n:statuswild_1"),
                              prior(exponential(.65), class = "sds")),
                    chains = 4, iter = 2000,
                    #sample_prior = "only",
                    save_pars = save_pars(all = TRUE),
                    backend = "cmdstanr",
                    threads = threading(2),
                    control=list(adapt_delta=0.99, max_treedepth = 15))

fit1_winter <- add_criterion(fit1_winter, c("loo", "loo_R2"), 
                              control=list(adapt_delta=0.8, max_treedepth = 11), 
                              backend = "cmdstanr",
                             reloo = TRUE)
                              
fit2_winter <- add_criterion(fit2_winter, c("loo", "loo_R2"),
                              control=list(adapt_delta=0.8, max_treedepth = 15), 
                              backend = "cmdstanr",
                              reloo = TRUE)

loo_compare(fit1_winter, fit2_winter)
                              
p1 = pp_check(fit2_winter, type = 'stat', stat = 'min') +
  ggtitle('Prior predictive distribution of minimum values') + 
  #geom_vline(xintercept = minvalue) +
  theme_classic()
p2 = pp_check(fit2_winter, type = 'stat', stat = 'max') +
  ggtitle('Prior predictive distribution of maximum values')  +
  #geom_vline(xintercept = maxvalue) +
  theme_classic()
p3 = pp_check(fit2_winter,type = 'stat', stat = 'mean') +
  ggtitle('Prior predictive distribution of average values') +
  #geom_vline(xintercept = meanvalue) +
  theme_classic()
p4 = pp_check(fit2_winter, ndraws = 100)
ggpubr::ggarrange(p1,p2,p3, p4)
                              
results=summary(fit2_winter)
results=data.frame(results$fixed)
results$covariate=rownames(results)
posteriors <- ggplot(results,aes(x=Estimate,y=covariate,color=I("blue"))) +
  geom_point() +
  geom_linerange(aes(xmin=l.95..CI,xmax=u.95..CI,color=I("blue"))) +
  geom_vline(xintercept = 0) +
  theme_classic()
posteriors

#Check 90% Confidence interval.
CI_90 <- posterior_interval(fit2_winter, prob = 0.90)
print(CI_90)


conditional_effects(fit1_winter)
conditional_effects(fit2_winter)

win_d <-  winter_data %>%
  distinct(ID, status, trial, trial_n, season)

# get the fitted draws
fitted_winter <- fitted(fit2_winter,
                 newdata = win_d) %>% 
  data.frame() %>% 
  bind_cols(win_d) %>% 
  mutate(ID = str_c("ID[", ID, "]"),
         status = factor(status),
         trial = factor(trial))

fitted_avg_winter <- fitted_winter %>%
  group_by(trial, status) %>%
  summarize(AvgEstimate = mean(Estimate),
            Lower = mean(Q2.5),
            Upper = mean(Q97.5))
logistic_function <- function(log_odds) {
  exp(log_odds) / (1 + exp(log_odds))
}
fitted_avg_winter_prob <- fitted_avg_winter %>%
  mutate(AvgProbability = logistic_function(AvgEstimate),
         LowerProbability = logistic_function(Lower),
         UpperProbability = logistic_function(Upper))
p_winter +
  scale_x_discrete(breaks = winter_data$trial, labels = winter_data$trial) +
  geom_line(data = fitted_avg_winter, aes(x = trial, y = AvgEstimate, color = status, group = status), linetype = "dashed") +
  geom_ribbon(data = fitted_avg_winter, aes(x = trial, ymin = Lower, ymax = Upper, fill = status, group = status), alpha= 0.2, inherit.aes = FALSE) 
                              
finalplot2 <- ggplot(winter_data, aes(x = trial, y = success, color = status)) +
  geom_jitter(width =  0.05, height = 0.05) +
  
  labs(x = "Trial", y = "Success rate") +
  labs(color="Status", fill="CI 95%") +
  
  scale_color_manual(values = c("#0AA9A9", "#0A33A9")) +
  scale_fill_manual(values = c("#0AA9A9", "#0A33A9")) + # Specify custom colors for geom_ribbon shading
  
  scale_x_discrete(breaks = winter_data$trial, labels = winter_data$trial) +
  geom_line(data = fitted_avg_winter, aes(x = trial, y = AvgEstimate, group = status)) +
  geom_ribbon(data = fitted_avg_winter, 
              aes(x = trial, ymin = Lower, ymax = Upper, fill = status, group = status), 
              alpha = 0.2, inherit.aes = FALSE)

print(finalplot2)



gridExtra::grid.arrange(final_plot, finalplot2)

h2e <- hypothesis(fit2e_winter, c("b_statuscaptive>0","b_statuscaptive<0",
                                  "b_statuswild>0","b_statuswild<0"),class="" )
print(h2e, digits = 6)


key_trials_w <- c(1, 6, 8, max(winter_data$trial_n))

key_fitted_avg_w <- fitted_avg_winter %>% 
  filter(trial %in% key_trials_w) %>% 
  mutate(Odds = exp(AvgEstimate),
         LowerOdds = exp(Lower),
         UpperOdds = exp(Upper))

##### FULL PLOT #####
colors <- c("Summer Wild" = "#D16103", "Winter Wild" = "#093E9A", "Winter Captive" = "#16738D")
fills <- c("Summer Wild" = "#D16103", "Winter Wild" = "#093E9A", "Winter Captive" = "#16738D")
colors <- c("Summer Wild" = "#F25202", "Winter Wild" = "#4482CE", "Winter Captive" = "#0AA9A9")
fills <- c("Summer Wild" = "#F25202", "Winter Wild" = "#4482CE", "Winter Captive" = "#0AA9A9")
labels <- c("Summer Wild", "Winter Wild", "Winter Captive")

captive_data <- fitted_avg_winter %>% filter(status == "captive")
plot_with_captive <- plot_with_fitted +
  geom_line(data = captive_data, aes(x = trial, y = AvgEstimate, group = status, color = status)) +
  geom_ribbon(data = captive_data, aes(x = trial, ymin = Lower, ymax = Upper, fill = status, group = status), alpha = 0.2, inherit.aes = FALSE)

wild_data$group_label <- ifelse(wild_data$season == "summer", "Summer Wild", "Winter Wild")
winter_data$group_label<- ifelse(winter_data$status == "wild", "Winter Wild", "Winter Captive")
fitted_avg$group_label <- ifelse(fitted_avg$season == "summer", "Summer Wild", "Winter Wild")
captive_data$group_label <- "Winter Captive"
captive_data <- captive_data %>% 
  ungroup()

ggplot() +
  geom_jitter(data = wild_data, aes(x = trial, y = success, color = group_label), 
              width = 0.05, height = 0.05) + 
  geom_jitter(data = subset(winter_data, group_label == "Winter Captive"), aes(x = trial, y = success, color = group_label), 
              width = 0.05, height = 0.05) +
  geom_line(data = fitted_avg, aes(x = trial, y = AvgEstimate, color = group_label, group = group_label), size = 1.5) +
  geom_ribbon(data = fitted_avg, aes(x = trial, ymin = Lower, ymax = Upper, fill = group_label, group = group_label), 
              alpha= 0.2) +
  geom_line(data = captive_data, aes(x = trial, y = AvgEstimate, color = group_label, group = group_label), size = 1.5) +
  geom_ribbon(data = captive_data, aes(x = trial, ymin = Lower, ymax = Upper, fill = group_label, group = group_label), 
              alpha = 0.2) +
  labs(x = "Trial", y = "Success rate") +
  labs(color = "Model", fill = "CI 95%", linetype = "Smooth") +
  scale_color_manual(values = c("Summer Wild" = "#D16103",
                                "Winter Captive" = "#0AA9A9",
                                "Winter Wild" = "#0A33A9"),
                     labels = c("Summer Wild", "Winter Captive", "Winter Wild")) +
  scale_fill_manual(values = c("Summer Wild" = "#D16103",
                               "Winter Captive" = "#0AA9A9",
                               "Winter Wild" = "#0A33A9"),
                    labels = c("Summer Wild", "Winter Captive", "Winter Wild")) +
  scale_x_discrete(breaks = wild_data$trial, labels = wild_data$trial)



###REINFORCEMENT LEARNING (SIDE)####
# unsuccess + left = 1 they chose left and were unsuccesful
# unsuccess + right = 2 they chose right and were unsuccessful
# success + right = 3 they chose right and were successful
# success + left = 4 they chose left and were successful

#if they were successful, will they choose the same side in the next trial?
od_data$combi <- 1 + od_data$left + 2*od_data$success

dat <- od_data %>% 
  select(ID, season, status, trial, success, left, combi) %>%
  mutate(
    prev_combi = lag(combi),
    #next_combi = lead(combi),
    consecutive_match = ifelse(
      trial == 1,
      NA,
      ((combi %in% c(1, 4)) & (prev_combi %in% c(2, 4))) | (
        (combi %in% c(2, 3)) & (prev_combi %in% c(1, 3))))) #%>%
#If they would decide based on which side is rewarded, it would be all true

dat <- dat %>%
  mutate(prev_combi = lag(combi),
    consecutive_match = ifelse(
      trial == 1,
      NA,
      ((combi %in% c(1, 4)) & (prev_combi %in% c(2, 4))) | (
        (combi %in% c(2, 3)) & (prev_combi %in% c(1, 3))))) %>%
  mutate(side = case_when(
    combi %in% c(1, 4) ~ "left",
    combi %in% c(2, 3) ~ "right"))
#select(trial, combi, consecutive_match)
#If they would decide based on which side is rewarded, it would be all true

ggplot(dat, aes(x = season, fill = consecutive_match)) +
  geom_bar(position = "dodge") +
  facet_wrap(~ status, ncol = 1) +
  labs(title = "Consecutive Match", x = "Season", y = "Count") +
  scale_fill_manual(values = c("TRUE" = "lightblue", "FALSE" = "lightgreen"))
ggplot(dat, aes(x = ID, fill = consecutive_match)) +
  geom_bar(position = "dodge") +
  labs(title = "Consecutive Match by ID", x = "ID", y = "Count") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = c("TRUE" = "lightblue", "FALSE" = "lightgreen"))

#check if there's significance in season and status
table_season <- table(dat$consecutive_match, dat$season)
contingencyTableBF(table_season, sampleType = "indepMulti", fixedMargin = "rows")
#BF= 0.29, evidence in favor of independence or no association.
#the seasons do not appear to have a substantial impact on whether consecutive matches occurred

table_status <- table(dat$consecutive_match, dat$status)
contingencyTableBF(table_status, sampleType = "indepMulti", fixedMargin = "rows")
#status does not have a big impact either

#let's check only the wild guys:
wild <- od_data %>% 
  filter(status == "wild") %>% 
  mutate(combi = 1 + left + 2 * success) %>% 
  select(ID, season, status, trial, success, left, combi)

wild <- wild %>%
  mutate(prev_combi = lag(combi),
         #next_combi = lead(combi),
         consecutive_match = ifelse(trial == 1, NA,
                                    ((combi %in% c(1, 4)) & (prev_combi %in% c(2, 4))) | ((combi %in% c(2, 3)) & (prev_combi %in% c(1, 3))))) %>% 
  mutate(side = case_when(
    combi %in% c(1, 4) ~ "left",
    combi %in% c(2, 3) ~ "right"))
ggplot(wild, aes(x = season, fill = consecutive_match)) +
  geom_bar(position = "dodge") +
  #facet_wrap(~ status, ncol = 1) +
  labs(title = "Consecutive Match by Season", x = "Season", y = "Count") +
  scale_fill_manual(values = c("TRUE" = "lightblue", "FALSE" = "lightgreen"))

#check if there's significance in season using Bayesian Chi-squared Test 
table_seasonwild <- table(wild$consecutive_match, wild$season)
contingencyTableBF(table_season, sampleType = "indepMulti", fixedMargin = "rows")
#BF= 0.29, evidence in favor of independence or no association.
#the seasons do not appear to have a substantial impact on whether consecutive matches occurred (within wild)

# and only winter
winter <- od_data %>% 
  filter(season == "winter") %>% 
  mutate(combi = 1 + left + 2 * success) %>% 
  select(ID, season, status, trial, success, left, combi)  %>% 
  mutate(prev_combi = lag(combi),
         #next_combi = lead(combi),
         consecutive_match = ifelse(trial == 1, NA,
                                    ((combi %in% c(1, 4)) & (prev_combi %in% c(2, 4))) | ((combi %in% c(2, 3)) & (prev_combi %in% c(1, 3))))) %>% 
  mutate(side = case_when(
    combi %in% c(1, 4) ~ "left",
    combi %in% c(2, 3) ~ "right"))
ggplot(winter, aes(x = season, fill = consecutive_match)) +
  geom_bar(position = "dodge") +
  #facet_wrap(~ status, ncol = 1) +
  labs(title = "Consecutive Match by Season", x = "Season", y = "Count") +
  scale_fill_manual(values = c("TRUE" = "lightblue", "FALSE" = "lightgreen"))

#check if there's significance in status using Bayesian Chi-squared Test 
table_statuswinter <- table(winter$consecutive_match, winter$status)
contingencyTableBF(table_statuswinter, sampleType = "indepMulti", fixedMargin = "rows")
#status does not have a big impact either
#are the occurrence of true and false significant?

####SIDE Preference####
dat <- dat %>%
  mutate(side = case_when(
    combi %in% c(1, 4) ~ "left",
    combi %in% c(2, 3) ~ "right"))

dat$side <- as.factor(dat$side)

ggplot(dat, aes(x = season, fill = side)) +
  geom_bar(position = "dodge") +
  facet_wrap(~ status, ncol = 1) +
  labs(title = "Side Bias by Season", x = "Season", y = "Count") +
  scale_fill_manual(values = c("left" = "lightblue", "right" = "orangered"))

ggplot(dat, aes(x = ID, fill = side)) +
  geom_bar(position = "dodge") +
  labs(title = "Side Bias by ID", x = "ID", y = "Count") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(dat, aes(x = trial, fill = side)) +
  geom_bar(position = "dodge") +
  labs(title = "Side Bias by Trial", x = "Trial", y = "Count") +
  facet_wrap(~ status + season, ncol = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Again, let's look at only the wild ones:
wild <- wild %>%
  mutate(side = case_when(
    combi %in% c(1, 4) ~ "left",
    combi %in% c(2, 3) ~ "right")) %>% 
  droplevels()

wild$side <- as.factor(wild$side)

ggplot(wild, aes(x = season, fill = side)) +
  geom_bar(position = "dodge") +
  facet_wrap(~ status, ncol = 1) +
  labs(title = "Side Bias by Season", x = "Season", y = "Count") +
  scale_fill_manual(values = c("left" = "lightblue", "right" = "orangered"))

ggplot(wild, aes(x = ID, fill = side)) +
  geom_bar(position = "dodge") +
  labs(title = "Side Bias by ID", x = "ID", y = "Count") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~season)
#It seems there's some side preference for some shrews, need to do some stats on this

a <- table(wild$season, wild$side)
contingencyTableBF(a, sampleType = "indepMulti", fixedMargin = "rows")
#If BF > 1: The data provide evidence in favor of the alternative hypothesis relative to the null.
#If BF = 1: The data do not favor either hypothesis
#In such cases, the data is not providing strong evidence for or against an association.

b <- table(winter$status, winter$side)
contingencyTableBF(b, sampleType = "indepMulti", fixedMargin = "rows")

ggplot(dat, aes(x = ID, fill = side)) +
  geom_bar(position = "dodge") +
  labs(title = "Side Bias by ID", x = "ID", y = "Count") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~season)

####LATENCY and TEST_TIME####
library(splines)

wild_latency <- wild_data %>% 
  mutate(latency = as.numeric(latency)) %>% 
  mutate(latency = latency/60) %>% 
  filter(!is.na(latency) & !is.na(test_time)) %>% 
  mutate(log_test_time = log(test_time))
str(winter_latency)
winter_latency <- od_data %>%
  filter(season == "winter") %>%
  mutate(latency = as.numeric(latency)) %>% 
  mutate(latency = latency/60) %>%
  filter(!is.na(latency) & !is.na(test_time)) %>% 
  mutate(log_test_time = log(test_time)) %>% 
  droplevels()

ggplot(wild_latency, aes(x = trial_n, y = latency, fill = season)) +
  geom_bar(stat = "identity") +
  labs(title = "Latency by Trial Number",
       x = "Trial Number",
       y = "Latency") +
  theme_minimal()
ggplot(wild_latency, aes(x = trial_n, y = latency)) +
  geom_bar(stat = "identity", aes(fill = season)) +
  labs(title = "Latency by Trial Number",
       x = "Trial Number",
       y = "Latency") +
  facet_wrap(~ season, ncol = 1) + 
  theme_minimal()

#There's a lot of zeros, and few high values. I try a hurdle model:
zero_counts <- table(wild_latency$latency == 0, wild_latency$season)
total_counts <- table(wild_latency$season)
zero_proportions <- zero_counts[2, ] / total_counts

cat("Proportion of zeros in summer:", zero_proportions["summer"], "\n")
cat("Proportion of zeros in winter:", zero_proportions["winter"], "\n")

latency_wildh <- brm(bf(latency ~ 0 + bs(trial_n, degree=3) * season + (1|ID),
                        hu ~ 1 + bs(trial_n, degree=3) * season),
                     data = wild_latency, 
                     family = hurdle_lognormal())

fitsh <- fitted(latency_wildh, newdata = nwd) %>% 
  data.frame() %>% 
  bind_cols(nwd) %>% 
  mutate(#ID     = str_c("ID[", ID, "]"),
    season = factor(season),
    trial = factor(trial))

fits_avgh <- fitsh %>%
  group_by(trial, season) %>%
  summarize(AvgEstimate = mean(Estimate),
            Lower = mean(Q2.5),
            Upper = mean(Q97.5))

ggplot(wild_latency, aes(x = trial, y = latency, color = season)) +
  geom_jitter() + 
  labs(x = "Trial", y = "Latency") +
  labs(color="Season Model", fill="CI 95%", linetype="Season Smooth") +
  scale_color_manual(values = c("#D16103","#0A33A9"), labels = c("Summer", "Winter")) +
  scale_fill_manual(values = c("#D16103","#0A33A9"), labels = c("Summer", "Winter")) +
  scale_x_discrete(breaks = wild_latency$trial, labels = wild_latency$trial) +
  geom_line(data = fits_avgh, aes(x = trial, y = AvgEstimate, group = season)) +
  geom_ribbon(data = fits_avgh, aes(x = trial, ymin = Lower, ymax = Upper, fill = season, group = season), alpha= 0.2, inherit.aes = FALSE)

#same with winter

#Check for zeros:

zero_counts <- table(winter_latency$latency == 0, winter_latency$status)
total_counts <- table(winter_latency$status)
zero_proportions <- zero_counts[2, ] / total_counts
str(winter_latency$status)
cat("Proportion of zeros in wild:", zero_proportions["wild"], "\n")
cat("Proportion of zeros in captive:", zero_proportions["captive"], "\n")

#no latency=0 for captive!

latency_winterh <- brm(bf(latency ~ 0 + bs(trial_n, degree=3) * status + (1|ID),
                          hu ~ 1 + bs(trial_n, degree=3) * status),
                       data = winter_latency, 
                       family = hurdle_lognormal())
nwind <-  winter_latency %>% 
  distinct(ID, status, trial, trial_n, season, latency)
fitwinh <- fitted(latency_winterh, newdata = nwind) %>% 
  data.frame() %>% 
  bind_cols(nwind) %>% 
  mutate(#ID     = str_c("ID[", ID, "]"),
    status = factor(status),
    trial = factor(trial))

fitwin_avgh <- fitwinh %>%
  group_by(trial, status) %>%
  summarize(AvgEstimate = mean(Estimate),
            Lower = mean(Q2.5),
            Upper = mean(Q97.5))

ggplot(winter_latency, aes(x = trial, y = latency, color = status)) +
  geom_jitter() + 
  labs(x = "Trial", y = "Latency") +
  labs(color="Status Model", fill="CI 95%", linetype="Status Smooth") +
  scale_color_manual(values = c("#0AA9A9","#0A33A9"), labels = c("Captive", "Wild")) +
  scale_fill_manual(values = c("#0AA9A9","#0A33A9"), labels = c("Captive", "Wild")) +
  scale_x_discrete(breaks = winter_latency$trial, labels = winter_latency$trial) +
  geom_line(data = fitwin_avgh, aes(x = trial, y = AvgEstimate, group = status)) +
  geom_ribbon(data = fitwin_avgh, aes(x = trial, ymin = Lower, ymax = Upper, fill = status, group = status), alpha= 0.2, inherit.aes = FALSE)

####TEST TIME ####
testime_wild <- brm(log_test_time ~ 0 +bs(trial_n, degree=3) * season + (1|ID), 
                    data = wild_latency, 
                    family = Gamma(link="log"))

nwdtt <-  wild_latency %>% 
  distinct(ID, status, trial, trial_n, season, log_test_time)

fit_tt <- fitted(testime_wild, newdata = nwdtt) %>% 
  data.frame() %>% 
  bind_cols(nwdtt) %>% 
  mutate(#ID     = str_c("ID[", ID, "]"),
    season = factor(season),
    trial = factor(trial))

fitt_avg <- fit_tt %>%
  group_by(trial, season) %>%
  summarize(AvgEstimate = mean(Estimate),
            Lower = mean(Q2.5),
            Upper = mean(Q97.5))

ggplot(wild_latency, aes(x = trial, y = log_test_time, color = season)) +
  geom_jitter() + 
  labs(x = "Trial", y = "Test Time") +
  labs(color="Season Model", fill="CI 95%", linetype="Season Smooth") +
  scale_color_manual(values = c("#D16103","#0A33A9"), labels = c("Summer", "Winter")) +
  scale_fill_manual(values = c("#D16103","#0A33A9"), labels = c("Summer", "Winter")) +
  scale_x_discrete(breaks = wild_latency$trial, labels = wild_latency$trial) +
  geom_line(data = fitt_avg, aes(x = trial, y = AvgEstimate, group = season)) +
  geom_ribbon(data = fitt_avg, aes(x = trial, ymin = Lower, ymax = Upper, fill = season, group = season), alpha= 0.2, inherit.aes = FALSE)

#same per winter
testime_winter <- brm(log_test_time ~ 0 +bs(trial_n, degree=3) * status + (1|ID), 
                      data = winter_latency, 
                      family = Gamma(link="log"))
#no strong evidence to suggest a significant difference in test_time between captive and wild statuses
nwindtt <-  winter_latency %>% 
  distinct(ID, status, trial, trial_n, season, log_test_time)

fit_wintt <- fitted(testime_winter, newdata = nwindtt) %>% 
  data.frame() %>% 
  bind_cols(nwindtt) %>% 
  mutate(#ID     = str_c("ID[", ID, "]"),
    status = factor(status),
    trial = factor(trial))

fittwin_avg <- fit_wintt %>%
  group_by(trial, status) %>%
  summarize(AvgEstimate = mean(Estimate),
            Lower = mean(Q2.5),
            Upper = mean(Q97.5))

ggplot(winter_latency, aes(x = trial, y = log_test_time, color = status)) +
  geom_jitter() + 
  labs(x = "Trial", y = "Test Time") +
  labs(color="Status Model", fill="CI 95%", linetype="Status Smooth") +
  scale_color_manual(values = c("#0AA9A9","#0A33A9"), labels = c("Captive", "Wild")) +
  scale_fill_manual(values = c("#0AA9A9","#0A33A9"), labels = c("Captive", "Wild")) +
  scale_x_discrete(breaks = winter_latency$trial, labels = winter_latency$trial) +
  geom_line(data = fittwin_avg, aes(x = trial, y = AvgEstimate, group = status)) +
  geom_ribbon(data = fittwin_avg, aes(x = trial, ymin = Lower, ymax = Upper, fill = status, group = status), alpha= 0.2, inherit.aes = FALSE)

####OTHERS####

#COD, or change of direction. we coded everytime the shrew turned 180 degrees and went back. cod =1 he switched one time before choosing a box, cod = 0 he went straight to the box (correct or incorrect.)
#is the number of time they switch direction related to the probability of success?
wild <- od_data %>%
  filter(status == "wild" & cod <= 10) %>%
  droplevels()
winter <- od_data %>% 
  filter(season == "winter" & cod <= 10) %>% 
  droplevels()
cod_plotwinter <- ggplot(winter, aes(x = trial, y = cod, color = factor(success))) +
  geom_point() +
  geom_smooth(mapping = aes(x = trial, y = cod, color = factor(success), group = factor(success)), se = FALSE) +
  labs(title = "COD Counts vs. Trial", x = "Trial", y = "COD Counts", color = "Success") +
  theme_minimal()
cod_plotwild <- ggplot(wild, aes(x = trial, y = cod, color = factor(success))) +
  geom_point() +
  geom_smooth(mapping = aes(x = trial, y = cod, color = factor(success), group = factor(success)), se = FALSE) +
  labs(title = "COD Counts vs. Trial", x = "Trial", y = "COD Counts", color = "Success") +
  theme_minimal()
grid.arrange(cod_plotwinter, cod_plotwild, ncol = 2)

fit_winter <- brm(data = winter,
                  family = bernoulli,
                  success ~ 0 + cod:status + latency:status + test_time:status + sitting:status,
                  prior = c(prior(normal(0, 1.5), class = b)
                            ),
                  chains = 4, iter = 2000,
                  save_pars = save_pars(all = TRUE),
                  backend = "cmdstanr",
                  threads = threading(2),
                  control=list(adapt_delta=0.9, max_treedepth = 10))
fit_winter

fit_wild <- brm(data = wild,
                family = bernoulli,
                success ~ 0  + cod:season + test_time:season + sitting:season,
                prior = c(prior(normal(0, 1), class = b)),
                chains = 4, iter = 2000,
                save_pars = save_pars(all = TRUE),
                backend = "cmdstanr",
                threads = threading(2),
                control=list(adapt_delta=0.8, max_treedepth = 10))

conditional_effects(cod_wild)

results=summary(fit1a_winter)
results=data.frame(results$fixed)
results$covariate=rownames(results)
posteriors=ggplot(results,aes(x=Estimate,y=covariate,color=I("blue")))+geom_point()+geom_linerange(aes(xmin=l.95..CI,xmax=u.95..CI,color=I("blue")))+
  geom_vline(xintercept = 0)+theme_classic()
posteriors


#### SMOOTH Difference ####

unique_trial_n <- unique(wild_data$trial_n)
pred_data <- data.frame(trial_n = rep(unique_trial_n, 2), 
                        season = c(rep("summer", length(unique_trial_n)), 
                                   rep("winter", length(unique_trial_n))))
linpred_matrix <- posterior_linpred(fit3, newdata = pred_data, re_formula = NA, transform = TRUE)
pred_summer_samples <- linpred_matrix[, 1:length(unique_trial_n)]
pred_winter_samples <- linpred_matrix[, (length(unique_trial_n) + 1):(2 * length(unique_trial_n))]
smooth_diff_samples <- pred_summer_samples - pred_winter_samples
hist(rowMeans(smooth_diff_samples), breaks=50, main="Distribution of Differences", xlab="Difference")
mean_diff_samples <- mean(rowMeans(smooth_diff_samples))
HPDI <- quantile(rowMeans(smooth_diff_samples), c(0.025, 0.975))
print(mean_diff_samples)
#For the season model (mean_diff_samples), the mean difference is 0.0654. However, in the plot, the difference changes sign around trial 6. 
#This means that one curve ( summer) is higher in the first six trials, and the other curve ( winter) becomes higher for subsequent trials.

print(HPDI)

med_diff_samples <- apply(smooth_diff_samples, 2, median)
lwr_samples <- apply(smooth_diff_samples, 2, quantile, probs = 0.025)
upr_samples <- apply(smooth_diff_samples, 2, quantile, probs = 0.975)

plot(unique_trial_n, med_diff_samples, type = "l", col = "red", xlab = "trial_n", ylab = "Difference", ylim = range(c(lwr_samples, upr_samples)))
polygon(c(unique_trial_n, rev(unique_trial_n)), c(lwr_samples, rev(upr_samples)), col = alpha("skyblue", 0.5), border = NA)
lines(unique_trial_n, med_diff_samples, col = "red")

### and winter

unique_trial_w <- unique(winter_data$trial_n)
pred_data_w <- data.frame(trial_n = rep(unique_trial_w, 2), 
                        status = c(rep("wild", length(unique_trial_w)), 
                                   rep("captive", length(unique_trial_w))))
linpred_matrix_w <- posterior_linpred(fit2_winter, newdata = pred_data_w, re_formula = NA, transform = TRUE)
pred_wild_samples <- linpred_matrix_w[, 1:length(unique_trial_w)]
pred_captive_samples <- linpred_matrix_w[, (length(unique_trial_w) + 1):(2 * length(unique_trial_w))]
smooth_diff_samples_w <- pred_wild_samples - pred_captive_samples

hist(rowMeans(smooth_diff_samples_w), breaks=50, main="Distribution of Differences", xlab="Difference")

mean_diff_samples_w <- mean(rowMeans(smooth_diff_samples_w))
HPDI_w <- quantile(rowMeans(smooth_diff_samples_w), c(0.025, 0.975))
print(mean_diff_samples_w)
#the mean difference is 0.0943. Since this difference is always positive in the plot, it suggests that one of the curves 
#(wild) is consistently higher than the other (captive) across all trials.
print(HPDI_w)

med_diff_samples_w <- apply(smooth_diff_samples_w, 2, median)
lwr_samples_w <- apply(smooth_diff_samples_w, 2, quantile, probs = 0.025)
upr_samples_w <- apply(smooth_diff_samples_w, 2, quantile, probs = 0.975)

plot(unique_trial_w, med_diff_samples_w, type = "l", col = "red", xlab = "trial_n", ylab = "Difference", ylim = range(c(lwr_samples_w, upr_samples_w)))
polygon(c(unique_trial_w, rev(unique_trial_w)), c(lwr_samples_w, rev(upr_samples_w)), col = alpha("skyblue", 0.5), border = NA)
lines(unique_trial_w, med_diff_samples_w, col = "red")

mean_abs_diff_season <- mean(abs(smooth_diff))


mean_abs_diff_status <- mean(abs(smooth_diff_samples_w))

#The mean_abs_diff_status is 0.1987, which is considerably higher than the mean_abs_diff_season of 0.0673. 
#This is a crucial point: the absolute difference takes into account the magnitude of the difference without 
#considering the direction. 
#A higher absolute difference for the status model suggests that the shapes of the wild and captive smooths 
#deviate from each other more than the summer and winter smooths.

