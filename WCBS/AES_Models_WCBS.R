#' WCBS - Analyse each response to AES scores plus environmental gradients
#' 
# library(MASS)
library(DHARMa)
# library(car)
library(dplyr)
library(tidyr)
library(ggplot2)
theme_set(theme_classic())
library(brms)

#' Get data ####
dir <- config::get()
scpath <- dir$directories$scoredata
pcapath <- dir$directories$pcadata

#' PCA scores
PCA <- read.csv(paste0(pcapath,"PCA scores for CS and LandSpAES squares.csv")) %>%
  select(-X) %>%
  distinct()

#' AES scores
AES <- read.csv(paste0(scpath, "Butts_Gradient_Scores.csv"))

#' WCBS response data from file
source("WCBS/Calculate_Responses_WCBS.R")

#' Manipulate data ####

#' Get AES scores for WCBS
wcbs_aes <- AES %>%
  filter(CELLCODE %in% wcbs_responses$buttsurv.GRIDREF_1km) %>%
  select(CELLCODE, starts_with("Sc")) %>%
  pivot_longer(starts_with("Sc"),
               names_to = c("scale","YEAR"),
               names_sep = "_",
               names_transform = list(YEAR = as.numeric)) %>%
  mutate(scale = recode(scale,
                        "Score1km" = "AES1KM",
                        "Score3km" = "AES3KM")) %>%
  distinct() %>%
  pivot_wider(names_from = scale,
              values_from = value)

#' Combine all data and scale
WCBS_all_data <- wcbs_responses %>%
  # inner_join(PCA, by = c("buttsurv.GRIDREF_1km" = "PLAN_NO","YEAR")) %>%
  inner_join(wcbs_aes, by = c("buttsurv.GRIDREF_1km" = "CELLCODE","YEAR")) %>%
  ungroup() %>% 
  filter(AES1KM < 50000) %>%
  mutate(YR = as.factor(YEAR),
         AES1KM = AES1KM/20000,
         AES3KM = AES3KM/10000,
         TRANSECT_LENGTH = TRANSECT_LENGTH/1000,
         N_VISITS_MAYTOAUGUST = N_VISITS_MAYTOAUGUST/10,
         YRnm = YEAR - 2016,
         YRnm2 = as.numeric(as.factor(YEAR)))
summary(WCBS_all_data)
# psych::multi.hist(select_if(WCBS_all_data, is.numeric))
# par(mfrow=c(1,1))

#' Richness models ####
#'
#' First we fit the models using a poisson response for richness, and then check
#' if the model assumptions (mean = variance) are satisfied using the DHARMa
#' package. 
#' Poisson models: 
# Site RE
mod_pr <- prior(normal(0,1), class = b) +
  prior(student_t(3,0,1), class = sd)
Rich_wcbs_mod <- brm(Richness ~ AES1KM*AES3KM + 
                        N_VISITS_MAYTOAUGUST + 
                        (1|SITENO),
                      data = WCBS_all_data, family = "poisson", prior = mod_pr,
                      cores = 4)
summary(Rich_wcbs_mod)
plot(Rich_wcbs_mod)
pp_check(Rich_wcbs_mod)
# not too far away, but still underpredicting at low counts and overpredicting at
# medium counts

# Temp AR
mod_pr <- prior(normal(0,1), class = b) +
  prior(uniform(-1,1), lb = -1, ub = 1, class = ar)
Rich_wcbs_mod_temp <- brm(Richness ~ AES1KM*AES3KM + 
                             N_VISITS_MAYTOAUGUST + 
                             ar(time = YRnm, gr = SITENO),
                           data = WCBS_all_data, family = "poisson", prior = mod_pr,
                           cores = 4)
summary(Rich_wcbs_mod_temp)
plot(Rich_wcbs_mod_temp)
pp_check(Rich_wcbs_mod_temp)
# not too far away, but still underpredicting at low counts and overpredicting at
# medium counts

mod_pr <- prior(normal(0,1), class = b) +
  prior(student_t(3,0,1), class = sd) +
  prior(uniform(-1,1), lb = -1, ub = 1, class = ar)
Rich_wcbs_mod_temp2 <- brm(Richness ~ AES1KM*AES3KM + 
                              N_VISITS_MAYTOAUGUST + 
                              ar(time = YRnm, gr = SITENO) + (1|SITENO),
                            data = WCBS_all_data, family = "poisson", prior = mod_pr,
                            cores = 4)
summary(Rich_wcbs_mod_temp2)
plot(Rich_wcbs_mod_temp2)
pp_check(Rich_wcbs_mod_temp2)
# If include the SITE random effect then the model cannot resolve the temporal
# autocorrelation - have to pick one or the other. No improvement to fit to data

# Check model residuals with DHARMa
model.check <- createDHARMa(
  simulatedResponse = t(posterior_predict(Rich_wcbs_mod_temp)),
  observedResponse = WCBS_all_data$Richness,
  fittedPredictedResponse = apply(t(posterior_epred(Rich_wcbs_mod_temp)), 1, mean),
  integerResponse = TRUE)

plot(model.check) # highly problematic for KS and dispersion

# Negative binomial models
mod_pr <- prior(normal(0,1), class = b) +
  prior(uniform(-1,1), lb = -1, ub = 1, class = ar) +
  prior(gamma(0.01,0.01), class = shape)
Rich_wcbs_mod_temp3 <- brm(Richness ~ AES1KM*AES3KM + 
                              N_VISITS_MAYTOAUGUST + 
                              ar(time = YRnm, gr = SITENO),
                            data = WCBS_all_data, family = "negbinomial", 
                            prior = mod_pr,
                            cores = 4)
summary(Rich_wcbs_mod_temp3)
plot(Rich_wcbs_mod_temp3)
pp_check(Rich_wcbs_mod_temp3)

model.check <- createDHARMa(
  simulatedResponse = t(posterior_predict(Rich_wcbs_mod_temp3)),
  observedResponse = WCBS_all_data$Richness,
  fittedPredictedResponse = apply(t(posterior_epred(Rich_wcbs_mod_temp3)), 1, mean),
  integerResponse = TRUE)

plot(model.check)
# still bad - KS, Dispersion and outlier test all failed

mod_pr <- prior(normal(0,1), class = b) +
  prior(student_t(3,0,1), class = sd) +
  prior(gamma(0.01,0.01), class = shape)
Rich_wcbs_mod_temp4 <- brm(Richness ~ AES1KM*AES3KM + 
                              N_VISITS_MAYTOAUGUST + 
                              (1|SITENO),
                            data = WCBS_all_data, family = "negbinomial", 
                            prior = mod_pr,
                            cores = 4)
summary(Rich_wcbs_mod_temp4)
plot(Rich_wcbs_mod_temp4)
pp_check(Rich_wcbs_mod_temp4)

model.check <- createDHARMa(
  simulatedResponse = t(posterior_predict(Rich_wcbs_mod_temp4)),
  observedResponse = WCBS_all_data$Richness,
  fittedPredictedResponse = apply(t(posterior_epred(Rich_wcbs_mod_temp4)), 1, mean),
  integerResponse = TRUE)

plot(model.check)
# Even worse!

# Revisit this when the PCA axes are included - there's a problem with not predicting
# the low numbers in some squares, which might be related to environmental conditions


#' ### Abundance models
# poisson models 

mod_pr <- prior(normal(0,1), class = b) +
  prior(student_t(3,0,1), class = sd)
Abun_wcbs_mod <- brm(Abundance ~ AES1KM*AES3KM + 
                        N_VISITS_MAYTOAUGUST + 
                        (1|SITENO),
                      data = WCBS_all_data, family = "poisson", prior = mod_pr,
                      cores = 4)
# does not converge - Rhat of intercept and group-level effects too high
# negative binomial - with site RE
Abun_wcbs_mod <- brm(Abundance ~ AES1KM*AES3KM + 
                        N_VISITS_MAYTOAUGUST + 
                        (1|SITENO),
                      data = WCBS_all_data, family = "negbinomial", prior = mod_pr,
                      cores = 4)
summary(Abun_wcbs_mod)
plot(Abun_wcbs_mod)
pp_check(Abun_wcbs_mod)
pp_check(Abun_wcbs_mod) +
  scale_x_continuous(limits = c(0,1000))
# reasonable recovery of data - worth revisiting when env covariates included
pp_check(Abun_wcbs_mod, type = "ecdf_overlay") + 
  scale_x_continuous(limits = c(0,1000))

# negative binomial with temp AR
mod_pr <- prior(normal(0,1), class = b) +
  prior(uniform(-1,1), lb = -1, ub = 1, class = ar)
Abun_wcbs_mod_temp <- brm(Abundance ~ AES1KM*AES3KM + 
                             N_VISITS_MAYTOAUGUST + 
                             ar(time = YRnm, gr = SITENO),
                           data = WCBS_all_data, family = "negbinomial", 
                           prior = mod_pr,
                           cores = 4)
# Rhat = 1.09 - shape parameter, ar[1] and sderr also unconverged. If fit with
# poisson then still doesn't converge, and all parameters bad, not just the shape
# parameter
summary(Abun_wcbs_mod_temp)

# Negative binomial with site RE and temp AR
mod_pr <- prior(normal(0,1), class = b) +
  prior(student_t(3,0,1), class = sd) +
  prior(uniform(-1,1), lb = -1, ub = 1, class = ar)
Abun_wcbs_mod_temp2 <- brm(Abundance ~ AES1KM*AES3KM + 
                              N_VISITS_MAYTOAUGUST + 
                              ar(time = YRnm, gr = SITENO) + (1|SITENO),
                            data = WCBS_all_data, family = "negbinomial", 
                            prior = mod_pr,
                            cores = 4)
# 4 divergent transitions and autoregressive part of model is not being identified
summary(Abun_wcbs_mod_temp2)
plot(Abun_wcbs_mod_temp2)
pp_check(Abun_wcbs_mod_temp2)
# check divergent transitions
pairs(Abun_wcbs_mod_temp2$fit, 
      pars = c("sd_SITENO__Intercept","shape","ar[1]",
               "sderr","b_Intercept"), 
      log = TRUE)
# all divergent transitions at high sderr and ar[1] 0.0-0.6

# May be possible to get rid of them through modifying adapt_delta but model is
# performing so badly at resolving the autocorrelation that I think it may not be
# worth trying that.


#' ### Diversity models
mod_pr <- prior(normal(0,1), class = b) +
  prior(student_t(3,0,1), class = sd)
Div_wcbs_mod <- brm(Shannon_diversity ~ AES1KM*AES3KM + 
                       N_VISITS_MAYTOAUGUST + 
                       (1|SITENO),
                     data = WCBS_all_data, prior = mod_pr,
                     cores = 4)
summary(Div_wcbs_mod)
plot(Div_wcbs_mod)
pp_check(Div_wcbs_mod)
# slight overprediction at low diversity and underprediction at high

mod_pr <- prior(normal(0,1), class = b) +
  prior(uniform(-1,1), lb = -1, ub = 1, class = ar)
Div_wcbs_mod_temp <- brm(Shannon_diversity ~ AES1KM*AES3KM + 
                            N_VISITS_MAYTOAUGUST + 
                            ar(time = YRnm, gr = SITENO),
                          data = WCBS_all_data, prior = mod_pr,
                          cores = 4)
summary(Div_wcbs_mod_temp)
plot(Div_wcbs_mod_temp)
pp_check(Div_wcbs_mod_temp)


mod_pr <- prior(normal(0,1), class = b) +
  prior(student_t(3,0,1), class = sd) +
  prior(uniform(-1,1), lb = -1, ub = 1, class = ar)
Div_wcbs_mod_temp2 <- brm(Shannon_diversity ~ AES1KM*AES3KM + 
                             N_VISITS_MAYTOAUGUST + 
                             ar(time = YRnm, gr = SITENO) + (1|SITENO),
                           data = WCBS_all_data, prior = mod_pr,
                           cores = 4)
# ESS warning
summary(Div_wcbs_mod_temp2)
plot(Div_wcbs_mod_temp2)
pp_check(Div_wcbs_mod_temp2)
# ar[1] resolved to be 0, and sd_SITENO is the same as the random effect on site
# number without the temporal autocorrelation. All 3 models resolve the data to the
# same level of accuracy
