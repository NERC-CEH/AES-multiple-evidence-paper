#' UKBMS - Analyse each response to AES scores plus environmental gradients
#' 
library(MASS)
library(DHARMa)
library(car)
library(dplyr)
library(tidyr)
library(ggplot2)
theme_set(theme_classic())
library(brms)

#' Get data ####
#' Environmental gradients represented as PCA scores, sourced from file on P
#' drive.

scpath <- dir$directories$scoredata

#' AES scores
AES <- read.csv(paste0(scpath, "Butts_Gradient_Scores.csv"))

#' UKBMS response data from file
source("Calculate_Responses_UKBMS.R")

#' Manipulate data ####
#' Get AES scores for UKBMS
ukbms_aes <- AES %>%
  filter(CELLCODE %in% UKBMS_RESPONSES$buttsurv.GRIDREF_1km) %>%
  filter(MaskStatus == "OK") %>%
  dplyr::select(CELLCODE, starts_with("Sc")) %>%
  tidyr::pivot_longer(starts_with("Sc"),
                      names_to = c("Scale","Year"),
                      names_sep = "_",
                      names_transform = list(Year = as.numeric)) %>%
  mutate(Scale = dplyr::recode(Scale,
                               "Score1km" = "AES1KM",
                               "Score3km" = "AES3KM")) %>%
  distinct() %>%
  tidyr::pivot_wider(names_from = Scale,
                     values_from = value)

#' Combine all data
UKBMS_all_data <- UKBMS_RESPONSES %>%
  left_join(ukbms_aes, by = c("buttsurv.GRIDREF_1km" = "CELLCODE","YEAR" = "Year")) %>%
  ungroup() %>%
  filter(AES1KM < 50000) %>%
  mutate(YR = as.factor(YEAR),
         AES1KM = AES1KM/20000,
         AES3KM = AES3KM/10000,
         TRANSECT_LENGTH_NEW = TRANSECT_LENGTH_NEW/1000,
         N_VISITS_MAYTOAUGUST = N_VISITS_MAYTOAUGUST/10,
         YRnm = YEAR - 2016,
         YRnm2 = as.numeric(as.factor(YEAR)))
summary(UKBMS_all_data)
# psych::multi.hist(select_if(UKBMS_all_data, is.numeric))
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
Rich_ukbms_mod <- brm(Richness ~ AES1KM*AES3KM + 
                        N_VISITS_MAYTOAUGUST + TRANSECT_LENGTH_NEW +
                        (1|SITENO),
                      data = UKBMS_all_data, family = "poisson", prior = mod_pr,
                      cores = 4)
summary(Rich_ukbms_mod)
plot(Rich_ukbms_mod)
pp_check(Rich_ukbms_mod)
# way out, overpredicting at low counts, underpredicting at medium counts and then
# underpredicting at highest counts

# Temp AR
mod_pr <- prior(normal(0,1), class = b) +
  prior(uniform(-1,1), lb = -1, ub = 1, class = ar)
Rich_ukbms_mod_temp <- brm(Richness ~ AES1KM*AES3KM + 
                             N_VISITS_MAYTOAUGUST + TRANSECT_LENGTH_NEW +
                             ar(time = YRnm, gr = SITENO),
                           data = UKBMS_all_data, family = "poisson", prior = mod_pr,
                           cores = 4)
summary(Rich_ukbms_mod_temp)
plot(Rich_ukbms_mod_temp)
pp_check(Rich_ukbms_mod_temp)
# Still not fitting data properly

mod_pr <- prior(normal(0,1), class = b) +
  prior(student_t(3,0,1), class = sd) +
  prior(uniform(-1,1), lb = -1, ub = 1, class = ar)
Rich_ukbms_mod_temp2 <- brm(Richness ~ AES1KM*AES3KM + 
                              N_VISITS_MAYTOAUGUST + TRANSECT_LENGTH_NEW +
                              ar(time = YRnm, gr = SITENO) + (1|SITENO),
                            data = UKBMS_all_data, family = "poisson", prior = mod_pr,
                            cores = 4)
summary(Rich_ukbms_mod_temp2)
plot(Rich_ukbms_mod_temp2)
pp_check(Rich_ukbms_mod_temp2)
# If include the SITE random effect then the model cannot resolve the temporal
# autocorrelation - have to pick one or the other. Data still isn't being fit
# properly

# Check model residuals with DHARMa
model.check <- createDHARMa(
  simulatedResponse = t(posterior_predict(Rich_ukbms_mod_temp)),
  observedResponse = UKBMS_all_data$Richness,
  fittedPredictedResponse = apply(t(posterior_epred(Rich_ukbms_mod_temp)), 1, mean),
  integerResponse = TRUE)

plot(model.check) # highly problematic

# Negative binomial models
mod_pr <- prior(normal(0,1), class = b) +
  prior(uniform(-1,1), lb = -1, ub = 1, class = ar) +
  prior(gamma(0.01,0.01), class = shape)
Rich_ukbms_mod_temp3 <- brm(Richness ~ AES1KM*AES3KM + 
                              N_VISITS_MAYTOAUGUST + TRANSECT_LENGTH_NEW +
                              ar(time = YRnm, gr = SITENO),
                            data = UKBMS_all_data, family = "negbinomial", 
                            prior = mod_pr,
                            cores = 4)
summary(Rich_ukbms_mod_temp3)
plot(Rich_ukbms_mod_temp3)
pp_check(Rich_ukbms_mod_temp3)

model.check <- createDHARMa(
  simulatedResponse = t(posterior_predict(Rich_ukbms_mod_temp3)),
  observedResponse = UKBMS_all_data$Richness,
  fittedPredictedResponse = apply(t(posterior_epred(Rich_ukbms_mod_temp3)), 1, mean),
  integerResponse = TRUE)

plot(model.check)
# still bad - KS, Dispersion and outlier test all failed

mod_pr <- prior(normal(0,1), class = b) +
  prior(student_t(3,0,1), class = sd) +
  prior(gamma(0.01,0.01), class = shape)
Rich_ukbms_mod_temp4 <- brm(Richness ~ AES1KM*AES3KM + 
                              N_VISITS_MAYTOAUGUST + TRANSECT_LENGTH_NEW +
                              (1|SITENO),
                            data = UKBMS_all_data, family = "negbinomial", 
                            prior = mod_pr,
                            cores = 4)
summary(Rich_ukbms_mod_temp4)
plot(Rich_ukbms_mod_temp4)
pp_check(Rich_ukbms_mod_temp4)

model.check <- createDHARMa(
  simulatedResponse = t(posterior_predict(Rich_ukbms_mod_temp4)),
  observedResponse = UKBMS_all_data$Richness,
  fittedPredictedResponse = apply(t(posterior_epred(Rich_ukbms_mod_temp4)), 1, mean),
  integerResponse = TRUE)

plot(model.check)
# Even worse!

# Revisit this when the PCA axes are included - there's a problem with not predicting
# the low numbers in some squares, which might be related to environmental conditions


#' ### Abundance models
# poisson models 

mod_pr <- prior(normal(0,1), class = b) +
  prior(student_t(3,0,1), class = sd)
Abun_ukbms_mod <- brm(Abundance ~ AES1KM*AES3KM + 
                        N_VISITS_MAYTOAUGUST + TRANSECT_LENGTH_NEW +
                        (1|SITENO),
                      data = UKBMS_all_data, family = "poisson", prior = mod_pr,
                      cores = 4)
# does not converge
# negative binomial - with site RE
Abun_ukbms_mod <- brm(Abundance ~ AES1KM*AES3KM + 
                        N_VISITS_MAYTOAUGUST + TRANSECT_LENGTH_NEW +
                        (1|SITENO),
                      data = UKBMS_all_data, family = "negbinomial", prior = mod_pr,
                      cores = 4)
summary(Abun_ukbms_mod)
plot(Abun_ukbms_mod)
pp_check(Abun_ukbms_mod)
pp_check(Abun_ukbms_mod) +
  scale_x_continuous(limits = c(0,10000))
# reasonable recovery of data - worth revisiting when env covariates included
pp_check(Abun_ukbms_mod, type = "ecdf_overlay") + 
  scale_x_continuous(limits = c(0,2000))

# negative binomial with temp AR
mod_pr <- prior(normal(0,1), class = b) +
  prior(uniform(-1,1), lb = -1, ub = 1, class = ar)
Abun_ukbms_mod_temp <- brm(Abundance ~ AES1KM*AES3KM + 
                             N_VISITS_MAYTOAUGUST + TRANSECT_LENGTH_NEW +
                             ar(time = YRnm, gr = SITENO),
                           data = UKBMS_all_data, family = "negbinomial", 
                           prior = mod_pr,
                           cores = 4)
# Rhat = 1.84 - shape parameter, ar[1] and sderr also unconverged. If fit with
# poisson then still doesn't converge, and all parameters bad, not just the shape
# parameter
summary(Abun_ukbms_mod_temp)
plot(Abun_ukbms_mod_temp)
pp_check(Abun_ukbms_mod_temp)

# Negative binomial with site RE and temp AR
mod_pr <- prior(normal(0,1), class = b) +
  prior(student_t(3,0,1), class = sd) +
  prior(uniform(-1,1), lb = -1, ub = 1, class = ar)
Abun_ukbms_mod_temp2 <- brm(Abundance ~ AES1KM*AES3KM + 
                              N_VISITS_MAYTOAUGUST + TRANSECT_LENGTH_NEW +
                              ar(time = YRnm, gr = SITENO) + (1|SITENO),
                            data = UKBMS_all_data, family = "negbinomial", 
                            prior = mod_pr,
                            cores = 4)
# 14 divergent transitions and autoregressive part of model is not being identified
summary(Abun_ukbms_mod_temp2)
plot(Abun_ukbms_mod_temp2)
pp_check(Abun_ukbms_mod_temp2)
# check divergent transitions
pairs(Abun_ukbms_mod_temp2$fit, 
      pars = c("sd_SITENO__Intercept","shape","ar[1]",
               "sderr","b_Intercept"), 
      log = TRUE)
# all divergent transitions at high sderr and ar[1] around 0.8

# May be possible to get rid of them through modifying adapt_delta but model is
# performing so badly at resolving the autocorrelation that I think it may not be
# worth trying that.


#' ### Diversity models
mod_pr <- prior(normal(0,1), class = b) +
  prior(student_t(3,0,1), class = sd)
Div_ukbms_mod <- brm(Shannon_diversity ~ AES1KM*AES3KM + 
                       N_VISITS_MAYTOAUGUST + TRANSECT_LENGTH_NEW +
                       (1|SITENO),
                     data = UKBMS_all_data, prior = mod_pr,
                     cores = 4)
summary(Div_ukbms_mod)
plot(Div_ukbms_mod)
pp_check(Div_ukbms_mod)


mod_pr <- prior(normal(0,1), class = b) +
  prior(uniform(-1,1), lb = -1, ub = 1, class = ar)
Div_ukbms_mod_temp <- brm(Shannon_diversity ~ AES1KM*AES3KM + 
                            N_VISITS_MAYTOAUGUST + TRANSECT_LENGTH_NEW +
                            ar(time = YRnm, gr = SITENO),
                          data = UKBMS_all_data, prior = mod_pr,
                          cores = 4)
summary(Div_ukbms_mod_temp)
plot(Div_ukbms_mod_temp)
pp_check(Div_ukbms_mod_temp)


mod_pr <- prior(normal(0,1), class = b) +
  prior(student_t(3,0,1), class = sd) +
  prior(uniform(-1,1), lb = -1, ub = 1, class = ar)
Div_ukbms_mod_temp2 <- brm(Shannon_diversity ~ AES1KM*AES3KM + 
                             N_VISITS_MAYTOAUGUST + TRANSECT_LENGTH_NEW +
                             ar(time = YRnm, gr = SITENO) + (1|SITENO),
                           data = UKBMS_all_data, prior = mod_pr,
                           cores = 4)
summary(Div_ukbms_mod_temp2)
plot(Div_ukbms_mod_temp2)
pp_check(Div_ukbms_mod_temp2)
# ar[1] resolved to be 0, and sd_SITENO is the same as the random effect on site
# number without the temporal autocorrelation. All 3 models run fine and resolve the
# data to the same level of accuracy