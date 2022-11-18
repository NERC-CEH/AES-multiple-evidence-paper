#' UKBMS - Analyse each response to AES scores plus environmental gradients
#' 
library(DHARMa)
library(dplyr)
library(tidyr)
library(ggplot2)
theme_set(theme_classic())
library(brms)

#' Get data ####
#' Environmental gradients represented as PCA scores, sourced from file on P
#' drive.
dir <- config::get()
scpath <- dir$directories$scoredata
pcapath <- dir$directories$pcadata

#' PCA scores
PCA <- read.csv(paste0(pcapath,"PCA scores for CS and LandSpAES squares.csv")) %>%
  select(-X) %>%
  distinct()

#' AES scores
AES <- read.csv(paste0(scpath, "Butts_Gradient_Scores.csv"))

#' UKBMS response data from file
source("UKBMS/Calculate_Responses_UKBMS.R")

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
  inner_join(PCA, by = c("buttsurv.GRIDREF_1km" = "PLAN_NO","YEAR")) %>%
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
                        N_VISITS_MAYTOAUGUST + TRANSECT_LENGTH_NEW + YR +
                        Climate_PC1 + Landscape_PC1 + Habitat_PC1 +
                        (1|SITENO),
                      data = UKBMS_all_data, family = "poisson", prior = mod_pr,
                      cores = 4)
summary(Rich_ukbms_mod)
plot(Rich_ukbms_mod)
pp_check(Rich_ukbms_mod)
# way out, overpredicting at low counts, underpredicting at medium counts and then
# overpredicting at highest counts
plot(conditional_effects(Rich_ukbms_mod, effects = "AES1KM:AES3KM",
                         int_conditions = list(AES3KM = c(0.1,0.25,0.75))),
     rug = TRUE, theme = ggplot2::theme_classic())

mod_pr <- prior(normal(0,1), class = b) +
  prior(student_t(3,0,1), class = sd)
Rich_ukbms_mod2 <- brm(Richness ~ AES1KM*AES3KM + 
                        N_VISITS_MAYTOAUGUST + TRANSECT_LENGTH_NEW + YR +
                        Climate_PC1 + Landscape_PC1 + Habitat_PC1 +
                        (1|SITENO),
                      data = UKBMS_all_data, family = "negbinomial", prior = mod_pr,
                      cores = 4)
summary(Rich_ukbms_mod2)
pp_check(Rich_ukbms_mod2)
# no improvement in fit to data, shape parameter 1452 +/- 279


#' ### Abundance models
# poisson models 

mod_pr <- prior(normal(0,1), class = b) +
  prior(student_t(3,0,1), class = sd)
Abun_ukbms_mod <- brm(Abundance ~ AES1KM*AES3KM + 
                        N_VISITS_MAYTOAUGUST + TRANSECT_LENGTH_NEW + YR +
                        Climate_PC1 + Landscape_PC1 + Habitat_PC1 +
                        (1|SITENO),
                      data = UKBMS_all_data, family = "poisson", prior = mod_pr,
                      cores = 4)
# does not converge - site effect Rhat 1.46, plus various fixed effects Rhat > 1.2
# negative binomial - with site RE
Abun_ukbms_mod <- brm(Abundance ~ AES1KM*AES3KM + 
                        N_VISITS_MAYTOAUGUST + TRANSECT_LENGTH_NEW + YR +
                        Climate_PC1 + Landscape_PC1 + Habitat_PC1 +
                        (1|SITENO),
                      data = UKBMS_all_data, family = "negbinomial", prior = mod_pr,
                      cores = 4)
summary(Abun_ukbms_mod)
plot(Abun_ukbms_mod)
pp_check(Abun_ukbms_mod)
pp_check(Abun_ukbms_mod) +
  scale_x_continuous(limits = c(0,10000))
# reasonable recovery of data 
pp_check(Abun_ukbms_mod, type = "ecdf_overlay") + 
  scale_x_continuous(limits = c(0,2000))
plot(conditional_effects(Abun_ukbms_mod, effects = "AES1KM:AES3KM",
                         int_conditions = list(AES3KM = c(0.1,0.25,0.75))),
     rug = TRUE, theme = ggplot2::theme_classic())



#' ### Diversity models
mod_pr <- prior(normal(0,1), class = b) +
  prior(student_t(3,0,1), class = sd)
Div_ukbms_mod <- brm(Shannon_diversity ~ AES1KM*AES3KM + 
                       N_VISITS_MAYTOAUGUST + TRANSECT_LENGTH_NEW + YR +
                       Climate_PC1 + Landscape_PC1 + Habitat_PC1 +
                       (1|SITENO),
                     data = UKBMS_all_data, prior = mod_pr,
                     cores = 4)
summary(Div_ukbms_mod)
plot(Div_ukbms_mod)
pp_check(Div_ukbms_mod)
plot(conditional_effects(Div_ukbms_mod, effects = "AES1KM:AES3KM",
                         int_conditions = list(AES3KM = c(0.1,0.25,0.75))),
     rug = TRUE, theme = ggplot2::theme_classic())

UKBMS_all_data$expDiversity <- exp(UKBMS_all_data$Shannon_diversity)
Div_ukbms_mod <- brm(expDiversity ~ AES1KM*AES3KM + 
                       N_VISITS_MAYTOAUGUST + TRANSECT_LENGTH_NEW + YR +
                       Climate_PC1 + Landscape_PC1 + Habitat_PC1 +
                       (1|SITENO),
                     data = UKBMS_all_data, prior = mod_pr,
                     cores = 4)
summary(Div_ukbms_mod)
plot(Div_ukbms_mod)
pp_check(Div_ukbms_mod) # much better fit to data
plot(conditional_effects(Div_ukbms_mod, effects = "AES1KM:AES3KM",
                         int_conditions = list(AES3KM = c(0.1,0.25,0.75))),
     rug = TRUE, theme = ggplot2::theme_classic())
