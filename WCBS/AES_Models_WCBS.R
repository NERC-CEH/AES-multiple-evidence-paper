#' WCBS - Analyse each response to AES scores plus environmental gradients
#' 
library(DHARMa)
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
  select(-X)

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
  inner_join(PCA, by = c("buttsurv.GRIDREF_1km" = "GRIDREF","YEAR")) %>%
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
                       N_VISITS_MAYTOAUGUST + YR +
                       Climate_PC1 + Landscape_PC1 + Habitat_PC1 +
                       (1|SITENO),
                     data = WCBS_all_data, family = "poisson", prior = mod_pr,
                     cores = 4)
summary(Rich_wcbs_mod)
plot(Rich_wcbs_mod)
pp_check(Rich_wcbs_mod)
# not too far away, but still slightly underpredicting at low counts and
# overpredicting at medium counts
plot(conditional_effects(Rich_wcbs_mod, effects = "AES1KM:AES3KM",
                         int_conditions = list(AES3KM = c(0.1,0.25,0.75))),
     rug = TRUE, theme = ggplot2::theme_classic())



# Check model residuals with DHARMa
model.check <- createDHARMa(
  simulatedResponse = t(posterior_predict(Rich_wcbs_mod)),
  observedResponse = WCBS_all_data$Richness,
  fittedPredictedResponse = apply(t(posterior_epred(Rich_wcbs_mod)), 1, mean),
  integerResponse = TRUE)

plot(model.check) # highly problematic for KS and dispersion

mod_pr <- prior(normal(0,1), class = b) +
  prior(student_t(3,0,1), class = sd)
Rich_wcbs_mod2 <- brm(Richness ~ AES1KM*AES3KM + 
                       N_VISITS_MAYTOAUGUST + YR +
                       Climate_PC1 + Landscape_PC1 + Habitat_PC1 +
                       (1|SITENO),
                     data = WCBS_all_data, family = "negbinomial", prior = mod_pr,
                     cores = 4)
summary(Rich_wcbs_mod2)
plot(Rich_wcbs_mod2)
pp_check(Rich_wcbs_mod2)
# no improvement in model fit to data, shape parameter is estimated to be 611 +/- 174

#' ### Abundance models
# poisson models 

mod_pr <- prior(normal(0,1), class = b) +
  prior(student_t(3,0,1), class = sd)
Abun_wcbs_mod <- brm(Abundance ~ AES1KM*AES3KM +
                        N_VISITS_MAYTOAUGUST + YR + 
                       Climate_PC1 + Landscape_PC1 + Habitat_PC1 +
                       (1|SITENO),
                     data = WCBS_all_data, family = "poisson", prior = mod_pr,
                     cores = 4)
summary(Abun_wcbs_mod)
# does not converge - Rhat of intercept and group-level effects too high
# negative binomial - with site RE
Abun_wcbs_mod <- brm(Abundance ~ AES1KM*AES3KM +
                       N_VISITS_MAYTOAUGUST + YR + 
                       Climate_PC1 + Landscape_PC1 + Habitat_PC1 +
                       (1|SITENO),
                     data = WCBS_all_data, family = "negbinomial", prior = mod_pr,
                     cores = 4)
summary(Abun_wcbs_mod)
plot(Abun_wcbs_mod)
pp_check(Abun_wcbs_mod)
pp_check(Abun_wcbs_mod) +
  scale_x_continuous(limits = c(0,1000))
pp_check(Abun_wcbs_mod, type = "ecdf_overlay") + 
  scale_x_continuous(limits = c(0,1000))
plot(conditional_effects(Abun_wcbs_mod, effects = "AES1KM:AES3KM",
                         int_conditions = list(AES3KM = c(0.1,0.25,0.75))),
     rug = TRUE, theme = ggplot2::theme_classic())


#' ### Diversity models
mod_pr <- prior(normal(0,1), class = b) +
  prior(student_t(3,0,1), class = sd)
Div_wcbs_mod <- brm(Shannon_diversity ~ AES1KM*AES3KM + 
                       N_VISITS_MAYTOAUGUST + YR +
                      Climate_PC1 + Landscape_PC1 + Habitat_PC1 +
                       (1|SITENO),
                     data = WCBS_all_data, prior = mod_pr,
                     cores = 4)
summary(Div_wcbs_mod)
plot(Div_wcbs_mod)
pp_check(Div_wcbs_mod)
# poorly fitted

WCBS_all_data$expDiversity <- exp(WCBS_all_data$Shannon_diversity)
Div_wcbs_mod <- brm(expDiversity ~ AES1KM*AES3KM + 
                      N_VISITS_MAYTOAUGUST + YR +
                      Climate_PC1 + Landscape_PC1 + Habitat_PC1 +
                      (1|SITENO),
                    data = WCBS_all_data, prior = mod_pr,
                    cores = 4)
summary(Div_wcbs_mod)
plot(Div_wcbs_mod)
pp_check(Div_wcbs_mod)
# good fit to data
plot(conditional_effects(Div_wcbs_mod, effects = "AES1KM:AES3KM",
                         int_conditions = list(AES3KM = c(0.1,0.25,0.75))),
     rug = TRUE, theme = ggplot2::theme_classic())
