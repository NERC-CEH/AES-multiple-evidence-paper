#' UKBMS - Analyse each response to AES scores plus environmental gradients
#' 
library(DHARMa)
library(dplyr)
library(tidyr)
library(ggplot2)
theme_set(theme_classic())
library(brms)


modpath <- getwd()


UKBMS_all_data <- read.csv("UKBMS butterfly data.csv")

UKBMS_all_data <- UKBMS_all_data %>%
  mutate(YR = as.factor(YEAR),
         AES1KM = (AES1KM-3000)/5000,
         AES3KM = (AES3KM-3000)/5000,
         TRANSECT_LENGTH_NEW = TRANSECT_LENGTH_NEW/1000,
         N_VISITS_MAYTOAUGUST = N_VISITS_MAYTOAUGUST/10,
         YRnm = YEAR - 2016,
         YRnm2 = as.numeric(as.factor(YEAR)))
summary(UKBMS_all_data)


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
                      cores = 4, file = paste0(modpath,"UKBMS_Richness_brm"))
summary(Rich_ukbms_mod)
plot(Rich_ukbms_mod)
pp_check(Rich_ukbms_mod)
# overpredicting at low counts, underpredicting at medium counts and then
# overpredicting at highest counts
plot(conditional_effects(Rich_ukbms_mod, effects = "AES1KM:AES3KM",
                         int_conditions = list(AES3KM = c(-0.5,-0.1,0.4))),
     rug = TRUE, theme = ggplot2::theme_classic(),
     rug_args = list(colour = "gray"))[[1]] +
  scale_x_continuous(breaks = c(-0.6,1.4,3.4,5.4,7.4),
                     labels = c(0,10,20,30,40),
                     expand = c(0,0)) +
  scale_colour_manual(values = c("#E69F00","#CC79A7","#0072B2"),
                       aesthetics = c("colour","fill"),
                       labels = c(5000,2500,500)) +
  labs(x = "AES1KM ('000s)")
ggsave("UKBMS Richness AES1km AES3km interaction.png",
       path = modpath, width = 15, height = 12, units = "cm", dpi = 600)

plot(conditional_effects(Rich_ukbms_mod, effects = "AES3KM:AES1KM",
                         int_conditions = list(AES1KM = c(-0.5,-0.1,0.4))),
     rug = TRUE, theme = ggplot2::theme_classic(),
     rug_args = list(colour = "gray"))[[1]] +
  scale_x_continuous(breaks = c(-0.6,1.4,3.4,5.4),
                     labels = c(0,10,20,30),
                     expand = c(0,0)) +
  scale_colour_manual(values = c("#E69F00","#CC79A7","#0072B2"),
                      aesthetics = c("colour","fill"),
                      labels = c(5000,2500,500)) +
  labs(x = "AES3KM ('000s)")
ggsave("UKBMS Richness AES3km AES1km interaction.png",
       path = modpath, width = 15, height = 12, units = "cm", dpi = 600)


# mod_pr <- prior(normal(0,1), class = b) +
#   prior(student_t(3,0,1), class = sd)
# Rich_ukbms_mod2 <- brm(Richness ~ AES1KM*AES3KM + 
#                         N_VISITS_MAYTOAUGUST + TRANSECT_LENGTH_NEW + YR +
#                         Climate_PC1 + Landscape_PC1 + Habitat_PC1 +
#                         (1|SITENO),
#                       data = UKBMS_all_data, family = "negbinomial", prior = mod_pr,
#                       cores = 4)
# summary(Rich_ukbms_mod2)
# pp_check(Rich_ukbms_mod2)
# no improvement in fit to data, shape parameter 1395 +/- 265


#' ### Abundance models
# poisson models 

# mod_pr <- prior(normal(0,1), class = b) +
#   prior(student_t(3,0,1), class = sd)
# Abun_ukbms_mod <- brm(Abundance ~ AES1KM*AES3KM + 
#                         N_VISITS_MAYTOAUGUST + TRANSECT_LENGTH_NEW + YR +
#                         Climate_PC1 + Landscape_PC1 + Habitat_PC1 +
#                         (1|SITENO),
#                       data = UKBMS_all_data, family = "poisson", prior = mod_pr,
#                       cores = 4)
# summary(Abun_ukbms_mod)
# does not converge - site effect Rhat 1.33, plus intercept 1.56 and transect length
# 1.81

# negative binomial - with site RE
Abun_ukbms_mod <- brm(Abundance ~ AES1KM*AES3KM + 
                        N_VISITS_MAYTOAUGUST + TRANSECT_LENGTH_NEW + YR +
                        Climate_PC1 + Landscape_PC1 + Habitat_PC1 +
                        (1|SITENO),
                      data = UKBMS_all_data, family = "negbinomial", prior = mod_pr,
                      cores = 4, file = paste0(modpath,"UKBMS_Abundance_brm"))
summary(Abun_ukbms_mod)
plot(Abun_ukbms_mod)
pp_check(Abun_ukbms_mod)
pp_check(Abun_ukbms_mod) +
  scale_x_continuous(limits = c(0,10000))
# reasonable recovery of data 
pp_check(Abun_ukbms_mod, type = "ecdf_overlay") + 
  scale_x_continuous(limits = c(0,2000))
plot(conditional_effects(Abun_ukbms_mod, effects = "AES1KM:AES3KM",
                         int_conditions = list(AES3KM = c(-0.5,-0.1,0.4))),
     rug = TRUE, theme = ggplot2::theme_classic(),
     rug_args = list(colour = "gray"))[[1]] +
  scale_x_continuous(breaks = c(-0.6,1.4,3.4,5.4,7.4),
                     labels = c(0,10,20,30,40),
                     expand = c(0,0)) +
  scale_colour_manual(values = c("#E69F00","#CC79A7","#0072B2"),
                      aesthetics = c("colour","fill"),
                      labels = c(5000,2500,500)) +
  labs(x = "AES1KM ('000s)")
ggsave("UKBMS Abundance AES1km AES3km interaction.png",
       path = modpath, width = 15, height = 12, units = "cm", dpi = 600)

plot(conditional_effects(Abun_ukbms_mod, effects = "AES3KM:AES1KM",
                         int_conditions = list(AES1KM = c(-0.5,-0.1,0.4))),
     rug = TRUE, theme = ggplot2::theme_classic(),
     rug_args = list(colour = "gray"))[[1]] +
  scale_x_continuous(breaks = c(-0.6,1.4,3.4,5.4),
                     labels = c(0,10,20,30),
                     expand = c(0,0)) +
  scale_colour_manual(values = c("#E69F00","#CC79A7","#0072B2"),
                      aesthetics = c("colour","fill"),
                      labels = c(5000,2500,500)) +
  labs(x = "AES3KM ('000s)")
ggsave("UKBMS Abundance AES3km AES1km interaction.png",
       path = modpath, width = 15, height = 12, units = "cm", dpi = 600)


#' ### Diversity models
mod_pr <- prior(normal(0,1), class = b) +
  prior(student_t(3,0,1), class = sd)
# Div_ukbms_mod <- brm(Shannon_diversity ~ AES1KM*AES3KM + 
#                        N_VISITS_MAYTOAUGUST + TRANSECT_LENGTH_NEW + YR +
#                        Climate_PC1 + Landscape_PC1 + Habitat_PC1 +
#                        (1|SITENO),
#                      data = UKBMS_all_data, prior = mod_pr,
#                      cores = 4)
# summary(Div_ukbms_mod)
# plot(Div_ukbms_mod)
# pp_check(Div_ukbms_mod)
# plot(conditional_effects(Div_ukbms_mod, effects = "AES1KM:AES3KM",
#                          int_conditions = list(AES3KM = c(0.1,0.25,0.75))),
#      rug = TRUE, theme = ggplot2::theme_classic())

UKBMS_all_data$expDiversity <- exp(UKBMS_all_data$Shannon_diversity)
Div_ukbms_mod <- brm(expDiversity ~ AES1KM*AES3KM + 
                       N_VISITS_MAYTOAUGUST + TRANSECT_LENGTH_NEW + YR +
                       Climate_PC1 + Landscape_PC1 + Habitat_PC1 +
                       (1|SITENO),
                     data = UKBMS_all_data, prior = mod_pr,
                     cores = 4, file = paste0(modpath,"UKBMS_Diversity_brm"))
summary(Div_ukbms_mod)
plot(Div_ukbms_mod)
pp_check(Div_ukbms_mod) # much better fit to data
plot(conditional_effects(Div_ukbms_mod, effects = "AES1KM:AES3KM",
                         int_conditions = list(AES3KM = c(-0.5,-0.1,0.4))),
     rug = TRUE, theme = ggplot2::theme_classic(),
     rug_args = list(colour = "gray"))[[1]] +
  scale_x_continuous(breaks = c(-0.6,1.4,3.4,5.4,7.4),
                     labels = c(0,10,20,30,40),
                     expand = c(0,0)) +
  scale_colour_manual(values = c("#E69F00","#CC79A7","#0072B2"),
                      aesthetics = c("colour","fill"),
                      labels = c(5000,2500,500)) +
  labs(x = "AES1KM ('000s)", y = "exp(Shannon diversity)")
ggsave("UKBMS Diversity AES1km AES3km interaction.png",
       path = modpath, width = 15, height = 12, units = "cm", dpi = 600)

plot(conditional_effects(Div_ukbms_mod, effects = "AES3KM:AES1KM",
                         int_conditions = list(AES1KM = c(-0.5,-0.1,0.4))),
     rug = TRUE, theme = ggplot2::theme_classic(),
     rug_args = list(colour = "gray"))[[1]] +
  scale_x_continuous(breaks = c(-0.6,1.4,3.4,5.4),
                     labels = c(0,10,20,30),
                     expand = c(0,0)) +
  scale_colour_manual(values = c("#E69F00","#CC79A7","#0072B2"),
                      aesthetics = c("colour","fill"),
                      labels = c(5000,2500,500)) +
  labs(x = "AES3KM ('000s)", y = "exp(Shannon diversity)")
ggsave("UKBMS Diversity AES3km AES1km interaction.png",
       path = modpath, width = 15, height = 12, units = "cm", dpi = 600)


