#' WCBS - Analyse each response to AES scores plus environmental gradients
#' 
library(dplyr)
library(tidyr)
library(ggplot2)
theme_set(theme_classic())
library(brms)


modpath <- getwd()


WCBS_all_data <- read.csv("WCBS butterfly data.csv")

WCBS_all_data <- WCBS_all_data %>%
  mutate(YR = as.factor(YEAR),
         AES1KM = (AES1KM-3000)/5000,
         AES3KM = (AES3KM-3000)/5000,
         TRANSECT_LENGTH = TRANSECT_LENGTH/1000,
         N_VISITS_MAYTOAUGUST = N_VISITS_MAYTOAUGUST/10,
         YRnm = YEAR - 2016,
         YRnm2 = as.numeric(as.factor(YEAR)))

#' Richness models ####
#'
#' First we fit the models using a poisson response for richness. 
#' Poisson models: 
# Site RE
mod_pr <- prior(normal(0,1), class = b) +
  prior(student_t(3,0,1), class = sd)
Rich_wcbs_mod <- brm(Richness ~ AES1KM*AES3KM + 
                       N_VISITS_MAYTOAUGUST + YR +
                       Climate_PC1 + Landscape_PC1 + Habitat_PC1 +
                       (1|SITENO),
                     data = WCBS_all_data, family = "poisson", prior = mod_pr,
                     cores = 4, file = paste0(modpath, "WCBS_Richness_brm"))
summary(Rich_wcbs_mod)
plot(Rich_wcbs_mod)
pp_check(Rich_wcbs_mod)
# not too far away, but still slightly underpredicting at low counts and
# overpredicting at medium counts
plot(conditional_effects(Rich_wcbs_mod, effects = "AES1KM:AES3KM",
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
ggsave("WCBS Richness AES1km AES3km interaction.png",
       path = modpath, width = 15, height = 12, units = "cm", dpi = 600)

plot(conditional_effects(Rich_wcbs_mod, effects = "AES3KM:AES1KM",
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
ggsave("WCBS Richness AES3km AES1km interaction.png",
       path = modpath, width = 15, height = 12, units = "cm", dpi = 600)


# mod_pr <- prior(normal(0,1), class = b) +
#   prior(student_t(3,0,1), class = sd)
# Rich_wcbs_mod2 <- brm(Richness ~ AES1KM*AES3KM + 
#                        N_VISITS_MAYTOAUGUST + YR +
#                        Climate_PC1 + Landscape_PC1 + Habitat_PC1 +
#                        (1|SITENO),
#                      data = WCBS_all_data, family = "negbinomial", prior = mod_pr,
#                      cores = 4)
# summary(Rich_wcbs_mod2)
# plot(Rich_wcbs_mod2)
# pp_check(Rich_wcbs_mod2)
# no improvement in model fit to data, shape parameter is estimated to be 614 +/- 179

#' ### Abundance models
# poisson models 

mod_pr <- prior(normal(0,1), class = b) +
  prior(student_t(3,0,1), class = sd)
# Abun_wcbs_mod <- brm(Abundance ~ AES1KM*AES3KM +
#                         N_VISITS_MAYTOAUGUST + YR + 
#                        Climate_PC1 + Landscape_PC1 + Habitat_PC1 +
#                        (1|SITENO),
#                      data = WCBS_all_data, family = "poisson", prior = mod_pr,
#                      cores = 4)
# summary(Abun_wcbs_mod)
# does not converge - Rhat of intercept and group-level effects too high
# negative binomial - with site RE
Abun_wcbs_mod <- brm(Abundance ~ AES1KM*AES3KM +
                       N_VISITS_MAYTOAUGUST + YR + 
                       Climate_PC1 + Landscape_PC1 + Habitat_PC1 +
                       (1|SITENO),
                     data = WCBS_all_data, family = "negbinomial", prior = mod_pr,
                     cores = 4, file = paste0(modpath,"WCBS_Abundance_brm"))
summary(Abun_wcbs_mod)
plot(Abun_wcbs_mod)
pp_check(Abun_wcbs_mod)
pp_check(Abun_wcbs_mod) +
  scale_x_continuous(limits = c(0,1000))
pp_check(Abun_wcbs_mod, type = "ecdf_overlay") + 
  scale_x_continuous(limits = c(0,1000))
plot(conditional_effects(Abun_wcbs_mod, effects = "AES1KM:AES3KM",
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
ggsave("WCBS Abundance AES1km AES3km interaction.png",
       path = modpath, width = 15, height = 12, units = "cm", dpi = 600)

plot(conditional_effects(Abun_wcbs_mod, effects = "AES3KM:AES1KM",
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
ggsave("WCBS Abundance AES3km AES1km interaction.png",
       path = modpath, width = 15, height = 12, units = "cm", dpi = 600)



#' ### Diversity models
# mod_pr <- prior(normal(0,1), class = b) +
#   prior(student_t(3,0,1), class = sd)
# Div_wcbs_mod <- brm(Shannon_diversity ~ AES1KM*AES3KM + 
#                        N_VISITS_MAYTOAUGUST + YR +
#                       Climate_PC1 + Landscape_PC1 + Habitat_PC1 +
#                        (1|SITENO),
#                      data = WCBS_all_data, prior = mod_pr,
#                      cores = 4)
# summary(Div_wcbs_mod)
# plot(Div_wcbs_mod)
# pp_check(Div_wcbs_mod)
# poorly fitted

WCBS_all_data$expDiversity <- exp(WCBS_all_data$Shannon_diversity)
Div_wcbs_mod <- brm(expDiversity ~ AES1KM*AES3KM + 
                      N_VISITS_MAYTOAUGUST + YR +
                      Climate_PC1 + Landscape_PC1 + Habitat_PC1 +
                      (1|SITENO),
                    data = WCBS_all_data, prior = mod_pr,
                    cores = 4, file = paste0(modpath,"WCBS_Diversity_brm"))
summary(Div_wcbs_mod)
plot(Div_wcbs_mod)
pp_check(Div_wcbs_mod)
# good fit to data
plot(conditional_effects(Div_wcbs_mod, effects = "AES1KM:AES3KM",
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
ggsave("WCBS Diversity AES1km AES3km interaction.png",
       path = modpath, width = 15, height = 12, units = "cm", dpi = 600)

plot(conditional_effects(Div_wcbs_mod, effects = "AES3KM:AES1KM",
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
ggsave("WCBS Diversity AES3km AES1km interaction.png",
       path = modpath, width = 15, height = 12, units = "cm", dpi = 600)
