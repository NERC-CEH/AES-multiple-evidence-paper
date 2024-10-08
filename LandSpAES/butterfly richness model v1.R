###Basic richness models

library(reshape)
library(reshape2)
library(vegan)
library(lme4)
library(effects)
library(tidyverse)
library(DHARMa)
library(dplyr)
library(tidyr)
library(ggplot2)
theme_set(theme_classic())
library(brms)


modpath <- getwd()

buttrich <- read.csv("LandSpAES butterfly richness.csv")


#scale predictors
buttrich$AES1KM <- (buttrich$AES1KM-3000)/5000
buttrich$AES3KM <- (buttrich$AES3KM-3000)/5000
buttrich$ROUND_NUMBER <- buttrich$ROUND_NUMBER/10

#treat survey year as factor
buttrich$SURVEY_YEAR <- factor(buttrich$SURVEY_YEAR)


#go for Poisson


mod1 <- glmer(RICHNESS_ID ~ AES1KM*AES3KM + ROUND_NUMBER + SURVEY_YEAR + 
                Climate_PC1 + Landscape_PC1 + Habitat_PC1 + 
                (1|SURVEY_SQUARE), data = buttrich, family = "poisson")

summary(mod1)


mod_pr <- prior(normal(0,1), class = b) +
  prior(student_t(3,0,1), class = sd)
Rich_LS_mod <- brm(RICHNESS_ID ~ AES1KM*AES3KM + 
                     ROUND_NUMBER + SURVEY_YEAR +
                       Climate_PC1 + Landscape_PC1 + Habitat_PC1 +
                       (1|SURVEY_SQUARE),
                     data = buttrich, family = "poisson", prior = mod_pr,
                     cores = 4, file = paste0(modpath, "LandSpAES_Richness_brm"))
summary(Rich_LS_mod)
plot(Rich_LS_mod)
pp_check(Rich_LS_mod)
# not too far away, but still slightly underpredicting at low counts and
# overpredicting at medium counts
plot(conditional_effects(Rich_LS_mod, effects = "AES1KM:AES3KM",
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
ggsave("LandSpAES Richness AES1km AES3km interaction.png",
       path = modpath, width = 15, height = 12, units = "cm", dpi = 600)

plot(conditional_effects(Rich_LS_mod, effects = "AES3KM:AES1KM",
                         int_conditions = list(AES1KM = c(-0.5, -0.1, 0.4))),
     rug = TRUE, theme = ggplot2::theme_classic(),
     rug_args = list(colour = "gray"))[[1]] +
  scale_x_continuous(breaks = c(-0.6,1.4,3.4,5.4),
                     labels = c(0,10,20,30),
                     expand = c(0,0)) +
  scale_colour_manual(values = c("#E69F00","#CC79A7","#0072B2"),
                      aesthetics = c("colour","fill"),
                      labels = c(5000,2500,500)) +
  labs(x = "AES3KM ('000s)")
ggsave("LandSpAES Richness AES3km AES1km interaction.png",
       path = modpath, width = 15, height = 12, units = "cm", dpi = 600)