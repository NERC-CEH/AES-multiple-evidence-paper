## Butterfly diversity models


library(reshape)
library(reshape2)
library(vegan)
library(nlme)
library(tidyverse)
library(effects)
library(dplyr)
library(tidyr)
library(ggplot2)
theme_set(theme_classic())
library(brms)


modpath <- getwd()

buttdiv <- read.csv("LandSpAES butterfly diversity.csv")

#check diversity is normally distributed
hist(exp(buttdiv$SHANNON_DIV))


#scale AES variables
buttdiv$AES1KM <- (buttdiv$AES1KM-3000)/5000
buttdiv$AES3KM <- (buttdiv$AES3KM-3000)/5000
buttdiv$ROUND_NUMBER <- buttdiv$ROUND_NUMBER/10

#treat survey year as factor
buttdiv$SURVEY_YEAR <- factor(buttdiv$SURVEY_YEAR)


boxplot(buttdiv$SHANNON_DIV ~ buttdiv$ROUND_NUMBER)

#model with the nlme package not lme4 as the response is normally distributed
mod1 <- lme(SHANNON_DIV ~ AES1KM*AES3KM + ROUND_NUMBER + SURVEY_YEAR + 
              Climate_PC1 + Landscape_PC1 + Habitat_PC1, 
            random = ~1|SURVEY_SQUARE, data = buttdiv)
summary(mod1)



mod_pr <- prior(normal(0,1), class = b) +
  prior(student_t(3,0,1), class = sd)
Div_LS_mod <- brm(SHANNON_DIV ~ AES1KM*AES3KM + 
                      ROUND_NUMBER + SURVEY_YEAR +
                      Climate_PC1 + Landscape_PC1 + Habitat_PC1 +
                      (1|SURVEY_SQUARE),
                    data = buttdiv, prior = mod_pr,
                    cores = 4)
summary(Div_LS_mod)
plot(Div_LS_mod)
pp_check(Div_LS_mod)
# fits fine - useful to use same transformation as other datasets

buttdiv$expSHANNON_DIV <- exp(buttdiv$SHANNON_DIV)
Div_LS_mod <- brm(expSHANNON_DIV ~ AES1KM*AES3KM + 
                    ROUND_NUMBER + SURVEY_YEAR +
                      Climate_PC1 + Landscape_PC1 + Habitat_PC1 +
                      (1|SURVEY_SQUARE),
                    data = buttdiv, prior = mod_pr,
                    cores = 4, file = paste0(modpath,"LandSpAES_Diversity_brm"))
summary(Div_LS_mod)
plot(Div_LS_mod)
pp_check(Div_LS_mod)
# good fit to data
plot(conditional_effects(Div_LS_mod, effects = "AES1KM:AES3KM",
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
ggsave("LandSpAES Diversity AES1km AES3km interaction.png",
       path = modpath, width = 15, height = 12, units = "cm", dpi = 600)

plot(conditional_effects(Div_LS_mod, effects = "AES3KM:AES1KM",
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
ggsave("LandSpAES Diversity AES3km AES1km interaction.png",
       path = modpath, width = 15, height = 12, units = "cm", dpi = 600)

