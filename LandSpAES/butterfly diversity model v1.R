## Butterfly diversity models


analysis_group <- "Headline analyses"
analysis_name <- "Butterfly diversity"
plot_name <- "Butterfly Shannon diversity index"


library(reshape)
library(reshape2)
library(vegan)
library(nlme)
library(tidyverse)
library(effects)
library(DHARMa)
library(dplyr)
library(tidyr)
library(ggplot2)
theme_set(theme_classic())
library(brms)


# folder setup for saving
dir <- config::get()
pcapath <- dir$directories$pcadata
modpath <- dir$directories$models

# data collation and aggregation
source("LandSpAES/collate butterfly data.R")

# functions
source("LandSpAES/functions_script.R")



PCA <- read.csv(paste0(pcapath,"PCA scores for CS and LandSpAES squares.csv")) %>%
  dplyr::select(-X)





##Now aggregate data for models



buttcount2$RICHNESS_ID <- buttcount2$BUTTERFLY_ID

#remove "other butterflies" id

buttcount2$RICHNESS_ID[buttcount2$RICHNESS_ID == 65] <- NA


length(buttcount2$RICHNESS_ID[is.na(buttcount2$RICHNESS_ID)])
# 1253/14266 no butterflies



#calculate vector of richness species
richness_bb <- unique(buttcount2$RICHNESS_ID[!is.na(buttcount2$RICHNESS_ID)])

#need unique combination of square and year...

buttcount2$SQ_YR <- paste(buttcount2$SURVEY_SQUARE, buttcount2$SURVEY_YEAR, sep = "_")

#create a species by site table for just the species with SPECIES_RICHNESS_SPECIES entries aggregated at the square level
c1 <- cast(buttcount2, SQ_YR ~ RICHNESS_ID, value = "BUTTERFLY_COUNT", fun.aggregate = "sum", fill = 0, subset = RICHNESS_ID %in% richness_bb)

#use the vegan package to calculate the Shannon diversity index
c1$shann <- diversity(c1[,2:35], index = "shannon")


#' Add the Shannon index to the buttcount2 table (note all observations in the same square will get the same entry)
#' 
buttcount2$SHANNON_DIV <- c1$shann[match(buttcount2$SQ_YR, c1$SQ_YR)]

#' Aggregate across rounds and transect section lengths taking the mean Shannon diversity index (can compare with values from c1 to check this is correct)
buttdiv <- aggregate(SHANNON_DIV ~ NCA + SURVEY_SQUARE + SURVEY_YEAR + AES1KM + AES3KM, data = buttcount2, FUN = function(x) mean(x), na.action = na.pass)

buttrounds <- aggregate(ROUND_NUMBER ~ SURVEY_SQUARE + SURVEY_YEAR, data = buttcount2, FUN = function(x) length(unique(x)))

buttdiv <- merge(buttdiv, buttrounds, by.x = c("SURVEY_SQUARE", "SURVEY_YEAR"), by.y = c("SURVEY_SQUARE", "SURVEY_YEAR"))

#sunshine
buttsun <- aggregate(PERC_SUN ~ SURVEY_SQUARE + SURVEY_YEAR, data = buttvariables, FUN = function(x) mean(x, na.rm = TRUE))

buttdiv <- merge(buttdiv, buttsun, by.x = c("SURVEY_SQUARE", "SURVEY_YEAR"), by.y = c("SURVEY_SQUARE", "SURVEY_YEAR"))

#temperature
butttemp <- aggregate(SURVEY_TEMP_SHADE ~ SURVEY_SQUARE + SURVEY_YEAR, data = buttvisit, FUN = function(x) mean(x, na.rm = TRUE))

buttdiv <- merge(buttdiv, butttemp, by.x = c("SURVEY_SQUARE", "SURVEY_YEAR"), by.y = c("SURVEY_SQUARE", "SURVEY_YEAR"))


#PCA scores

#need to remove half square for matching to PCA scores
buttdiv$SURVEY_SQUARE <- sapply(strsplit(buttdiv$SURVEY_SQUARE, "\\/"),function(x) x[1])

buttdiv <- merge(buttdiv, PCA, by.x = c("SURVEY_SQUARE", "SURVEY_YEAR"), by.y = c("GRIDREF", "YEAR"))



#check diversity is normally distributed
hist(exp(buttdiv$SHANNON_DIV))





#scale AES variables
buttdiv$AES1KM <- scale(buttdiv$AES1KM)
buttdiv$AES3KM <- scale(buttdiv$AES3KM)

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
                         int_conditions = list(AES3KM = c(0.05,0.25,0.75))),
     rug = TRUE, theme = ggplot2::theme_classic(),
     rug_args = list(colour = "gray"))[[1]] +
  scale_x_continuous(breaks = c(0,0.5,1,1.5,2.0),
                     labels = c(0,10,20,30,40),
                     expand = c(0,0)) +
  scale_colour_manual(values = c("#E69F00","#CC79A7","#0072B2"),
                      aesthetics = c("colour","fill"),
                      labels = c(7500,2500,500)) +
  labs(x = "AES1KM ('000s)", y = "exp(Shannon diversity)")
ggsave("LandSpAES Diversity AES1km AES3km interaction.png",
       path = modpath, width = 15, height = 12, units = "cm", dpi = 600)


plot(conditional_effects(Div_LS_mod, effects = "AES3KM:AES1KM",
                         int_conditions = list(AES1KM = c(0.05,0.25,0.75))),
     rug = TRUE, theme = ggplot2::theme_classic(),
     rug_args = list(colour = "gray"))[[1]] +
  scale_x_continuous(breaks = c(0,0.5,1,1.5,2.0),
                     labels = c(0,10,20,30,40),
                     expand = c(0,0)) +
  scale_colour_manual(values = c("#E69F00","#CC79A7","#0072B2"),
                      aesthetics = c("colour","fill"),
                      labels = c(7500,2500,500)) +
  labs(x = "AES3KM ('000s)", y = "exp(Shannon diversity)")
ggsave("LandSpAES Diversity AES3km AES1km interaction.png",
       path = modpath, width = 15, height = 12, units = "cm", dpi = 600)

