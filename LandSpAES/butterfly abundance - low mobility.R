##butterfly trait abundance models


homedir <- getwd()

analysis_group <- "Butterfly trait analyses"
analysis_name <- "Butterfly abundance - low mobility"
plot_name <- "Total count of low mobility butterflies"



library(lme4)
library(effects)
library(tidyverse)
library(DHARMa)
library(dplyr)
library(tidyr)
library(ggplot2)
theme_set(theme_classic())
library(brms)



# load file paths
dir <- config::get()
pcapath <- dir$directories$pcadata
modpath <- dir$directories$models

# source data
source("LandSpAES/collate butterfly data.R")

# source function
source("LandSpAES/functions_script.R")

#PCA scores

PCA <- read.csv(paste0(pcapath,"PCA scores for CS and LandSpAES squares.csv")) %>%
  dplyr::select(-X)


#Select trait grouping


db <- config::get("AES")

con <- odbcConnect(db$DSN, db$uid, db$pwd, believeNRows = FALSE)

butttraits <- sqlFetch(con, "TBL_BUTTERFLY_SPECIES")

close(con)


buttcount2$TRAIT <- butttraits$WINGSPAN_CATEGORY[match(buttcount2$BUTTERFLY_ID, butttraits$BUTTERFLY_ID)]

buttcount2$BUTTERFLY_COUNT <- ifelse(buttcount2$TRAIT == 1, buttcount2$BUTTERFLY_COUNT, NA)




#' **Calculate BUTTERFLY abundance - aggregated**

buttabund <- aggregate(BUTTERFLY_COUNT ~ NCA + SURVEY_SQUARE + SURVEY_YEAR + AES1KM + AES3KM, data = buttcount2, FUN = function(x) sum(x))
#90 observations

buttrounds <- aggregate(ROUND_NUMBER ~ SURVEY_SQUARE + SURVEY_YEAR, data = buttcount2, FUN = function(x) length(unique(x)))

buttabund <- merge(buttabund, buttrounds, by.x = c("SURVEY_SQUARE", "SURVEY_YEAR"), by.y = c("SURVEY_SQUARE", "SURVEY_YEAR"))

#sunshine
buttsun <- aggregate(PERC_SUN ~ SURVEY_SQUARE + SURVEY_YEAR, data = buttvariables, FUN = function(x) mean(x, na.rm = TRUE))

buttabund <- merge(buttabund, buttsun, by.x = c("SURVEY_SQUARE", "SURVEY_YEAR"), by.y = c("SURVEY_SQUARE", "SURVEY_YEAR"))

#temperature
butttemp <- aggregate(SURVEY_TEMP_SHADE ~ SURVEY_SQUARE + SURVEY_YEAR, data = buttvisit, FUN = function(x) mean(x, na.rm = TRUE))

buttabund <- merge(buttabund, butttemp, by.x = c("SURVEY_SQUARE", "SURVEY_YEAR"), by.y = c("SURVEY_SQUARE", "SURVEY_YEAR"))

#PCA scores

#need to remove half square for matching to PCA scores
buttabund$SURVEY_SQUARE <- sapply(strsplit(buttabund$SURVEY_SQUARE, "\\/"),function(x) x[1])

buttabund <- merge(buttabund, PCA, by.x = c("SURVEY_SQUARE", "SURVEY_YEAR"), by.y = c("GRIDREF", "YEAR"))


#' Look at distribution of abundance variable

fit_distrib(buttabund, 'BUTTERFLY_COUNT')

# fits negative binomial fairly well

#' Rescale predictors and fit model
#' 
buttabund$AES1KM <- (buttabund$AES1KM-3000)/5000
buttabund$AES3KM <- (buttabund$AES3KM-3000)/5000
buttabund$ROUND_NUMBER <- buttabund$ROUND_NUMBER/10

#treat survey year as factor
buttabund$SURVEY_YEAR <- factor(buttabund$SURVEY_YEAR)


# models
mod1 <- glmer.nb(BUTTERFLY_COUNT ~ AES1KM*AES3KM + ROUND_NUMBER +  SURVEY_YEAR + Climate_PC1 + Landscape_PC1 + Habitat_PC1 + (1|SURVEY_SQUARE), data = buttabund)

summary(mod1)


mod_pr <- prior(normal(0,1), class = b) +
  prior(student_t(3,0,1), class = sd)
Abun_LS_mod <- brm(BUTTERFLY_COUNT ~ AES1KM*AES3KM +
                     ROUND_NUMBER + SURVEY_YEAR + 
                     Climate_PC1 + Landscape_PC1 + Habitat_PC1 +
                     (1|SURVEY_SQUARE),
                   data = buttabund, family = "negbinomial", prior = mod_pr,
                   cores = 4, file = paste0(modpath,"LandSpAES_Lowmob_Abundance_brm"))
summary(Abun_LS_mod)


plot(Abun_LS_mod)
pp_check(Abun_LS_mod)
pp_check(Abun_LS_mod) +
  scale_x_continuous(limits = c(0,1000))
pp_check(Abun_LS_mod, type = "ecdf_overlay") + 
  scale_x_continuous(limits = c(0,1000))
plot(conditional_effects(Abun_LS_mod, effects = "AES1KM:AES3KM",
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
ggsave("LandSpAES Low mobility Abundance AES1km AES3km interaction.png",
       path = modpath, width = 15, height = 12, units = "cm", dpi = 600)

plot(conditional_effects(Abun_LS_mod, effects = "AES3KM:AES1KM",
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
ggsave("LandSpAES Low mobility Abundance AES3km AES1km interaction.png",
       path = modpath, width = 15, height = 12, units = "cm", dpi = 600)













