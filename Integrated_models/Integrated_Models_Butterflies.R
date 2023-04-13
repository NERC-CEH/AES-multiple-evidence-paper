#' ## Butterfly integrated models

library(dplyr)
library(tidyr)
library(ggplot2)
theme_set(theme_classic())
library(brms)


my_col <- unname(palette.colors()[c(9,8)])

#' Get data ####
#' Environmental gradients represented as PCA scores, sourced from file on P
#' drive.
dir <- config::get()
scpath <- dir$directories$scoredata
pcapath <- dir$directories$pcadata
modpath <- dir$directories$models

#' PCA scores
PCA <- read.csv(paste0(pcapath,"PCA scores for CS and LandSpAES squares.csv")) %>%
  select(-X)

#' AES scores
AES <- read.csv(paste0(scpath, "Butts_Gradient_Scores.csv"))


#' source individual scheme models

source("UKBMS/Calculate_Responses_UKBMS.R")

source("WCBS/Calculate_Responses_WCBS.R")

source("LandSpaes/collate butterfly data.R")

source("LandSpAES/functions_script.R")

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

#' Combine all UKBMS data
UKBMS_all_data <- UKBMS_RESPONSES %>%
  inner_join(PCA, by = c("buttsurv.GRIDREF_1km" = "GRIDREF","YEAR")) %>%
  left_join(ukbms_aes, by = c("buttsurv.GRIDREF_1km" = "CELLCODE","YEAR" = "Year")) %>%
  ungroup() %>%
  filter(AES1KM < 50000) %>%
  mutate(YR = as.factor(YEAR),
         AES1KM = (AES1KM-3000)/5000,
         AES3KM = (AES3KM-3000)/5000,
         TRANSECT_LENGTH_NEW = TRANSECT_LENGTH_NEW/1000,
         N_VISITS_MAYTOAUGUST = N_VISITS_MAYTOAUGUST/10,
         YRnm = YEAR - 2016,
         YRnm2 = as.numeric(as.factor(YEAR)))

#' WCBS
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

#' Combine all WCBS data and scale
WCBS_all_data <- wcbs_responses %>%
  inner_join(PCA, by = c("buttsurv.GRIDREF_1km" = "GRIDREF","YEAR")) %>%
  inner_join(wcbs_aes, by = c("buttsurv.GRIDREF_1km" = "CELLCODE","YEAR")) %>%
  ungroup() %>% 
  filter(AES1KM < 50000) %>%
  mutate(YR = as.factor(YEAR),
         AES1KM = (AES1KM-3000)/5000,
         AES3KM = (AES3KM-3000)/5000,
         TRANSECT_LENGTH = TRANSECT_LENGTH/1000,
         N_VISITS_MAYTOAUGUST = N_VISITS_MAYTOAUGUST/10,
         YRnm = YEAR - 2016,
         YRnm2 = as.numeric(as.factor(YEAR)))

# LandSpAES
buttabund <- aggregate(BUTTERFLY_COUNT ~ NCA + SURVEY_SQUARE + SURVEY_YEAR + AES1KM + AES3KM, data = buttcount2, FUN = function(x) sum(x))
#

buttcount2$RICHNESS_ID <- buttcount2$BUTTERFLY_ID

#remove "other butterflies" id

buttcount2$RICHNESS_ID[buttcount2$RICHNESS_ID == 65] <- NA


buttrich <- aggregate(RICHNESS_ID ~ NCA + SURVEY_SQUARE + SURVEY_YEAR + AES1KM + AES3KM, data = buttcount2, FUN = function(x) length(unique(x[!is.na(x)])), na.action = na.pass)
#

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

butt_vars <- merge(buttabund, buttrich, by = c("NCA","SURVEY_SQUARE","SURVEY_YEAR","AES1KM","AES3KM"))
butt_vars <- merge(butt_vars, buttdiv, by = c("NCA","SURVEY_SQUARE","SURVEY_YEAR","AES1KM","AES3KM"))

buttrounds <- aggregate(ROUND_NUMBER ~ SURVEY_SQUARE + SURVEY_YEAR, data = buttcount2, FUN = function(x) length(unique(x)))

butt_vars <- merge(butt_vars, buttrounds, by.x = c("SURVEY_SQUARE", "SURVEY_YEAR"), by.y = c("SURVEY_SQUARE", "SURVEY_YEAR"))

#sunshine
buttsun <- aggregate(PERC_SUN ~ SURVEY_SQUARE + SURVEY_YEAR, data = buttvariables, FUN = function(x) mean(x, na.rm = TRUE))

butt_vars <- merge(butt_vars, buttsun, by.x = c("SURVEY_SQUARE", "SURVEY_YEAR"), by.y = c("SURVEY_SQUARE", "SURVEY_YEAR"))

#temperature
butttemp <- aggregate(SURVEY_TEMP_SHADE ~ SURVEY_SQUARE + SURVEY_YEAR, data = buttvisit, FUN = function(x) mean(x, na.rm = TRUE))

butt_vars <- merge(butt_vars, butttemp, by.x = c("SURVEY_SQUARE", "SURVEY_YEAR"), by.y = c("SURVEY_SQUARE", "SURVEY_YEAR"))

#PCA scores

#need to remove half square for matching to PCA scores
butt_vars$SURVEY_SQUARE <- sapply(strsplit(buttabund$SURVEY_SQUARE, "\\/"),function(x) x[1])

butt_vars <- merge(butt_vars, PCA, by.x = c("SURVEY_SQUARE", "SURVEY_YEAR"), by.y = c("GRIDREF", "YEAR"))


#' Rescale predictors and fit model
butt_vars$AES1KM <- (butt_vars$AES1KM-3000)/5000
butt_vars$AES3KM <- (butt_vars$AES3KM-3000)/5000
butt_vars$ROUND_NUMBER <- butt_vars$ROUND_NUMBER/10

#treat survey year as factor
butt_vars$SURVEY_YEAR <- factor(butt_vars$SURVEY_YEAR)

#
all_data <- butt_vars %>%
  dplyr::rename(N_VISITS_MAYTOAUGUST = ROUND_NUMBER, YR = SURVEY_YEAR,
         Abundance = BUTTERFLY_COUNT, Richness = RICHNESS_ID, 
         Shannon_diversity = SHANNON_DIV, SITENO = SURVEY_SQUARE) %>%
  mutate(SURVEY = "LandSpAES", TRANSECT_LENGTH_NEW = 2) %>%
  full_join(mutate(UKBMS_all_data, SURVEY = "UKBMS", 
                   SITENO = as.character(SITENO))) %>%
  full_join(mutate(WCBS_all_data, SURVEY = "WCBS", TRANSECT_LENGTH_NEW = 2,
                   SITENO = as.character(SITENO)))

#version without UKBMS for Abundance mod
all_data2 <- butt_vars %>%
  dplyr::rename(N_VISITS_MAYTOAUGUST = ROUND_NUMBER, YR = SURVEY_YEAR,
                Abundance = BUTTERFLY_COUNT, Richness = RICHNESS_ID, 
                Shannon_diversity = SHANNON_DIV, SITENO = SURVEY_SQUARE) %>%
  mutate(SURVEY = "LandSpAES", TRANSECT_LENGTH_NEW = 2) %>%
  full_join(mutate(WCBS_all_data, SURVEY = "WCBS", TRANSECT_LENGTH_NEW = 2,
                   SITENO = as.character(SITENO)))

#transform Shannon diversity

all_data$expShannon_diversity <- exp(all_data$Shannon_diversity)

##fit model in glmmTMB

library(glmmTMB)
abund_mod <- glmmTMB(Abundance ~ AES1KM * AES3KM + Climate_PC1 + Habitat_PC1 + Landscape_PC1 + 
                     YR + N_VISITS_MAYTOAUGUST + 
                     (1|SITENO) + (1|SURVEY), data = all_data2, family = "nbinom2")


rich_mod <- glmmTMB(Richness ~ AES1KM * AES3KM + Climate_PC1 + Habitat_PC1 + Landscape_PC1 + 
                       YR + N_VISITS_MAYTOAUGUST + TRANSECT_LENGTH_NEW +
                       (1|SITENO) + (1|SURVEY), data = all_data, family = "poisson")


div_mod <- glmmTMB(expShannon_diversity ~ AES1KM * AES3KM + Climate_PC1 + Habitat_PC1 + Landscape_PC1 + 
                      YR + N_VISITS_MAYTOAUGUST + TRANSECT_LENGTH_NEW +
                      (1|SITENO) + (1|SURVEY), data = all_data, family = "gaussian")


#save model outputs

saveRDS(abund_mod, file = "Integrated_Abundance_model.RDS")
saveRDS(rich_mod, file = "Integrated_Richness_model.RDS")
saveRDS(div_mod, file = "Integrated_Diversity_model.RDS")




##read in LandSpAES models for plotting


Rich_LS_mod <- readRDS(paste0(modpath, "LandSpAES_Richness_brm.RDS"))

Abund_LS_mod <- readRDS(paste0(modpath, "LandSpAES_Abundance_brm.RDS"))

Div_LS_mod <- readRDS(paste0(modpath, "LandSpAES_Diversity_brm.RDS"))






##plot abundance

p1km <- ggpredict(abund_mod, "AES1KM[-0.60:8.21, by = 0.5]",
                condition = c("AES3KM" = "0.122",
                              "N_VISITS_MAYTOAUGUST" = "0.4",
                              "YR" = "2018",
                              "Climate_PC1" = "-0.04",
                              "Landscape_PC1" = "-0.05",
                              "Habitat_PC1" = "0.26"))
p1 <- ggpredict(Abund_LS_mod, "AES1KM[-0.60:7.03, by = 0.5]",
                condition = c("AES3KM" = "0.122",
                              "ROUND_NUMBER" = "0.4",
                              "SURVEY_YEAR" = "2018",
                              "Climate_PC1" = "-0.04",
                              "Landscape_PC1" = "-0.05",
                              "Habitat_PC1" = "0.26"))
p1km$group <- "Integrated"
p1$group <- "LandSpAES"

p3km <- ggpredict(abund_mod, "AES3KM[-0.60:7.38, by = 0.5]",
                condition = c("AES1KM" = "0.362",
                              "N_VISITS_MAYTOAUGUST" = "0.4",
                              "YR" = "2018",
                              "Climate_PC1" = "-0.04",
                              "Landscape_PC1" = "-0.05",
                              "Habitat_PC1" = "0.26"))
p3 <- ggpredict(Abund_LS_mod, "AES3KM[-0.60:2.42, by = 0.5]",
                condition = c("AES1KM" = "0.362",
                              "ROUND_NUMBER" = "0.4",
                              "SURVEY_YEAR" = "2018",
                              "Climate_PC1" = "-0.04",
                              "Landscape_PC1" = "-0.05",
                              "Habitat_PC1" = "0.26"))
p3km$group <- "Integrated"
p3$group <- "LandSpAES"

p1km$std.error <- NULL
p3km$std.error <- NULL
p1kmc <- do.call(rbind, list(p1km, p1)) %>%
  plyr::mutate(x = (x*5000 + 3000)/1000)
p3kmc <- do.call(rbind, list(p3km, p3)) %>%
  plyr::mutate(x = (x*5000 + 3000)/1000)

abund_1km <- ggplot(p1kmc, aes(x = x, y = predicted, colour = group, fill = group)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.3, colour = NA) +
  geom_line() +
  coord_cartesian(clip = "off") +
  scale_fill_manual(aesthetics = c("fill","colour"),
                    values = my_col,
                    name = "Survey") +
  scale_y_continuous(limits = c(0,1100), expand = c(0,0)) +
  labs(x = "AES 1km ('000s)", y = "Predicted Butterfly Abundance") +
  NULL
abund_1km

ggsave(paste0(modpath,"Combined plots/Butterfly abundance 1km integrated.png"), height = 800, width = 1000, units = "mm", scale = 0.15)


abund_3km <- ggplot(p3kmc, aes(x = x, y = predicted, colour = group, fill = group)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.3, colour = NA) +
  geom_line() +
  coord_cartesian(clip = "off") +
  scale_fill_manual(aesthetics = c("fill","colour"),
                    values = my_col,
                    name = "Survey") +
  scale_y_continuous(limits = c(0,1100), expand = c(0,0)) +
  labs(x = "AES 3km ('000s)", y = "Predicted Butterfly Abundance") +
  NULL
abund_3km

ggsave(paste0(modpath,"Combined plots/Butterfly abundance 3km integrated.png"), height = 800, width = 1000, units = "mm", scale = 0.15)

abund_1km + abund_3km + plot_layout(guides='collect') &
  theme(legend.position='bottom')
ggsave(paste0(modpath,"Combined plots/Butterfly abundance combined 1km and 3km integrated.png"), height = 800, width = 1200, units = "mm", scale = 0.15)


##plot richness

p1km <- ggpredict(rich_mod, "AES1KM[-0.60:8.21, by = 0.5]",
                  condition = c("AES3KM" = "0.122",
                                "N_VISITS_MAYTOAUGUST" = "0.4",
                                "YR" = "2018",
                                "TRANSECT_LENGTH_NEW" = "2",
                                "Climate_PC1" = "-0.04",
                                "Landscape_PC1" = "-0.05",
                                "Habitat_PC1" = "0.26"))
p1 <- ggpredict(Rich_LS_mod, "AES1KM[-0.60:7.03, by = 0.5]",
                condition = c("AES3KM" = "0.122",
                              "ROUND_NUMBER" = "0.4",
                              "SURVEY_YEAR" = "2018",
                              "Climate_PC1" = "-0.04",
                              "Landscape_PC1" = "-0.05",
                              "Habitat_PC1" = "0.26"))

p3km <- ggpredict(rich_mod, "AES3KM[-0.60:7.38, by = 0.5]",
                  condition = c("AES1KM" = "0.362",
                                "N_VISITS_MAYTOAUGUST" = "0.4",
                                "YR" = "2018",
                                "TRANSECT_LENGTH_NEW" = "2",
                                "Climate_PC1" = "-0.04",
                                "Landscape_PC1" = "-0.05",
                                "Habitat_PC1" = "0.26"))
p3 <- ggpredict(Rich_LS_mod, "AES3KM[-0.60:2.42, by = 0.5]",
                condition = c("AES1KM" = "0.362",
                              "ROUND_NUMBER" = "0.4",
                              "SURVEY_YEAR" = "2018",
                              "Climate_PC1" = "-0.04",
                              "Landscape_PC1" = "-0.05",
                              "Habitat_PC1" = "0.26"))

p1km$group <- "Integrated"
p1$group <- "LandSpAES"
p3km$group <- "Integrated"
p3$group <- "LandSpAES"

p1km$std.error <- NULL
p3km$std.error <- NULL
p1kmc <- do.call(rbind, list(p1km, p1)) %>%
  plyr::mutate(x = (x*5000 + 3000)/1000)
p3kmc <- do.call(rbind, list(p3km, p3)) %>%
  plyr::mutate(x = (x*5000 + 3000)/1000)


rich_1km <- ggplot(p1kmc, aes(x = x, y = predicted, colour = group, fill = group)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.3, colour = NA) +
  geom_line() +
  coord_cartesian(clip = "off") +
  scale_fill_manual(aesthetics = c("fill","colour"),
                    values = my_col,
                    name = "Survey") +
  scale_y_continuous(limits = c(0,30), expand = c(0,0)) +
  labs(x = "AES 1km ('000s)", y = "Predicted Butterfly Richness") +
  NULL
rich_1km

ggsave(paste0(modpath,"Combined plots/Butterfly richness 1km integrated.png"), height = 800, width = 1000, units = "mm", scale = 0.15)


rich_3km <- ggplot(p3kmc, aes(x = x, y = predicted, colour = group, fill = group)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.3, colour = NA) +
  geom_line() +
  coord_cartesian(clip = "off") +
  scale_fill_manual(aesthetics = c("fill","colour"),
                    values = my_col,
                    name = "Survey") +
  scale_y_continuous(limits = c(0,30), expand = c(0,0)) +
  labs(x = "AES 3km ('000s)", y = "Predicted Butterfly Richness") +
  NULL
rich_3km

ggsave(paste0(modpath,"Combined plots/Butterfly richness 3km integrated.png"), height = 800, width = 1000, units = "mm", scale = 0.15)

rich_1km + rich_3km + plot_layout(guides='collect') &
  theme(legend.position='bottom')
ggsave(paste0(modpath,"Combined plots/Butterfly richness combined 1km and 3km integrated.png"), height = 800, width = 1200, units = "mm", scale = 0.15)



##plot diversity

p1km <- ggpredict(div_mod, "AES1KM[-0.60:8.21, by = 0.5]",
                  condition = c("AES3KM" = "0.122",
                                "N_VISITS_MAYTOAUGUST" = "0.4",
                                "YR" = "2018",
                                "TRANSECT_LENGTH_NEW" = "2",
                                "Climate_PC1" = "-0.04",
                                "Landscape_PC1" = "-0.05",
                                "Habitat_PC1" = "0.26"))
p1 <- ggpredict(Div_LS_mod, "AES1KM[-0.60:7.03, by = 0.5]",
                condition = c("AES3KM" = "0.122",
                              "ROUND_NUMBER" = "0.4",
                              "SURVEY_YEAR" = "2018",
                              "Climate_PC1" = "-0.04",
                              "Landscape_PC1" = "-0.05",
                              "Habitat_PC1" = "0.26"))

p3km <- ggpredict(div_mod, "AES3KM[-0.60:7.38, by = 0.5]",
                  condition = c("AES1KM" = "0.362",
                                "N_VISITS_MAYTOAUGUST" = "0.4",
                                "YR" = "2018",
                                "TRANSECT_LENGTH_NEW" = "2",
                                "Climate_PC1" = "-0.04",
                                "Landscape_PC1" = "-0.05",
                                "Habitat_PC1" = "0.26"))
p3 <- ggpredict(Div_LS_mod, "AES3KM[-0.60:2.42, by = 0.5]",
                condition = c("AES1KM" = "0.362",
                              "ROUND_NUMBER" = "0.4",
                              "SURVEY_YEAR" = "2018",
                              "Climate_PC1" = "-0.04",
                              "Landscape_PC1" = "-0.05",
                              "Habitat_PC1" = "0.26"))

p1km$group <- "Integrated"
p1$group <- "LandSpAES"
p3km$group <- "Integrated"
p3$group <- "LandSpAES"

p1km$std.error <- NULL
p3km$std.error <- NULL
p1kmc <- do.call(rbind, list(p1km, p1)) %>%
  plyr::mutate(x = (x*5000 + 3000)/1000)
p3kmc <- do.call(rbind, list(p3km, p3)) %>%
  plyr::mutate(x = (x*5000 + 3000)/1000)


div_1km <- ggplot(p1kmc, aes(x = x, y = predicted, colour = group, fill = group)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.3, colour = NA) +
  geom_line() +
  coord_cartesian(clip = "off") +
  scale_fill_manual(aesthetics = c("fill","colour"),
                    values = my_col,
                    name = "Survey") +
  scale_y_continuous(limits = c(0,15), expand = c(0,0)) +
  labs(x = "AES 1km ('000s)", y = "Predicted Butterfly Diversity") +
  NULL
div_1km

ggsave(paste0(modpath,"Combined plots/Butterfly diversity 1km integrated.png"), height = 800, width = 1000, units = "mm", scale = 0.15)


div_3km <- ggplot(p3kmc, aes(x = x, y = predicted, colour = group, fill = group)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.3, colour = NA) +
  geom_line() +
  coord_cartesian(clip = "off") +
  scale_fill_manual(aesthetics = c("fill","colour"),
                    values = my_col,
                    name = "Survey") +
  scale_y_continuous(limits = c(0,15), expand = c(0,0)) +
  labs(x = "AES 3km ('000s)", y = "Predicted Butterfly diversity") +
  NULL
div_3km

ggsave(paste0(modpath,"Combined plots/Butterfly diversity 3km integrated.png"), height = 800, width = 1000, units = "mm", scale = 0.15)

div_1km + div_3km + plot_layout(guides='collect') &
  theme(legend.position='bottom')
ggsave(paste0(modpath,"Combined plots/Butterfly diversity combined 1km and 3km integrated.png"), height = 800, width = 1200, units = "mm", scale = 0.15)


abund_1km + abund_3km + div_1km + div_3km + rich_1km + rich_3km + plot_layout(guides = 'collect', ncol = 2) &
  theme(legend.position = "bottom")
ggsave(paste0(modpath,"Combined plots/All integrated plots combined.png"), height = 1600, width = 1100, units = "mm", scale = 0.15)




