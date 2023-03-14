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
modpath <- dir$directories$models

#' PCA scores
PCA <- read.csv(paste0(pcapath,"PCA scores for CS and LandSpAES squares.csv")) %>%
  select(-X)

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
summary(UKBMS_all_data)
# psych::multi.hist(select_if(UKBMS_all_data, is.numeric))
# par(mfrow=c(1,1))


#' ### Abundance models

mod_pr <- prior(normal(0,1), class = b) +
  prior(student_t(3,0,1), class = sd)


#LOW MOBILITY

# negative binomial - with site RE
Abun_ukbms_mod <- brm(Low_mobility_abund ~ AES1KM*AES3KM + 
                        N_VISITS_MAYTOAUGUST + TRANSECT_LENGTH_NEW + YR +
                        Climate_PC1 + Landscape_PC1 + Habitat_PC1 +
                        (1|SITENO),
                      data = UKBMS_all_data, family = "negbinomial", prior = mod_pr,
                      cores = 4, file = paste0(modpath,"UKBMS_LM_Abundance_brm"))
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
ggsave("UKBMS Low mobility Abundance AES1km AES3km interaction.png",
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
ggsave("UKBMS Low mobility Abundance AES3km AES1km interaction.png",
       path = modpath, width = 15, height = 12, units = "cm", dpi = 600)



#MEDIUM MOBILITY

# negative binomial - with site RE
Abun_ukbms_mod <- brm(Med_mobility_abund ~ AES1KM*AES3KM + 
                        N_VISITS_MAYTOAUGUST + TRANSECT_LENGTH_NEW + YR +
                        Climate_PC1 + Landscape_PC1 + Habitat_PC1 +
                        (1|SITENO),
                      data = UKBMS_all_data, family = "negbinomial", prior = mod_pr,
                      cores = 4, file = paste0(modpath,"UKBMS_MM_Abundance_brm"))
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
ggsave("UKBMS Medium mobility Abundance AES1km AES3km interaction.png",
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
ggsave("UKBMS Medium mobility Abundance AES3km AES1km interaction.png",
       path = modpath, width = 15, height = 12, units = "cm", dpi = 600)



# negative binomial - with site RE
Abun_ukbms_mod <- brm(High_mobility_abund ~ AES1KM*AES3KM + 
                        N_VISITS_MAYTOAUGUST + TRANSECT_LENGTH_NEW + YR +
                        Climate_PC1 + Landscape_PC1 + Habitat_PC1 +
                        (1|SITENO),
                      data = UKBMS_all_data, family = "negbinomial", prior = mod_pr,
                      cores = 4, file = paste0(modpath,"UKBMS_HM_Abundance_brm"))
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
ggsave("UKBMS High mobility Abundance AES1km AES3km interaction.png",
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
ggsave("UKBMS High mobility Abundance AES3km AES1km interaction.png",
       path = modpath, width = 15, height = 12, units = "cm", dpi = 600)

