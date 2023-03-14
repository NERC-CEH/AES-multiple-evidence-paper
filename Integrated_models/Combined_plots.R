
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
library(ggeffects)
library(patchwork)
library(gtools)

# folder setup for saving
dir <- config::get()
modpath <- dir$directories$models

my_col <- unname(palette.colors()[c(8,3,4)])


Rich_LS_mod <- readRDS(paste0(modpath, "LandSpAES_Richness_brm.RDS"))
Rich_WCBS_mod <- readRDS(paste0(modpath, "WCBS_Richness_brm.RDS"))
Rich_UKBMS_mod <- readRDS(paste0(modpath, "UKBMS_Richness_brm.RDS"))

Abund_LS_mod <- readRDS(paste0(modpath, "LandSpAES_Abundance_brm.RDS"))
Abund_WCBS_mod <- readRDS(paste0(modpath, "WCBS_Abundance_brm.RDS"))
Abund_UKBMS_mod <- readRDS(paste0(modpath, "UKBMS_Abundance_brm.RDS"))

Div_LS_mod <- readRDS(paste0(modpath, "LandSpAES_Diversity_brm.RDS"))
Div_WCBS_mod <- readRDS(paste0(modpath, "WCBS_Diversity_brm.RDS"))
Div_UKBMS_mod <- readRDS(paste0(modpath, "UKBMS_Diversity_brm.RDS"))


#result tables AES effects only (for reference in BES talk)


m1 <- summary(Rich_LS_mod)$fixed[c(2:3,11),]
m1$Model <- "Rich_LS"
m2 <- summary(Rich_UKBMS_mod)$fixed[c(2:3,12),]
m2$Model <- "Rich_UKBMS"
m3 <- summary(Rich_WCBS_mod)$fixed[c(2:3,11),]
m3$Model <- "Rich_WCBS"
m4 <- summary(Abund_LS_mod)$fixed[c(2:3,11),]
m4$Model <- "Abund_LS"
m5 <- summary(Abund_UKBMS_mod)$fixed[c(2:3,12),]
m5$Model <- "Abund_UKBMS"
m6 <- summary(Abund_WCBS_mod)$fixed[c(2:3,11),]
m6$Model <- "Abund_WCBS"
m7 <- summary(Div_LS_mod)$fixed[c(2:3,11),]
m7$Model <- "Div_LS"
m8 <- summary(Div_UKBMS_mod)$fixed[c(2:3,12),]
m8$Model <- "Div_UKBMS"
m9 <- summary(Div_WCBS_mod)$fixed[c(2:3,11),]
m9$Model <- "Div_WCBS"

write.csv(rbind(m1,m2,m3,m4,m5,m6,m7,m8,m9), "Results_summary.csv")




range(Rich_LS_mod$data$AES1KM)
range(Rich_WCBS_mod$data$AES1KM)
range(Rich_UKBMS_mod$data$AES1KM)

## add AES score histograms for BES


AES_data <- smartbind(Rich_LS_mod$data, Rich_WCBS_mod$data, Rich_UKBMS_mod$data)
AES_data$SURVEY <- c(rep("LandSpAES", nrow(Rich_LS_mod$data)),
                         rep("WCBS", nrow(Rich_WCBS_mod$data)),
                     rep("UKBMS", nrow(Rich_UKBMS_mod$data)))

ggplot(AES_data, aes(x = AES1KM, fill = SURVEY))+
  geom_density(alpha = .25)
ggplot(AES_data, aes(x = AES3KM, fill = SURVEY))+
  geom_density(alpha = .25)+
  xlab("Landscape AES score") +
  ylab("Density")+
  scale_fill_manual(values = c("#F8766D", "#619CFF", "#00BA38"))

###




mean(c(Rich_LS_mod$data$AES3KM, Rich_WCBS_mod$data$AES3KM, Rich_UKBMS_mod$data$AES3KM))

mean(c(Rich_LS_mod$data$Climate_PC1, Rich_WCBS_mod$data$Climate_PC1, Rich_UKBMS_mod$data$Climate_PC1))
mean(c(Rich_LS_mod$data$Landscape_PC1, Rich_WCBS_mod$data$Landscape_PC1, Rich_UKBMS_mod$data$Landscape_PC1))
mean(c(Rich_LS_mod$data$Habitat_PC1, Rich_WCBS_mod$data$Habitat_PC1, Rich_UKBMS_mod$data$Habitat_PC1))

#scaled value for 4 rounds = 0.4

p1 <- ggpredict(Rich_LS_mod, "AES1KM[-0.60:7.03, by = 0.5]",
                condition = c("AES3KM" = "0.122",
                              "ROUND_NUMBER" = "0.4",
                              "SURVEY_YEAR" = "2018",
                              "Climate_PC1" = "-0.04",
                              "Landscape_PC1" = "-0.05",
                              "Habitat_PC1" = "0.26"))

p2 <- ggpredict(Rich_WCBS_mod, "AES1KM[-0.60:7.97, by = 0.5]",
                condition = c("AES3KM" = "0.122",
                              "N_VISITS_MAYTOAUGUST" = "0.4",
                              "YR" = "2018",
                              "Climate_PC1" = "-0.04",
                              "Landscape_PC1" = "-0.05",
                              "Habitat_PC1" = "0.26"))

p3 <- ggpredict(Rich_UKBMS_mod, "AES1KM[-0.60:8.21, by = 0.5]",
                condition = c("AES3KM" = "0.122",
                              "N_VISITS_MAYTOAUGUST" = "0.4",
                              "YR" = "2018",
					"TRANSECT_LENGTH_NEW" = "2",
                              "Climate_PC1" = "-0.04",
                              "Landscape_PC1" = "-0.05",
                              "Habitat_PC1" = "0.26"))

p1$group <- "LandSpAES"
p2$group <- "WCBS"
p3$group <- "UKBMS"


p4 <- do.call(rbind, list(p2, p1, p3)) %>%
  mutate(x = (x*5000 + 3000)/1000)
rich_1km <- ggplot(p4, aes(x = x, y = predicted, colour = group, fill = group)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.3, colour = NA) +
  geom_line() +
  # geom_rug(data = p1_raw, aes(x = x, y = response, colour = group),
  #          sides = "b") +
  # geom_rug(data = p2_raw, aes(x = x, y = response, colour = group),
  #          sides = "t") +
  # geom_rug(data = p3_raw, aes(x = x, y = response, colour = group),
  #          sides = "t", outside = TRUE) +
  coord_cartesian(clip = "off") +
  scale_fill_manual(aesthetics = c("fill","colour"),
                    values = my_col,
                    name = "Survey") +
  scale_y_continuous(limits = c(0,30), expand = c(0,0)) +
  # scale_x_continuous(limits = c(0,74000)) +
  labs(x = "AES 1km ('000s)", y = "Predicted Butterfly Richness") +
  NULL
rich_1km

ggsave(paste0(modpath,"Combined plots/Butterfly richness 1km individual schemes.png"), height = 800, width = 1000, units = "mm", scale = 0.15)

##3km plots


range(Rich_LS_mod$data$AES3KM)
range(Rich_WCBS_mod$data$AES3KM)
range(Rich_UKBMS_mod$data$AES3KM)

mean(c(Rich_LS_mod$data$AES1KM, Rich_WCBS_mod$data$AES1KM, Rich_UKBMS_mod$data$AES1KM))


p1 <- ggpredict(Rich_LS_mod, "AES3KM[-0.60:2.42, by = 0.5]",
                condition = c("AES1KM" = "0.362",
                              "ROUND_NUMBER" = "0.4",
                              "SURVEY_YEAR" = "2018",
                              "Climate_PC1" = "-0.04",
                              "Landscape_PC1" = "-0.05",
                              "Habitat_PC1" = "0.26"))

p2 <- ggpredict(Rich_WCBS_mod, "AES3KM[-0.60:6.66, by = 0.5]",
                condition = c("AES1KM" = "0.362",
                              "N_VISITS_MAYTOAUGUST" = "0.4",
                              "YR" = "2018",
                              "Climate_PC1" = "-0.04",
                              "Landscape_PC1" = "-0.05",
                              "Habitat_PC1" = "0.26"))

p3 <- ggpredict(Rich_UKBMS_mod, "AES3KM[-0.60:7.38, by = 0.5]",
                condition = c("AES1KM" = "0.362",
                              "N_VISITS_MAYTOAUGUST" = "0.4",
                              "YR" = "2018",
					"TRANSECT_LENGTH_NEW" = "2",
                              "Climate_PC1" = "-0.04",
                              "Landscape_PC1" = "-0.05",
                              "Habitat_PC1" = "0.26"))

p1$group <- "LandSpAES"
p2$group <- "WCBS"
p3$group <- "UKBMS"


p4 <- do.call(rbind, list(p2, p1, p3)) %>%
  mutate(x = (x*5000 + 3000)/1000)
rich_3km <- ggplot(p4, aes(x = x, y = predicted, colour = group, fill = group)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.3, colour = NA) +
  geom_line() +
  # geom_rug(data = p1_raw, aes(x = x, y = response, colour = group),
  #          sides = "b") +
  # geom_rug(data = p2_raw, aes(x = x, y = response, colour = group),
  #          sides = "t") +
  # geom_rug(data = p3_raw, aes(x = x, y = response, colour = group),
  #          sides = "t", outside = TRUE) +
  coord_cartesian(clip = "off") +
  scale_fill_manual(aesthetics = c("fill","colour"),
                    values = my_col,
                    name = "Survey") +
  scale_y_continuous(limits = c(0,30), expand = c(0,0)) +
  # scale_x_continuous(limits = c(0,74000)) +
  labs(x = "AES 3km ('000s)", y = "Predicted Butterfly Richness") +
  NULL
rich_3km

ggsave(paste0(modpath,"Combined plots/Butterfly richness 3km individual schemes.png"), height = 800, width = 1000, units = "mm", scale = 0.15)

rich_1km + rich_3km + plot_layout(guides='collect') &
  theme(legend.position='bottom')
ggsave(paste0(modpath,"Combined plots/Butterfly richness combined 1km and 3km individual schemes.png"), height = 800, width = 1200, units = "mm", scale = 0.15)



###Abundance###



range(Abund_LS_mod$data$AES1KM)
range(Abund_WCBS_mod$data$AES1KM)
range(Abund_UKBMS_mod$data$AES1KM)




mean(c(Abund_LS_mod$data$AES3KM, Abund_WCBS_mod$data$AES3KM, Abund_UKBMS_mod$data$AES3KM))

mean(c(Abund_LS_mod$data$Climate_PC1, Abund_WCBS_mod$data$Climate_PC1, Abund_UKBMS_mod$data$Climate_PC1))
mean(c(Abund_LS_mod$data$Landscape_PC1, Abund_WCBS_mod$data$Landscape_PC1, Abund_UKBMS_mod$data$Landscape_PC1))
mean(c(Abund_LS_mod$data$Habitat_PC1, Abund_WCBS_mod$data$Habitat_PC1, Abund_UKBMS_mod$data$Habitat_PC1))

#scaled value for 4 rounds = 0.4

p1 <- ggpredict(Abund_LS_mod, "AES1KM[-0.60:7.03, by = 0.5]",
                condition = c("AES3KM" = "0.122",
                              "ROUND_NUMBER" = "0.4",
                              "SURVEY_YEAR" = "2018",
                              "Climate_PC1" = "-0.04",
                              "Landscape_PC1" = "-0.05",
                              "Habitat_PC1" = "0.26"))

p2 <- ggpredict(Abund_WCBS_mod, "AES1KM[-0.60:7.97, by = 0.5]",
                condition = c("AES3KM" = "0.122",
                              "N_VISITS_MAYTOAUGUST" = "0.4",
                              "YR" = "2018",
                              "Climate_PC1" = "-0.04",
                              "Landscape_PC1" = "-0.05",
                              "Habitat_PC1" = "0.26"))

p3 <- ggpredict(Abund_UKBMS_mod, "AES1KM[-0.60:8.21, by = 0.5]",
                condition = c("AES3KM" = "0.122",
                              "N_VISITS_MAYTOAUGUST" = "0.4",
                              "YR" = "2018",
					"TRANSECT_LENGTH_NEW" = "2",
                              "Climate_PC1" = "-0.04",
                              "Landscape_PC1" = "-0.05",
                              "Habitat_PC1" = "0.26"))

p1$group <- "LandSpAES"
p2$group <- "WCBS"
p3$group <- "UKBMS"


p4 <- do.call(rbind, list(p2, p1, p3)) #%>%
  mutate(x = (x*5000 + 3000)/1000)
abund_1km <- ggplot(p4, aes(x = x, y = predicted, colour = group, fill = group)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.3, colour = NA) +
  geom_line() +
  # geom_rug(data = p1_raw, aes(x = x, y = response, colour = group),
  #          sides = "b") +
  # geom_rug(data = p2_raw, aes(x = x, y = response, colour = group),
  #          sides = "t") +
  # geom_rug(data = p3_raw, aes(x = x, y = response, colour = group),
  #          sides = "t", outside = TRUE) +
  coord_cartesian(clip = "off") +
  scale_fill_manual(aesthetics = c("fill","colour"),
                    values = my_col,
                    name = "Survey") +
  scale_y_continuous(limits = c(0,1000), expand = c(0,0)) +
  # scale_x_continuous(limits = c(0,74000)) +
  labs(x = "AES 1km ('000s)", y = "Predicted Butterfly Abundance") +
  NULL
abund_1km

ggsave(paste0(modpath,"Combined plots/Butterfly abundance 1km individual schemes.png"), height = 800, width = 1000, units = "mm", scale = 0.15)

##3km plots


range(Abund_LS_mod$data$AES3KM)
range(Abund_WCBS_mod$data$AES3KM)
range(Abund_UKBMS_mod$data$AES3KM)

mean(c(Abund_LS_mod$data$AES1KM, Abund_WCBS_mod$data$AES1KM, Abund_UKBMS_mod$data$AES1KM))


p1 <- ggpredict(Abund_LS_mod, "AES3KM[-0.60:2.42, by = 0.5]",
                condition = c("AES1KM" = "0.362",
                              "ROUND_NUMBER" = "0.4",
                              "SURVEY_YEAR" = "2018",
                              "Climate_PC1" = "-0.04",
                              "Landscape_PC1" = "-0.05",
                              "Habitat_PC1" = "0.26"))

p2 <- ggpredict(Abund_WCBS_mod, "AES3KM[-0.60:6.66, by = 0.5]",
                condition = c("AES1KM" = "0.362",
                              "N_VISITS_MAYTOAUGUST" = "0.4",
                              "YR" = "2018",
                              "Climate_PC1" = "-0.04",
                              "Landscape_PC1" = "-0.05",
                              "Habitat_PC1" = "0.26"))

p3 <- ggpredict(Abund_UKBMS_mod, "AES3KM[-0.60:7.38, by = 0.5]",
                condition = c("AES1KM" = "0.362",
                              "N_VISITS_MAYTOAUGUST" = "0.4",
                              "YR" = "2018",
					"TRANSECT_LENGTH_NEW" = "2",
                              "Climate_PC1" = "-0.04",
                              "Landscape_PC1" = "-0.05",
                              "Habitat_PC1" = "0.26"))

p1$group <- "LandSpAES"
p2$group <- "WCBS"
p3$group <- "UKBMS"


p4 <- do.call(rbind, list(p2, p1, p3)) #%>%
  mutate(x = (x*5000 + 3000)/1000)
abund_3km <- ggplot(p4, aes(x = x, y = predicted, colour = group, fill = group)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.3, colour = NA) +
  geom_line() +
  # geom_rug(data = p1_raw, aes(x = x, y = response, colour = group),
  #          sides = "b") +
  # geom_rug(data = p2_raw, aes(x = x, y = response, colour = group),
  #          sides = "t") +
  # geom_rug(data = p3_raw, aes(x = x, y = response, colour = group),
  #          sides = "t", outside = TRUE) +
  coord_cartesian(clip = "off") +
  scale_fill_manual(aesthetics = c("fill","colour"),
                    values = my_col,
                    name = "Survey") +
  scale_y_continuous(limits = c(0,1000), expand = c(0,0)) +
  # scale_x_continuous(limits = c(0,74000)) +
  labs(x = "AES 3km ('000s)", y = "Predicted Butterfly Abundance") +
  NULL
abund_3km

ggsave(paste0(modpath,"Combined plots/Butterfly abundance 3km individual schemes.png"), height = 800, width = 1000, units = "mm", scale = 0.15)

abund_1km + abund_3km + plot_layout(guides='collect') &
  theme(legend.position='bottom')
ggsave(paste0(modpath,"Combined plots/Butterfly abundance combined 1km and 3km individual schemes.png"), height = 800, width = 1200, units = "mm", scale = 0.15)



###Diversity###



range(Div_LS_mod$data$AES1KM)
range(Div_WCBS_mod$data$AES1KM)
range(Div_UKBMS_mod$data$AES1KM)




mean(c(Div_LS_mod$data$AES3KM, Div_WCBS_mod$data$AES3KM, Div_UKBMS_mod$data$AES3KM))

mean(c(Div_LS_mod$data$Climate_PC1, Div_WCBS_mod$data$Climate_PC1, Div_UKBMS_mod$data$Climate_PC1))
mean(c(Div_LS_mod$data$Landscape_PC1, Div_WCBS_mod$data$Landscape_PC1, Div_UKBMS_mod$data$Landscape_PC1))
mean(c(Div_LS_mod$data$Habitat_PC1, Div_WCBS_mod$data$Habitat_PC1, Div_UKBMS_mod$data$Habitat_PC1))

#scaled value for 4 rounds = 0.4

p1 <- ggpredict(Div_LS_mod, "AES1KM[-0.60:7.03, by = 0.5]",
                condition = c("AES3KM" = "0.122",
                              "ROUND_NUMBER" = "0.4",
                              "SURVEY_YEAR" = "2018",
                              "Climate_PC1" = "-0.04",
                              "Landscape_PC1" = "-0.05",
                              "Habitat_PC1" = "0.26"))

p2 <- ggpredict(Div_WCBS_mod, "AES1KM[-0.60:7.97, by = 0.5]",
                condition = c("AES3KM" = "0.122",
                              "N_VISITS_MAYTOAUGUST" = "0.4",
                              "YR" = "2018",
                              "Climate_PC1" = "-0.04",
                              "Landscape_PC1" = "-0.05",
                              "Habitat_PC1" = "0.26"))

p3 <- ggpredict(Div_UKBMS_mod, "AES1KM[-0.60:8.21, by = 0.5]",
                condition = c("AES3KM" = "0.122",
                              "N_VISITS_MAYTOAUGUST" = "0.4",
                              "YR" = "2018",
					"TRANSECT_LENGTH_NEW" = "2",
                              "Climate_PC1" = "-0.04",
                              "Landscape_PC1" = "-0.05",
                              "Habitat_PC1" = "0.26"))

p1$group <- "LandSpAES"
p2$group <- "WCBS"
p3$group <- "UKBMS"


p4 <- do.call(rbind, list(p2, p1, p3)) #%>%
  mutate(x = (x*5000 + 3000)/1000)
div_1km <- ggplot(p4, aes(x = x, y = predicted, colour = group, fill = group)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.3, colour = NA) +
  geom_line() +
  # geom_rug(data = p1_raw, aes(x = x, y = response, colour = group),
  #          sides = "b") +
  # geom_rug(data = p2_raw, aes(x = x, y = response, colour = group),
  #          sides = "t") +
  # geom_rug(data = p3_raw, aes(x = x, y = response, colour = group),
  #          sides = "t", outside = TRUE) +
  coord_cartesian(clip = "off") +
  scale_fill_manual(aesthetics = c("fill","colour"),
                    values = my_col,
                    name = "Survey") +
  scale_y_continuous(limits = c(0,15), expand = c(0,0)) +
  # scale_x_continuous(limits = c(0,74000)) +
  labs(x = "AES 1km ('000s)", y = "Predicted Butterfly Diversity") +
  NULL
div_1km

ggsave(paste0(modpath,"Combined plots/Butterfly diversity 1km individual schemes.png"), height = 800, width = 1000, units = "mm", scale = 0.15)

##3km plots


range(Div_LS_mod$data$AES3KM)
range(Div_WCBS_mod$data$AES3KM)
range(Div_UKBMS_mod$data$AES3KM)

mean(c(Div_LS_mod$data$AES1KM, Div_WCBS_mod$data$AES1KM, Div_UKBMS_mod$data$AES1KM))


p1 <- ggpredict(Div_LS_mod, "AES3KM[-0.60:2.42, by = 0.5]",
                condition = c("AES1KM" = "0.362",
                              "ROUND_NUMBER" = "0.4",
                              "SURVEY_YEAR" = "2018",
                              "Climate_PC1" = "-0.04",
                              "Landscape_PC1" = "-0.05",
                              "Habitat_PC1" = "0.26"))

p2 <- ggpredict(Div_WCBS_mod, "AES3KM[-0.60:6.66, by = 0.5]",
                condition = c("AES1KM" = "0.362",
                              "N_VISITS_MAYTOAUGUST" = "0.4",
                              "YR" = "2018",
                              "Climate_PC1" = "-0.04",
                              "Landscape_PC1" = "-0.05",
                              "Habitat_PC1" = "0.26"))

p3 <- ggpredict(Div_UKBMS_mod, "AES3KM[-0.60:7.38, by = 0.5]",
                condition = c("AES1KM" = "0.362",
                              "N_VISITS_MAYTOAUGUST" = "0.4",
                              "YR" = "2018",
					"TRANSECT_LENGTH_NEW" = "2",
                              "Climate_PC1" = "-0.04",
                              "Landscape_PC1" = "-0.05",
                              "Habitat_PC1" = "0.26"))

p1$group <- "LandSpAES"
p2$group <- "WCBS"
p3$group <- "UKBMS"

p4 <- do.call(rbind, list(p2, p1, p3)) #%>%
  mutate(x = (x*5000 + 3000)/1000)
div_3km <- ggplot(p4, aes(x = x, y = predicted, colour = group, fill = group)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.3, colour = NA) +
  geom_line() +
  # geom_rug(data = p1_raw, aes(x = x, y = response, colour = group),
  #          sides = "b") +
  # geom_rug(data = p2_raw, aes(x = x, y = response, colour = group),
  #          sides = "t") +
  # geom_rug(data = p3_raw, aes(x = x, y = response, colour = group),
  #          sides = "t", outside = TRUE) +
  coord_cartesian(clip = "off") +
  scale_fill_manual(aesthetics = c("fill","colour"),
                    values = my_col,
                    name = "Survey") +
  scale_y_continuous(limits = c(0,15), expand = c(0,0)) +
  # scale_x_continuous(limits = c(0,74000)) +
  labs(x = "AES 3km ('000s)", y = "Predicted Butterfly Diversity") +
  NULL
div_3km

ggsave(paste0(modpath,"Combined plots/Butterfly diversity 3km individual schemes.png"), height = 800, width = 1000, units = "mm", scale = 0.15)

div_1km + div_3km + plot_layout(guides='collect') &
  theme(legend.position='bottom')
ggsave(paste0(modpath,"Combined plots/Butterfly diversity combined 1km and 3km individual schemes.png"), height = 800, width = 1200, units = "mm", scale = 0.15)


#####Add z test equivalent by calculating difference between posteriors and then assessing overlap with 1

