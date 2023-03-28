
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


Low_LS_mod <- readRDS(paste0(modpath, "LandSpAES_Lowmob_Abundance_brm.RDS"))
Med_LS_mod <- readRDS(paste0(modpath, "LandSpAES_Medmob_Abundance_brm.RDS"))
High_LS_mod <- readRDS(paste0(modpath, "LandSpAES_Highmob_Abundance_brm.RDS"))

Low_WCBS_mod <- readRDS(paste0(modpath, "WCBS_LM_Abundance_brm.RDS"))
Med_WCBS_mod <- readRDS(paste0(modpath, "WCBS_MM_Abundance_brm.RDS"))
High_WCBS_mod <- readRDS(paste0(modpath, "WCBS_HM_Abundance_brm.RDS"))

Low_UKBMS_mod <- readRDS(paste0(modpath, "UKBMS_LM_Abundance_brm.RDS"))
Med_UKBMS_mod <- readRDS(paste0(modpath, "UKBMS_MM_Abundance_brm.RDS"))
High_UKBMS_mod <- readRDS(paste0(modpath, "UKBMS_HM_Abundance_brm.RDS"))



#result tables AES effects only (for reference in BES talk)


m1 <- summary(Low_LS_mod)$fixed[c(2:3,11),]
m1$Model <- "Low_LS"
m2 <- summary(Low_UKBMS_mod)$fixed[c(2:3,12),]
m2$Model <- "Low_UKBMS"
m3 <- summary(Low_WCBS_mod)$fixed[c(2:3,12),]
m3$Model <- "Low_WCBS"
m4 <- summary(Med_LS_mod)$fixed[c(2:3,11),]
m4$Model <- "Med_LS"
m5 <- summary(Med_UKBMS_mod)$fixed[c(2:3,12),]
m5$Model <- "Med_UKBMS"
m6 <- summary(Med_WCBS_mod)$fixed[c(2:3,12),]
m6$Model <- "Med_WCBS"
m7 <- summary(High_LS_mod)$fixed[c(2:3,11),]
m7$Model <- "High_LS"
m8 <- summary(High_UKBMS_mod)$fixed[c(2:3,12),]
m8$Model <- "High_UKBMS"
m9 <- summary(High_WCBS_mod)$fixed[c(2:3,12),]
m9$Model <- "High_WCBS"

write.csv(rbind(m1,m2,m3,m4,m5,m6,m7,m8,m9), "Results_summary_mobility.csv")




range(Low_LS_mod$data$AES1KM)
range(Low_WCBS_mod$data$AES1KM)
range(Low_UKBMS_mod$data$AES1KM)

## add AES score histograms for BES


AES_data <- smartbind(Low_LS_mod$data, Low_WCBS_mod$data, Low_UKBMS_mod$data)
AES_data$SURVEY <- c(rep("LandSpAES", nrow(Low_LS_mod$data)),
                     rep("WCBS", nrow(Low_WCBS_mod$data)),
                     rep("UKBMS", nrow(Low_UKBMS_mod$data)))

ggplot(AES_data, aes(x = AES1KM, fill = SURVEY))+
  geom_density(alpha = .25)
ggplot(AES_data, aes(x = AES3KM, fill = SURVEY))+
  geom_density(alpha = .25)+
  xlab("Landscape AES score") +
  ylab("Density")+
  scale_fill_manual(values = c("#F8766D", "#619CFF", "#00BA38"))

###




mean(c(Low_LS_mod$data$AES3KM, Low_WCBS_mod$data$AES3KM, Low_UKBMS_mod$data$AES3KM))

mean(c(Low_LS_mod$data$Climate_PC1, Low_WCBS_mod$data$Climate_PC1, Low_UKBMS_mod$data$Climate_PC1))
mean(c(Low_LS_mod$data$Landscape_PC1, Low_WCBS_mod$data$Landscape_PC1, Low_UKBMS_mod$data$Landscape_PC1))
mean(c(Low_LS_mod$data$Habitat_PC1, Low_WCBS_mod$data$Habitat_PC1, Low_UKBMS_mod$data$Habitat_PC1))

#scaled value for 4 rounds = 0.4

p1 <- ggpredict(Low_LS_mod, "AES1KM[-0.60:7.03, by = 0.5]",
                condition = c("AES3KM" = "0.122",
                              "ROUND_NUMBER" = "0.4",
                              "SURVEY_YEAR" = "2018",
                              "Climate_PC1" = "-0.04",
                              "Landscape_PC1" = "-0.05",
                              "Habitat_PC1" = "0.26"))

p2 <- ggpredict(Low_WCBS_mod, "AES1KM[-0.60:7.97, by = 0.5]",
                condition = c("AES3KM" = "0.122",
                              "N_VISITS_MAYTOAUGUST" = "0.4",
                              "YR" = "2018",
                              "Climate_PC1" = "-0.04",
                              "Landscape_PC1" = "-0.05",
                              "Habitat_PC1" = "0.26"))

p3 <- ggpredict(Low_UKBMS_mod, "AES1KM[-0.60:8.21, by = 0.5]",
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
low_1km <- ggplot(p4, aes(x = x, y = predicted, colour = group, fill = group)) + 
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
  scale_y_continuous(limits = c(0,500), expand = c(0,0)) +
  # scale_x_continuous(limits = c(0,74000)) +
  labs(x = "AES 1km ('000s)", y = "Predicted Butterfly Abundance - low mobility") +
  NULL
low_1km

ggsave(paste0(modpath,"Combined plots/Butterfly abundance 1km low mobility.png"), height = 800, width = 1000, units = "mm", scale = 0.15)

##3km plots


range(Low_LS_mod$data$AES3KM)
range(Low_WCBS_mod$data$AES3KM)
range(Low_UKBMS_mod$data$AES3KM)

mean(c(Low_LS_mod$data$AES1KM, Low_WCBS_mod$data$AES1KM, Low_UKBMS_mod$data$AES1KM))


p1 <- ggpredict(Low_LS_mod, "AES3KM[-0.60:2.42, by = 0.5]",
                condition = c("AES1KM" = "0.362",
                              "ROUND_NUMBER" = "0.4",
                              "SURVEY_YEAR" = "2018",
                              "Climate_PC1" = "-0.04",
                              "Landscape_PC1" = "-0.05",
                              "Habitat_PC1" = "0.26"))

p2 <- ggpredict(Low_WCBS_mod, "AES3KM[-0.60:6.66, by = 0.5]",
                condition = c("AES1KM" = "0.362",
                              "N_VISITS_MAYTOAUGUST" = "0.4",
                              "YR" = "2018",
                              "Climate_PC1" = "-0.04",
                              "Landscape_PC1" = "-0.05",
                              "Habitat_PC1" = "0.26"))

p3 <- ggpredict(Low_UKBMS_mod, "AES3KM[-0.60:7.38, by = 0.5]",
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
low_3km <- ggplot(p4, aes(x = x, y = predicted, colour = group, fill = group)) + 
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
  scale_y_continuous(limits = c(0,500), expand = c(0,0)) +
  # scale_x_continuous(limits = c(0,74000)) +
  labs(x = "AES 3km ('000s)", y = "Predicted Butterfly Abundance - low mobility") +
  NULL
low_3km

ggsave(paste0(modpath,"Combined plots/Butterfly abundance 3km low mobility.png"), height = 800, width = 1000, units = "mm", scale = 0.15)

low_1km + low_3km + plot_layout(guides='collect') &
  theme(legend.position='bottom')
ggsave(paste0(modpath,"Combined plots/Butterfly abundance combined 1km and 3km low mobility.png"), height = 800, width = 1200, units = "mm", scale = 0.15)



###Abundance###



range(Med_LS_mod$data$AES1KM)
range(Med_WCBS_mod$data$AES1KM)
range(Med_UKBMS_mod$data$AES1KM)




mean(c(Med_LS_mod$data$AES3KM, Med_WCBS_mod$data$AES3KM, Med_UKBMS_mod$data$AES3KM))

mean(c(Med_LS_mod$data$Climate_PC1, Med_WCBS_mod$data$Climate_PC1, Med_UKBMS_mod$data$Climate_PC1))
mean(c(Med_LS_mod$data$Landscape_PC1, Med_WCBS_mod$data$Landscape_PC1, Med_UKBMS_mod$data$Landscape_PC1))
mean(c(Med_LS_mod$data$Habitat_PC1, Med_WCBS_mod$data$Habitat_PC1, Med_UKBMS_mod$data$Habitat_PC1))

#scaled value for 4 rounds = 0.4

p1 <- ggpredict(Med_LS_mod, "AES1KM[-0.60:7.03, by = 0.5]",
                condition = c("AES3KM" = "0.122",
                              "ROUND_NUMBER" = "0.4",
                              "SURVEY_YEAR" = "2018",
                              "Climate_PC1" = "-0.04",
                              "Landscape_PC1" = "-0.05",
                              "Habitat_PC1" = "0.26"))

p2 <- ggpredict(Med_WCBS_mod, "AES1KM[-0.60:7.97, by = 0.5]",
                condition = c("AES3KM" = "0.122",
                              "N_VISITS_MAYTOAUGUST" = "0.4",
                              "YR" = "2018",
                              "Climate_PC1" = "-0.04",
                              "Landscape_PC1" = "-0.05",
                              "Habitat_PC1" = "0.26"))

p3 <- ggpredict(Med_UKBMS_mod, "AES1KM[-0.60:8.21, by = 0.5]",
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
med_1km <- ggplot(p4, aes(x = x, y = predicted, colour = group, fill = group)) + 
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
  labs(x = "AES 1km ('000s)", y = "Predicted Butterfly Abundance - medium mobility") +
  NULL
med_1km

ggsave(paste0(modpath,"Combined plots/Butterfly abundance medium 1km mobility.png"), height = 800, width = 1000, units = "mm", scale = 0.15)

##3km plots


range(Med_LS_mod$data$AES3KM)
range(Med_WCBS_mod$data$AES3KM)
range(Med_UKBMS_mod$data$AES3KM)

mean(c(Med_LS_mod$data$AES1KM, Med_WCBS_mod$data$AES1KM, Med_UKBMS_mod$data$AES1KM))


p1 <- ggpredict(Med_LS_mod, "AES3KM[-0.60:2.42, by = 0.5]",
                condition = c("AES1KM" = "0.362",
                              "ROUND_NUMBER" = "0.4",
                              "SURVEY_YEAR" = "2018",
                              "Climate_PC1" = "-0.04",
                              "Landscape_PC1" = "-0.05",
                              "Habitat_PC1" = "0.26"))

p2 <- ggpredict(Med_WCBS_mod, "AES3KM[-0.60:6.66, by = 0.5]",
                condition = c("AES1KM" = "0.362",
                              "N_VISITS_MAYTOAUGUST" = "0.4",
                              "YR" = "2018",
                              "Climate_PC1" = "-0.04",
                              "Landscape_PC1" = "-0.05",
                              "Habitat_PC1" = "0.26"))

p3 <- ggpredict(Med_UKBMS_mod, "AES3KM[-0.60:7.38, by = 0.5]",
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
med_3km <- ggplot(p4, aes(x = x, y = predicted, colour = group, fill = group)) + 
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
  labs(x = "AES 3km ('000s)", y = "Predicted Butterfly Abundance - medium mobility") +
  NULL
med_3km

ggsave(paste0(modpath,"Combined plots/Butterfly abundance 3km medium mobility.png"), height = 800, width = 1000, units = "mm", scale = 0.15)

med_1km + med_3km + plot_layout(guides='collect') &
  theme(legend.position='bottom')
ggsave(paste0(modpath,"Combined plots/Butterfly abundance combined 1km and 3km medium mobility.png"), height = 800, width = 1200, units = "mm", scale = 0.15)



###Diversity###



range(High_LS_mod$data$AES1KM)
range(High_WCBS_mod$data$AES1KM)
range(High_UKBMS_mod$data$AES1KM)




mean(c(High_LS_mod$data$AES3KM, High_WCBS_mod$data$AES3KM, High_UKBMS_mod$data$AES3KM))

mean(c(High_LS_mod$data$Climate_PC1, High_WCBS_mod$data$Climate_PC1, High_UKBMS_mod$data$Climate_PC1))
mean(c(High_LS_mod$data$Landscape_PC1, High_WCBS_mod$data$Landscape_PC1, High_UKBMS_mod$data$Landscape_PC1))
mean(c(High_LS_mod$data$Habitat_PC1, High_WCBS_mod$data$Habitat_PC1, High_UKBMS_mod$data$Habitat_PC1))

#scaled value for 4 rounds = 0.4

p1 <- ggpredict(High_LS_mod, "AES1KM[-0.60:7.03, by = 0.5]",
                condition = c("AES3KM" = "0.122",
                              "ROUND_NUMBER" = "0.4",
                              "SURVEY_YEAR" = "2018",
                              "Climate_PC1" = "-0.04",
                              "Landscape_PC1" = "-0.05",
                              "Habitat_PC1" = "0.26"))

p2 <- ggpredict(High_WCBS_mod, "AES1KM[-0.60:7.97, by = 0.5]",
                condition = c("AES3KM" = "0.122",
                              "N_VISITS_MAYTOAUGUST" = "0.4",
                              "YR" = "2018",
                              "Climate_PC1" = "-0.04",
                              "Landscape_PC1" = "-0.05",
                              "Habitat_PC1" = "0.26"))

p3 <- ggpredict(High_UKBMS_mod, "AES1KM[-0.60:8.21, by = 0.5]",
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
high_1km <- ggplot(p4, aes(x = x, y = predicted, colour = group, fill = group)) + 
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
  scale_y_continuous(limits = c(0,200), expand = c(0,0)) +
  # scale_x_continuous(limits = c(0,74000)) +
  labs(x = "AES 1km ('000s)", y = "Predicted Butterfly Abundance - high mobility") +
  NULL
high_1km

ggsave(paste0(modpath,"Combined plots/Butterfly abundance 1km high mobility.png"), height = 800, width = 1000, units = "mm", scale = 0.15)

##3km plots


range(High_LS_mod$data$AES3KM)
range(High_WCBS_mod$data$AES3KM)
range(High_UKBMS_mod$data$AES3KM)

mean(c(High_LS_mod$data$AES1KM, High_WCBS_mod$data$AES1KM, High_UKBMS_mod$data$AES1KM))


p1 <- ggpredict(High_LS_mod, "AES3KM[-0.60:2.42, by = 0.5]",
                condition = c("AES1KM" = "0.362",
                              "ROUND_NUMBER" = "0.4",
                              "SURVEY_YEAR" = "2018",
                              "Climate_PC1" = "-0.04",
                              "Landscape_PC1" = "-0.05",
                              "Habitat_PC1" = "0.26"))

p2 <- ggpredict(High_WCBS_mod, "AES3KM[-0.60:6.66, by = 0.5]",
                condition = c("AES1KM" = "0.362",
                              "N_VISITS_MAYTOAUGUST" = "0.4",
                              "YR" = "2018",
                              "Climate_PC1" = "-0.04",
                              "Landscape_PC1" = "-0.05",
                              "Habitat_PC1" = "0.26"))

p3 <- ggpredict(High_UKBMS_mod, "AES3KM[-0.60:7.38, by = 0.5]",
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
high_3km <- ggplot(p4, aes(x = x, y = predicted, colour = group, fill = group)) + 
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
  scale_y_continuous(limits = c(0,200), expand = c(0,0)) +
  # scale_x_continuous(limits = c(0,74000)) +
  labs(x = "AES 3km ('000s)", y = "Predicted Butterfly Abundance - high mobility") +
  NULL
high_3km

ggsave(paste0(modpath,"Combined plots/Butterfly abundance 3km high mobility.png"), height = 800, width = 1000, units = "mm", scale = 0.15)

high_1km + high_3km + plot_layout(guides='collect') &
  theme(legend.position='bottom')
ggsave(paste0(modpath,"Combined plots/Butterfly abundance combined 1km and 3km high mobility.png"), height = 800, width = 1200, units = "mm", scale = 0.15)


#####Add z test equivalent by calculating difference between posteriors and then assessing overlap with 1



#try idea of z test to test similarity
ztest <- function(mod1, mod2, term){
  p <- vector(); z <- vector()
  for (i in term){
    b1 <- summary(mod1)$fixed[row.names(summary(mod1)$fixed) == i,1]
    b2 <- summary(mod2)$fixed[row.names(summary(mod2)$fixed) == i,1]
    se1 <- sqrt(diag(vcov(mod1)))[names(sqrt(diag(vcov(mod1)))) == i]
    se2 <- sqrt(diag(vcov(mod2)))[names(sqrt(diag(vcov(mod2)))) == i]
    z[which(term == i)] <- (b1-b2)/(sqrt((se1^2) + (se2^2)))
    p[which(term == i)] <- 2*pnorm(abs(z[which(term == i)]), lower.tail = FALSE)
  }
  z.df <- data.frame(term = term, z = z, p = p)
  return(z.df)
}

#richness
ztest(Low_LS_mod, Low_WCBS_mod, c("AES1KM", "AES3KM", "AES1KM:AES3KM"))#OK
ztest(Low_LS_mod, Low_UKBMS_mod, c("AES1KM", "AES3KM", "AES1KM:AES3KM"))#OK
ztest(Low_UKBMS_mod, Low_WCBS_mod, c("AES1KM", "AES3KM", "AES1KM:AES3KM"))#OK

#abundance
ztest(Med_LS_mod, Med_WCBS_mod, c("AES1KM", "AES3KM", "AES1KM:AES3KM"))#OK
ztest(Med_LS_mod, Med_UKBMS_mod, c("AES1KM", "AES3KM", "AES1KM:AES3KM"))#OK
ztest(Med_UKBMS_mod, Med_WCBS_mod, c("AES1KM", "AES3KM", "AES1KM:AES3KM"))#OK

#diversity
ztest(High_LS_mod, High_WCBS_mod, c("AES1KM", "AES3KM", "AES1KM:AES3KM"))#OK
ztest(High_LS_mod, High_UKBMS_mod, c("AES1KM", "AES3KM", "AES1KM:AES3KM"))#OK
ztest(High_UKBMS_mod, High_WCBS_mod, c("AES1KM", "AES3KM", "AES1KM:AES3KM"))#OK
#richness and abundance definitely OK, diversity borderline
#may need to remove UKBMS from abundance and diversity?



