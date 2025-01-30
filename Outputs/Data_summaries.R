#Data summaries
library(dplyr)

buttabund <- read.csv("LandSpAES butterfly abundance.csv")
buttdiv <- read.csv("LandSpAES butterfly diversity.csv")
buttrich <- read.csv("LandSpAES butterfly richness.csv")

UKBMS_all_data <- read.csv("UKBMS butterfly data.csv")

WCBS_all_data <- read.csv("WCBS butterfly data.csv")

LandSpAES_summ1 <- buttabund %>% 
  summarise(Mean = mean(BUTTERFLY_COUNT),
            Median = median(BUTTERFLY_COUNT),
            Standard_deviation = sd(BUTTERFLY_COUNT),
            Min = min(BUTTERFLY_COUNT),
            Max = max(BUTTERFLY_COUNT)) %>%
  mutate(Response = "Abundance", Dataset = "LandSpAES")
  
LandSpAES_summ2 <- buttdiv %>% 
  summarise(Mean = mean(SHANNON_DIV),
            Median = median(SHANNON_DIV),
            Standard_deviation = sd(SHANNON_DIV),
            Min = min(SHANNON_DIV),
            Max = max(SHANNON_DIV)) %>%
  mutate(Response = "Diversity", Dataset = "LandSpAES")

LandSpAES_summ3 <- buttrich %>% 
  summarise(Mean = mean(RICHNESS_ID),
            Median = median(RICHNESS_ID),
            Standard_deviation = sd(RICHNESS_ID),
            Min = min(RICHNESS_ID),
            Max = max(RICHNESS_ID)) %>%
  mutate(Response = "Richness", Dataset = "LandSpAES")


LS <- do.call("rbind", list(LandSpAES_summ1, LandSpAES_summ2, LandSpAES_summ3))

UKBMS_summ1 <- UKBMS_all_data %>% 
  summarise(Mean = mean(Abundance),
            Median = median(Abundance),
            Standard_deviation = sd(Abundance),
            Min = min(Abundance),
            Max = max(Abundance)) %>%
  mutate(Response = "Abundance", Dataset = "UKBMS")

UKBMS_summ2 <- UKBMS_all_data %>% 
  summarise(Mean = mean(Shannon_diversity),
            Median = median(Shannon_diversity),
            Standard_deviation = sd(Shannon_diversity),
            Min = min(Shannon_diversity),
            Max = max(Shannon_diversity)) %>%
  mutate(Response = "Diversity", Dataset = "UKBMS")

UKBMS_summ3 <- UKBMS_all_data %>% 
  summarise(Mean = mean(Richness),
            Median = median(Richness),
            Standard_deviation = sd(Richness),
            Min = min(Richness),
            Max = max(Richness)) %>%
  mutate(Response = "Richness", Dataset = "UKBMS")

UKBMS <- do.call("rbind", list(UKBMS_summ1, UKBMS_summ2, UKBMS_summ3))

WCBS_summ1 <- WCBS_all_data %>% 
  summarise(Mean = mean(Abundance),
            Median = median(Abundance),
            Standard_deviation = sd(Abundance),
            Min = min(Abundance),
            Max = max(Abundance)) %>%
  mutate(Response = "Abundance", Dataset = "WCBS")

WCBS_summ2 <- WCBS_all_data %>% 
  summarise(Mean = mean(Shannon_diversity),
            Median = median(Shannon_diversity),
            Standard_deviation = sd(Shannon_diversity),
            Min = min(Shannon_diversity),
            Max = max(Shannon_diversity)) %>%
  mutate(Response = "Diversity", Dataset = "WCBS")

WCBS_summ3 <- WCBS_all_data %>% 
  summarise(Mean = mean(Richness),
            Median = median(Richness),
            Standard_deviation = sd(Richness),
            Min = min(Richness),
            Max = max(Richness)) %>%
  mutate(Response = "Richness", Dataset = "WCBS")

WCBS <- do.call("rbind", list(WCBS_summ1, WCBS_summ2, WCBS_summ3))

all_summaries <- do.call("rbind", list(LS, UKBMS, WCBS))

all_summaries <- all_summaries[,c(6,7,1:5)]

all_summaries <- all_summaries %>%
  arrange(Response, Dataset)

all_summaries[,3:7] <- round(all_summaries[,3:7],2)


write.csv(all_summaries, "Data summaries for SI.csv")
