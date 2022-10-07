#' WCBS - Analyse each response to AES scores plus environmental gradients
#' 
library(MASS)
library(DHARMa)
library(car)
library(dplyr)
library(tidyr)
library(ggplot2)
theme_set(theme_classic())

#' Get data ####
#' Environmental gradients represented as PCA scores, sourced from file on P
#' drive.

fpath <- config_path

load(paste0(fpath,"Workfiles/Covariate selection/PCA scores per 1km square v2.Rdata"))

#' AES scores
AES <- read.csv(paste0(fpath, "Data/AES uptake/Outputs/All_1kmCells_Scored.csv"))

#' WCBS response data from file
source("Calculate_Responses_WCBS.R")

#' Manipulate data ####
#' Get PCA scores for WCBS
wcbs_pca_scores <- all_scores[all_scores[,1] %in% wcbs_responses$GRIDREF_1KM,]
colnames(wcbs_pca_scores)[1] <- "GRIDREF_1KM"
wcbs_pca_scores <- as.data.frame(wcbs_pca_scores)

#' Get AES scores for WCBS
wcbs_aes <- AES %>%
  filter(CELLCODE %in% wcbs_responses$GRIDREF_1KM) %>%
  select(CELLCODE, starts_with("Sc")) %>%
  tidyr::pivot_longer(starts_with("Sc"),
                      names_to = c("YEAR","scale"),
                      names_sep = "\\.") %>%
  mutate(YEAR = recode(YEAR,
                       "Sc17" = 2017,
                       "Sc18" = 2018,
                       "Sc19" = 2019,
                       "Sc20" = 2020),
         scale = recode(scale,
                        "1km" = "AES1KM",
                        "3km" = "AES3KM")) %>%
  tidyr::pivot_wider(names_from = scale,
                     values_from = value) %>%
  filter(YEAR != 2020)

#' Combine all data
WCBS_all_data <- inner_join(wcbs_responses, wcbs_pca_scores) %>%
  inner_join(wcbs_aes, by = c("GRIDREF_1KM" = "CELLCODE","YEAR")) %>%
  ungroup() %>% 
  mutate(across(starts_with("Comp"),as.numeric)) %>%
  filter(AES1KM < 50000) %>%
  mutate(YEAR = as.factor(YEAR),
         AES1KM = as.numeric(scale(AES1KM, center = 3340.947, scale = 6966.606)),
         AES3KM = as.numeric(scale(AES3KM, center = 3245.172, scale = 3745.525)))
summary(WCBS_all_data)
# psych::multi.hist(select_if(WCBS_all_data, is.numeric))
# par(mfrow=c(1,1))

#' Richness models ####
#'
#' First we fit the models using a poisson response for richness, and then check
#' if the model assumptions (mean = variance) are satisfied using the DHARMa
#' package. 
#' Poisson models: 
Rich_wcbs_mod <- glm(Richness ~ AES1KM*AES3KM + N_VISITS_MAYTOAUGUST + YEAR +
                       Comp.1 + Comp.2, 
                     data = WCBS_all_data, family = "poisson")

# test_wcbs_mod <- simulateResiduals(Rich_wcbs_mod, plot = TRUE)
# testDispersion(test_wcbs_mod)
# dispersion issues significant

#' Now we fit negative binomial models with the MASS package and see if these
#' fit the data. 
#' Negative binomial models:
Rich_wcbs_mod <- glm.nb(Richness ~ AES1KM*AES3KM + N_VISITS_MAYTOAUGUST + YEAR +
                          Comp.1 + Comp.2,
                        data = WCBS_all_data)
#' Check model assumptions - dispersion first, then look at variance inflation
#' to see if there is a collinearity problem.
test_wcbs_mod <- simulateResiduals(Rich_wcbs_mod, plot = TRUE)
# testOutliers(test_wcbs_mod, type = "bootstrap")

vif(Rich_wcbs_mod)
# no issues

summary(Rich_wcbs_mod)


#' Now we check to see if adding any of the remaining PCA axes improves model
#' fit (as measured by AIC) substantially. Looking for a reduction in AIC of at
#' least 4 here.

# All PCA axes
addterm(Rich_wcbs_mod, 
        scope = Richness ~ AES1KM*AES3KM + N_VISITS_MAYTOAUGUST + YEAR +
          Comp.1 + Comp.2 + Comp.3 + Comp.4 + Comp.5 + Comp.6 +
          Comp.7 + Comp.8 + Comp.9 + Comp.10 + Comp.11 + Comp.12 +
          Comp.13 + Comp.14 + Comp.15 + Comp.16 + Comp.17 + Comp.18 + 
          Comp.19 + Comp.20 + Comp.21 + Comp.22 + Comp.23 + Comp.24 +
          Comp.25 + Comp.26 + Comp.27 + Comp.28, sort = TRUE)

#' Butterfly richness model
Rich_wcbs_mod2 <- glm.nb(Richness ~ AES1KM*AES3KM + N_VISITS_MAYTOAUGUST + YEAR +
                           Comp.1 + Comp.2 + Comp.21, 
                         data = WCBS_all_data)
# test_wcbs_mod <- simulateResiduals(Rich_wcbs_mod2, plot = TRUE, n = 1000)
vif(Rich_wcbs_mod2)

summary(Rich_wcbs_mod2)

addterm(Rich_wcbs_mod2,
        scope = Richness ~ AES1KM*AES3KM + N_VISITS_MAYTOAUGUST + YEAR +
          Comp.1 + Comp.2 + Comp.3 + Comp.4 + Comp.5 + Comp.6 +
          Comp.7 + Comp.8 + Comp.9 + Comp.10 + Comp.11 + Comp.12 +
          Comp.13 + Comp.14 + Comp.15 + Comp.16 + Comp.17 + Comp.18 + 
          Comp.19 + Comp.20 + Comp.21 + Comp.22 + Comp.23 + Comp.24 +
          Comp.25 + Comp.26 + Comp.27 + Comp.28, 
        sort = TRUE)


Rich_wcbs_mod3 <- glm.nb(Richness ~ AES1KM*AES3KM + N_VISITS_MAYTOAUGUST + YEAR +
                           Comp.1 + Comp.2 + Comp.21 + Comp.6, 
                         data = WCBS_all_data)
# test_wcbs_mod <- simulateResiduals(Rich_wcbs_mod3, plot = TRUE)
vif(Rich_wcbs_mod3)

summary(Rich_wcbs_mod3)

addterm(Rich_wcbs_mod3,
        scope = Richness ~ AES1KM*AES3KM + N_VISITS_MAYTOAUGUST + YEAR +
          Comp.1 + Comp.2 + Comp.3 + Comp.4 + Comp.5 + Comp.6 +
          Comp.7 + Comp.8 + Comp.9 + Comp.10 + Comp.11 + Comp.12 +
          Comp.13 + Comp.14 + Comp.15 + Comp.16 + Comp.17 + Comp.18 + 
          Comp.19 + Comp.20 + Comp.21 + Comp.22 + Comp.23 + Comp.24 +
          Comp.25 + Comp.26 + Comp.27 + Comp.28, 
        sort = TRUE)


Rich_wcbs_mod4 <- glm.nb(Richness ~ AES1KM*AES3KM + N_VISITS_MAYTOAUGUST + YEAR +
                           Comp.1 + Comp.2 + Comp.6 + Comp.21 + Comp.23, 
                         data = WCBS_all_data)
# test_wcbs_mod <- simulateResiduals(Rich_wcbs_mod4, plot = TRUE)
vif(Rich_wcbs_mod4)

summary(Rich_wcbs_mod4)

addterm(Rich_wcbs_mod4,
        scope = Richness ~ AES1KM*AES3KM + N_VISITS_MAYTOAUGUST + YEAR +
          Comp.1 + Comp.2 + Comp.3 + Comp.4 + Comp.5 + Comp.6 +
          Comp.7 + Comp.8 + Comp.9 + Comp.10 + Comp.11 + Comp.12 +
          Comp.13 + Comp.14 + Comp.15 + Comp.16 + Comp.17 + Comp.18 + 
          Comp.19 + Comp.20 + Comp.21 + Comp.22 + Comp.23 + Comp.24 +
          Comp.25 + Comp.26 + Comp.27 + Comp.28, 
        sort = TRUE)

Rich_wcbs_mod5 <- glm.nb(Richness ~ AES1KM*AES3KM + N_VISITS_MAYTOAUGUST + YEAR +
                           Comp.1 + Comp.2 + Comp.6 + Comp.21 + Comp.23 + Comp.15, 
                         data = WCBS_all_data)
# test_wcbs_mod <- simulateResiduals(Rich_wcbs_mod5, plot = TRUE)
vif(Rich_wcbs_mod5)

summary(Rich_wcbs_mod5)

addterm(Rich_wcbs_mod5,
        scope = Richness ~ AES1KM*AES3KM + N_VISITS_MAYTOAUGUST + YEAR +
          Comp.1 + Comp.2 + Comp.3 + Comp.4 + Comp.5 + Comp.6 +
          Comp.7 + Comp.8 + Comp.9 + Comp.10 + Comp.11 + Comp.12 +
          Comp.13 + Comp.14 + Comp.15 + Comp.16 + Comp.17 + Comp.18 + 
          Comp.19 + Comp.20 + Comp.21 + Comp.22 + Comp.23 + Comp.24 +
          Comp.25 + Comp.26 + Comp.27 + Comp.28, 
        sort = TRUE)

Rich_wcbs_model_results <- do.call(rbind, list(
  data.frame(Variable = names(coef(Rich_wcbs_mod)),
             Coef = coef(Rich_wcbs_mod),
             SE = sqrt(diag(vcov(Rich_wcbs_mod))),
             pval = summary(Rich_wcbs_mod)$coefficients[,4],
             Model = "Rich_wcbs_mod",
             AIC = AIC(Rich_wcbs_mod)),
  data.frame(Variable = names(coef(Rich_wcbs_mod2)),
             Coef = coef(Rich_wcbs_mod2),
             SE = sqrt(diag(vcov(Rich_wcbs_mod2))),
             pval = summary(Rich_wcbs_mod2)$coefficients[,4],
             Model = "Rich_wcbs_mod2",
             AIC = AIC(Rich_wcbs_mod2)),
  data.frame(Variable = names(coef(Rich_wcbs_mod3)),
             Coef = coef(Rich_wcbs_mod3),
             SE = sqrt(diag(vcov(Rich_wcbs_mod3))),
             pval = summary(Rich_wcbs_mod3)$coefficients[,4],
             Model = "Rich_wcbs_mod3",
             AIC = AIC(Rich_wcbs_mod3)),
  data.frame(Variable = names(coef(Rich_wcbs_mod4)),
             Coef = coef(Rich_wcbs_mod4),
             SE = sqrt(diag(vcov(Rich_wcbs_mod4))),
             pval = summary(Rich_wcbs_mod4)$coefficients[,4],
             Model = "Rich_wcbs_mod4",
             AIC = AIC(Rich_wcbs_mod4)),
  data.frame(Variable = names(coef(Rich_wcbs_mod5)),
             Coef = coef(Rich_wcbs_mod5),
             SE = sqrt(diag(vcov(Rich_wcbs_mod5))),
             pval = summary(Rich_wcbs_mod5)$coefficients[,4],
             Model = "Rich_wcbs_mod5",
             AIC = AIC(Rich_wcbs_mod5))
)) %>%
  pivot_longer(Coef:pval) %>%
  pivot_wider(names_from = c(Variable,name),
              names_sep = "_",
              values_from = value,
              values_fill = NA) %>%
  mutate(Family = "Negative binomial")




#' ### Abundance models
# poisson models 
Abund_wcbs_mod <- glm(Abundance ~ AES1KM*AES3KM + N_VISITS_MAYTOAUGUST + YEAR +
                        Comp.1 + Comp.2, 
                      data = WCBS_all_data, family = "poisson")

summary(Abund_wcbs_mod)

# test_wcbs_mod <- simulateResiduals(Abund_wcbs_mod, plot = TRUE)
# There are many issues with the poisson models here

# negative binomial models
Abund_wcbs_mod <- glm.nb(Abundance ~ AES1KM*AES3KM +N_VISITS_MAYTOAUGUST + YEAR +
                           Comp.1 + Comp.2, 
                         data = WCBS_all_data)

summary(Abund_wcbs_mod)

#' Check model assumptions
# test_wcbs_mod <- simulateResiduals(Abund_wcbs_mod, plot = TRUE)
# no issues
vif(Abund_wcbs_mod)

#' Check to see if any PCA axes scores should be added
addterm(Abund_wcbs_mod, 
        scope = Abundance ~ AES1KM*AES3KM + N_VISITS_MAYTOAUGUST + YEAR +
          Comp.1 + Comp.2 + Comp.3 + Comp.4 + Comp.5 + Comp.6 +
          Comp.7 + Comp.8 + Comp.9 + Comp.10 + Comp.11 + Comp.12 +
          Comp.13 + Comp.14 + Comp.15 + Comp.16 + Comp.17 + Comp.18 + 
          Comp.19 + Comp.20 + Comp.21 + Comp.22 + Comp.23 + Comp.24 +
          Comp.25 + Comp.26 + Comp.27 + Comp.28, sort = TRUE)

#' abundance model with more PCA axes
Abund_wcbs_mod2 <- glm.nb(Abundance ~ AES1KM*AES3KM + N_VISITS_MAYTOAUGUST + YEAR +
                            Comp.1 + Comp.2 + Comp.6, 
                          data = WCBS_all_data)
# test_wcbs_mod <- simulateResiduals(Abund_wcbs_mod2, plot = TRUE)
vif(Abund_wcbs_mod2)
summary(Abund_wcbs_mod2)

addterm(Abund_wcbs_mod2,
        scope = Abundance ~ AES1KM*AES3KM + N_VISITS_MAYTOAUGUST + YEAR +
          Comp.1 + Comp.2 + Comp.3 + Comp.4 + Comp.5 + Comp.6 +
          Comp.7 + Comp.8 + Comp.9 + Comp.10 + Comp.11 + Comp.12 +
          Comp.13 + Comp.14 + Comp.15 + Comp.16 + Comp.17 + Comp.18 + 
          Comp.19 + Comp.20 + Comp.21 + Comp.22 + Comp.23 + Comp.24 +
          Comp.25 + Comp.26 + Comp.27 + Comp.28, 
        sort = TRUE)

Abund_wcbs_mod3 <- glm.nb(Abundance ~ AES1KM*AES3KM + N_VISITS_MAYTOAUGUST + YEAR +
                            Comp.1 + Comp.2 + Comp.6 + Comp.19, 
                          data = WCBS_all_data)
# test_wcbs_mod <- simulateResiduals(Abund_wcbs_mod3, plot = TRUE)
vif(Abund_wcbs_mod3)
summary(Abund_wcbs_mod3)

addterm(Abund_wcbs_mod3,
        scope = Abundance ~ AES1KM*AES3KM + N_VISITS_MAYTOAUGUST + YEAR +
          Comp.1 + Comp.2 + Comp.3 + Comp.4 + Comp.5 + Comp.6 +
          Comp.7 + Comp.8 + Comp.9 + Comp.10 + Comp.11 + Comp.12 +
          Comp.13 + Comp.14 + Comp.15 + Comp.16 + Comp.17 + Comp.18 + 
          Comp.19 + Comp.20 + Comp.21 + Comp.22 + Comp.23 + Comp.24 +
          Comp.25 + Comp.26 + Comp.27 + Comp.28, 
        sort = TRUE)

Abund_wcbs_mod4 <- glm.nb(Abundance ~ AES1KM*AES3KM + N_VISITS_MAYTOAUGUST + YEAR +
                            Comp.1 + Comp.2 + Comp.19 + Comp.6 + Comp.15, 
                          data = WCBS_all_data)
# test_wcbs_mod <- simulateResiduals(Abund_wcbs_mod4, plot = TRUE)
vif(Abund_wcbs_mod4)
summary(Abund_wcbs_mod4)

addterm(Abund_wcbs_mod4,
        scope = Abundance ~ AES1KM*AES3KM + N_VISITS_MAYTOAUGUST + YEAR +
          Comp.1 + Comp.2 + Comp.3 + Comp.4 + Comp.5 + Comp.6 +
          Comp.7 + Comp.8 + Comp.9 + Comp.10 + Comp.11 + Comp.12 +
          Comp.13 + Comp.14 + Comp.15 + Comp.16 + Comp.17 + Comp.18 + 
          Comp.19 + Comp.20 + Comp.21 + Comp.22 + Comp.23 + Comp.24 +
          Comp.25 + Comp.26 + Comp.27 + Comp.28, 
        sort = TRUE)

Abund_wcbs_mod5 <- glm.nb(Abundance ~ AES1KM*AES3KM + N_VISITS_MAYTOAUGUST + YEAR +
                            Comp.1 + Comp.2 + Comp.19 + Comp.6 + Comp.15 + Comp.20, 
                          data = WCBS_all_data)
# test_wcbs_mod <- simulateResiduals(Abund_wcbs_mod5, plot = TRUE)
vif(Abund_wcbs_mod5)
summary(Abund_wcbs_mod5)

addterm(Abund_wcbs_mod5,
        scope = Abundance ~ AES1KM*AES3KM + N_VISITS_MAYTOAUGUST + YEAR +
          Comp.1 + Comp.2 + Comp.3 + Comp.4 + Comp.5 + Comp.6 +
          Comp.7 + Comp.8 + Comp.9 + Comp.10 + Comp.11 + Comp.12 +
          Comp.13 + Comp.14 + Comp.15 + Comp.16 + Comp.17 + Comp.18 + 
          Comp.19 + Comp.20 + Comp.21 + Comp.22 + Comp.23 + Comp.24 +
          Comp.25 + Comp.26 + Comp.27 + Comp.28, 
        sort = TRUE)
# still big reductions possible but stopping for now

Abund_wcbs_model_results <- do.call(rbind, list(
  data.frame(Variable = names(coef(Abund_wcbs_mod)),
             Coef = coef(Abund_wcbs_mod),
             SE = sqrt(diag(vcov(Abund_wcbs_mod))),
             pval = summary(Abund_wcbs_mod)$coefficients[,4],
             Model = "Abund_wcbs_mod",
             AIC = AIC(Abund_wcbs_mod)),
  data.frame(Variable = names(coef(Abund_wcbs_mod2)),
             Coef = coef(Abund_wcbs_mod2),
             SE = sqrt(diag(vcov(Abund_wcbs_mod2))),
             pval = summary(Abund_wcbs_mod2)$coefficients[,4],
             Model = "Abund_wcbs_mod2",
             AIC = AIC(Abund_wcbs_mod2)),
  data.frame(Variable = names(coef(Abund_wcbs_mod3)),
             Coef = coef(Abund_wcbs_mod3),
             SE = sqrt(diag(vcov(Abund_wcbs_mod3))),
             pval = summary(Abund_wcbs_mod3)$coefficients[,4],
             Model = "Abund_wcbs_mod3",
             AIC = AIC(Abund_wcbs_mod3)),
  data.frame(Variable = names(coef(Abund_wcbs_mod4)),
             Coef = coef(Abund_wcbs_mod4),
             SE = sqrt(diag(vcov(Abund_wcbs_mod4))),
             pval = summary(Abund_wcbs_mod4)$coefficients[,4],
             Model = "Abund_wcbs_mod4",
             AIC = AIC(Abund_wcbs_mod4)),
  data.frame(Variable = names(coef(Abund_wcbs_mod5)),
             Coef = coef(Abund_wcbs_mod5),
             SE = sqrt(diag(vcov(Abund_wcbs_mod5))),
             pval = summary(Abund_wcbs_mod5)$coefficients[,4],
             Model = "Abund_wcbs_mod5",
             AIC = AIC(Abund_wcbs_mod5))
)) %>%
  pivot_longer(Coef:pval) %>%
  pivot_wider(names_from = c(Variable,name),
              names_sep = "_",
              values_from = value,
              values_fill = NA) %>%
  mutate(Family = "Negative binomial")



#' ### Diversity models
Div_wcbs_mod <- lm(exp(Shannon_diversity) ~ AES1KM*AES3KM + N_VISITS_MAYTOAUGUST + YEAR +
                     Comp.1 + Comp.2, 
                   data = WCBS_all_data)

summary(Div_wcbs_mod)
vif(Div_wcbs_mod)
# par(mfrow=c(2,2))
# plot(Div_wcbs_mod)
# par(mfrow=c(1,1))

# Check if should add terms
addterm(Div_wcbs_mod,
        scope = exp(Shannon_diversity) ~ AES1KM*AES3KM + N_VISITS_MAYTOAUGUST + YEAR +
          Comp.1 + Comp.2 + Comp.3 + Comp.4 + Comp.5 + Comp.6 +
          Comp.7 + Comp.8 + Comp.9 + Comp.10 + Comp.11 + Comp.12 +
          Comp.13 + Comp.14 + Comp.15 + Comp.16 + Comp.17 + Comp.18 + 
          Comp.19 + Comp.20 + Comp.21 + Comp.22 + Comp.23 + Comp.24 +
          Comp.25 + Comp.26 + Comp.27 + Comp.28, sort = TRUE)

#' Shannon diversity model with more PCA axes
Div_wcbs_mod2 <- lm(exp(Shannon_diversity) ~ AES1KM*AES3KM + N_VISITS_MAYTOAUGUST + YEAR +
                      Comp.1 + Comp.2 + Comp.15, 
                    data = WCBS_all_data)
summary(Div_wcbs_mod2)
vif(Div_wcbs_mod2)

# par(mfrow=c(2,2));plot(Div_wcbs_mod2);par(mfrow=c(1,1))


addterm(Div_wcbs_mod2,
        scope = exp(Shannon_diversity) ~ AES1KM*AES3KM + N_VISITS_MAYTOAUGUST + YEAR +
          Comp.1 + Comp.2 + Comp.3 + Comp.4 + Comp.5 + Comp.6 +
          Comp.7 + Comp.8 + Comp.9 + Comp.10 + Comp.11 + Comp.12 +
          Comp.13 + Comp.14 + Comp.15 + Comp.16 + Comp.17 + Comp.18 + 
          Comp.19 + Comp.20 + Comp.21 + Comp.22 + Comp.23 + Comp.24 +
          Comp.25 + Comp.26 + Comp.27 + Comp.28, sort = TRUE)


Div_wcbs_mod3 <- lm(exp(Shannon_diversity) ~ AES1KM*AES3KM + N_VISITS_MAYTOAUGUST + YEAR +
                      Comp.1 + Comp.2 + Comp.15 + Comp.13, 
                    data = WCBS_all_data)
summary(Div_wcbs_mod3)
vif(Div_wcbs_mod3)

# par(mfrow=c(2,2));plot(Div_wcbs_mod3);par(mfrow=c(1,1))


addterm(Div_wcbs_mod3,
        scope = exp(Shannon_diversity) ~ AES1KM*AES3KM + N_VISITS_MAYTOAUGUST + YEAR +
          Comp.1 + Comp.2 + Comp.3 + Comp.4 + Comp.5 + Comp.6 +
          Comp.7 + Comp.8 + Comp.9 + Comp.10 + Comp.11 + Comp.12 +
          Comp.13 + Comp.14 + Comp.15 + Comp.16 + Comp.17 + Comp.18 + 
          Comp.19 + Comp.20 + Comp.21 + Comp.22 + Comp.23 + Comp.24 +
          Comp.25 + Comp.26 + Comp.27 + Comp.28, sort = TRUE)


Div_wcbs_mod4 <- lm(exp(Shannon_diversity) ~ AES1KM*AES3KM + N_VISITS_MAYTOAUGUST + YEAR +
                      Comp.1 + Comp.2 + Comp.15 + Comp.13 + Comp.7, 
                    data = WCBS_all_data)
summary(Div_wcbs_mod4)
vif(Div_wcbs_mod4)

# par(mfrow=c(2,2));plot(Div_wcbs_mod4);par(mfrow=c(1,1))


addterm(Div_wcbs_mod4,
        scope = exp(Shannon_diversity) ~ AES1KM*AES3KM + N_VISITS_MAYTOAUGUST + YEAR +
          Comp.1 + Comp.2 + Comp.3 + Comp.4 + Comp.5 + Comp.6 +
          Comp.7 + Comp.8 + Comp.9 + Comp.10 + Comp.11 + Comp.12 +
          Comp.13 + Comp.14 + Comp.15 + Comp.16 + Comp.17 + Comp.18 + 
          Comp.19 + Comp.20 + Comp.21 + Comp.22 + Comp.23 + Comp.24 +
          Comp.25 + Comp.26 + Comp.27 + Comp.28, sort = TRUE)


Div_wcbs_mod5 <- lm(exp(Shannon_diversity) ~ AES1KM*AES3KM + N_VISITS_MAYTOAUGUST + YEAR +
                      Comp.1 + Comp.2 + Comp.15 + Comp.13 + Comp.7 + Comp.3, 
                    data = WCBS_all_data)
summary(Div_wcbs_mod5)
vif(Div_wcbs_mod5)

# par(mfrow=c(2,2));plot(Div_wcbs_mod5);par(mfrow=c(1,1))


div_wcbs_model_results <- do.call(rbind, list(
  data.frame(Variable = names(coef(Div_wcbs_mod)),
             Coef = coef(Div_wcbs_mod),
             SE = sqrt(diag(vcov(Div_wcbs_mod))),
             pval = summary(Div_wcbs_mod)$coefficients[,4],
             Model = "Div_wcbs_mod",
             AIC = AIC(Div_wcbs_mod)),
  data.frame(Variable = names(coef(Div_wcbs_mod2)),
             Coef = coef(Div_wcbs_mod2),
             SE = sqrt(diag(vcov(Div_wcbs_mod2))),
             pval = summary(Div_wcbs_mod2)$coefficients[,4],
             Model = "Div_wcbs_mod2",
             AIC = AIC(Div_wcbs_mod2)),
  data.frame(Variable = names(coef(Div_wcbs_mod3)),
             Coef = coef(Div_wcbs_mod3),
             SE = sqrt(diag(vcov(Div_wcbs_mod3))),
             pval = summary(Div_wcbs_mod3)$coefficients[,4],
             Model = "Div_wcbs_mod3",
             AIC = AIC(Div_wcbs_mod3)),
  data.frame(Variable = names(coef(Div_wcbs_mod4)),
             Coef = coef(Div_wcbs_mod4),
             SE = sqrt(diag(vcov(Div_wcbs_mod4))),
             pval = summary(Div_wcbs_mod4)$coefficients[,4],
             Model = "Div_wcbs_mod4",
             AIC = AIC(Div_wcbs_mod4)),
  data.frame(Variable = names(coef(Div_wcbs_mod5)),
             Coef = coef(Div_wcbs_mod5),
             SE = sqrt(diag(vcov(Div_wcbs_mod5))),
             pval = summary(Div_wcbs_mod5)$coefficients[,4],
             Model = "Div_wcbs_mod5",
             AIC = AIC(Div_wcbs_mod5))
)) %>%
  pivot_longer(Coef:pval) %>%
  pivot_wider(names_from = c(Variable,name),
              names_sep = "_",
              values_from = value,
              values_fill = NA) %>%
  mutate(Family = "Normal")


wcbs_model_results <- full_join(Rich_wcbs_model_results, Abund_wcbs_model_results) %>%
  full_join(div_wcbs_model_results) %>%
  mutate(Survey = "WCBS",
         Taxa = "Butterfly",
         Response = sapply(strsplit(Model, "_"),"[",1)) %>%
  mutate(Response = recode(Response,
                           "Rich" = "Richness",
                           "Div" = "Shannon diversity (exp)",
                           "Abund" = "Abundance"))
# Plots ####
# create dummy data for prediction
# pred_data <- data.frame(
#   AES1KM = rep(seq(from = min(WCBS_all_data$AES1KM, na.rm=TRUE),
#                    to = max(WCBS_all_data$AES1KM, na.rm=TRUE),
#                    length.out = 100), 3),
#   AES3KM = c(rep(-0.74,100),rep(-0.25,100),rep(1.4,100)),
#   YEAR = factor(2018, levels = levels(WCBS_all_data$YEAR)),
#   N_VISITS_MAYTOAUGUST = 4, Comp.1 = 0, Comp.2 = 0, 
#   Comp.26 = 0, Comp.6 = 0, Comp.13 = 0, Comp.23 = 0, 
#   Comp.18 = 0, Comp.11 = 0, Comp.25 = 0, Comp.15 = 0,
#   Comp.8 = 0)
# 
# # predict richness values
# pred_data_richness <- cbind(pred_data, predict(Rich_wcbs_mod5, pred_data, type = "link",
#                                                se.fit=TRUE)) %>%
#   mutate(LL = exp(fit - se.fit),
#          UL = exp(fit + se.fit)) %>%
#   mutate(AES3KM = recode(AES3KM,
#                          `1.4` = "High",
#                          `-0.25` = "Medium",
#                          `-0.74` = "Low"),
#          fit = exp(fit),
#          AES1KM = 6301*AES1KM + 3110.32)
# 
# # plot richness values with data
# wcbs_r_plot <- WCBS_all_data %>%
#   mutate(AES3KM = cut(AES3KM, c(-Inf,-0.66,0.39,Inf),
#                       c("Low","Medium","High")),
#          AES1KM = 6301*AES1KM + 3110.32) %>%
#   filter(!is.na(AES3KM)) %>%
#   ggplot(aes(x = AES1KM, y = Richness, colour = AES3KM, fill = AES3KM)) +
#   geom_point() +
#   geom_ribbon(data = pred_data_richness,
#               aes(x = AES1KM, y = fit, ymin = LL, ymax = UL,
#                   fill = AES3KM, colour = NULL), alpha = .25) +
#   geom_line(data = pred_data_richness,
#             aes(x = AES1KM, y = fit, colour = AES3KM)) +
#   labs(x = "AES 1km", y = "Butterfly richness (WCBS)")
# 
# 
# # predict abundance values
# pred_data_abundance <- cbind(pred_data, predict(Abund_wcbs_mod5, pred_data, type = "link",
#                                                 se.fit=TRUE)) %>%
#   mutate(LL = exp(fit - se.fit),
#          UL = exp(fit + se.fit)) %>%
#   mutate(AES3KM = recode(AES3KM,
#                          `1.4` = "High",
#                          `-0.25` = "Medium",
#                          `-0.74` = "Low"),
#          fit = exp(fit),
#          AES1KM = 6301*AES1KM + 3110.32)
# 
# # plot abundance predictions with data
# wcbs_a_plot <- WCBS_all_data %>%
#   mutate(AES3KM = cut(AES3KM, c(-Inf,-0.66,0.39,Inf),
#                       c("Low","Medium","High")),
#          AES1KM = 6301*AES1KM + 3110.32) %>%
#   filter(!is.na(AES3KM)) %>%
#   ggplot(aes(x = AES1KM, y = Abundance, colour = AES3KM, fill = AES3KM)) +
#   geom_point() +
#   geom_ribbon(data = pred_data_abundance,
#               aes(x = AES1KM, y = fit, ymin = LL, ymax = UL,
#                   fill = AES3KM, colour = NULL), alpha = .25) +
#   geom_line(data = pred_data_abundance,
#             aes(x = AES1KM, y = fit, colour = AES3KM)) +
#   labs(x = "AES 1km", y = "Butterfly abundance (WCBS)")
# 
# # predict Shannon diversity values for different AES scores
# pred_data_diversity <- cbind(pred_data, predict(Div_wcbs_mod5, pred_data, type = "response",
#                                                 se.fit=TRUE)) %>%
#   mutate(LL = log(fit - se.fit),
#          UL = log(fit + se.fit)) %>%
#   mutate(AES3KM = recode(AES3KM,
#                          `1.4` = "High",
#                          `-0.25` = "Medium",
#                          `-0.74` = "Low"),
#          fit = log(fit),
#          AES1KM = 6301*AES1KM + 3110.32)
# 
# # plot Shannon diversity predictions with data
# wcbs_d_plot <- WCBS_all_data %>%
#   mutate(AES3KM = cut(AES3KM, c(-Inf,-0.66,0.39,Inf),
#                       c("Low","Medium","High")),
#          AES1KM = 6301*AES1KM + 3110.32) %>%
#   filter(!is.na(AES3KM)) %>%
#   ggplot(aes(x = AES1KM, y = Shannon_diversity, colour = AES3KM, fill = AES3KM)) +
#   geom_point() +
#   geom_ribbon(data = pred_data_diversity,
#               aes(x = AES1KM, y = fit, ymin = LL, ymax = UL,
#                   fill = AES3KM, colour = NULL), alpha = .25) +
#   geom_line(data = pred_data_diversity,
#             aes(x = AES1KM, y = fit, colour = AES3KM)) +
#   labs(x = "AES 1km", y = "Butterfly diversity (WCBS)")
