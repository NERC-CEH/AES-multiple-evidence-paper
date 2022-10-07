#' ## Butterfly integrated models

library(gtools)
library(lme4)
library(brms)
library(scales)


#' source individual scheme models

source("AES_Models_UKBMS.R")

source("AES_Models_WCBS.R")

source("LandSpaes_analysis/LandSpaes_butterfly_env_data_models.R")

source("validation_function.R")




#try idea of z test to test similarity
ztest <- function(mod1, mod2, term){
  p <- vector(); z <- vector()
  for (i in term){
  b1 <- mod1$coefficients[names(mod1$coefficients) == i]
  b2 <- mod2$coefficients[names(mod2$coefficients) == i]
  se1 <- sqrt(diag(vcov(mod1)))[names(sqrt(diag(vcov(mod1)))) == i]
  se2 <- sqrt(diag(vcov(mod2)))[names(sqrt(diag(vcov(mod2)))) == i]
  z[which(term == i)] <- (b1-b2)/(sqrt((se1^2) + (se2^2)))
  p[which(term == i)] <- 2*pnorm(abs(z[which(term == i)]), lower.tail = FALSE)
  }
  z.df <- data.frame(term = term, z = z, p = p)
  return(z.df)
}

#richness
ztest(buttrich_mod1, Rich_wcbs_mod5, c("AES1KM", "AES3KM", "AES1KM:AES3KM"))#OK
ztest(buttrich_mod1, Rich_ukbms_mod5, c("AES1KM", "AES3KM", "AES1KM:AES3KM"))#OK
ztest(Rich_ukbms_mod5, Rich_wcbs_mod5, c("AES1KM", "AES3KM", "AES1KM:AES3KM"))#3km and int diff

#abundance
ztest(buttabund_mod2, Abund_wcbs_mod5, c("AES1KM", "AES3KM", "AES1KM:AES3KM"))#OK
ztest(buttabund_mod2, Abund_ukbms_mod5, c("AES1KM", "AES3KM", "AES1KM:AES3KM"))#diff 3km
ztest(Abund_ukbms_mod5, Abund_wcbs_mod5, c("AES1KM", "AES3KM", "AES1KM:AES3KM"))#diff 1 and 3km

#diversity
ztest(buttdiv_mod2, Div_wcbs_mod5, c("AES1KM", "AES3KM", "AES1KM:AES3KM"))#borderline OK
ztest(buttdiv_mod2, Div_ukbms_mod5, c("AES1KM", "AES3KM", "AES1KM:AES3KM"))#diff 1km
ztest(Div_ukbms_mod5, Div_wcbs_mod5, c("AES1KM", "AES3KM", "AES1KM:AES3KM"))#OK

#richness and abundance definitely OK, diversity borderline
#may need to remove UKBMS from abundance and diversity?






#' Overall similar coefficients, although slopes are smaller in UKBMS 
#' 
#' Intercepts are fairly similar, slightly lower in WCBS
#' 
#' Terms in integrated model:
#' 
#' - Richness (or RICHNESS_ID)
#' - AES1KM
#' - AES3KM
#' - N_VISITS_MAYTOAUGUST (or ROUND_NUMBER)
#' - YEAR (or SURVEY_YEAR)
#' - PCA axes (specific only to each submodel)
#' 
#' 

#' ## Richness

#' Compare results from individual scheme models

#' LandSpAES model

summary(buttrich_mod1)

#' WCBS model

summary(Rich_wcbs_mod5)

#' UKBMS model

summary(Rich_ukbms_mod5)

#for consistency rename LM0465 data columns
names(buttrich_env)[c(1,30,31,32,33,34,35)] <- c("SURVEY_SQUARE", "YEAR", "NCA", "AES1KM", "AES3KM", "Richness", "N_VISITS_MAYTOAUGUST")

#' For each scheme, set the value of the PCA axes for components not selected in the individual scheme models to NA. This *should* mean that axes relationships will be shared across schemes only if the same axes are selected. 
#' 
#' Checked individual schemes models that no strongly constrasting patterns in relation to PCA axes - overall where scores are shared the relationships are broadly consistent between schemes. Doesn't matter too much if there are slight differences as we are not looking to reduce uncertainty in these estimated relationships

WCBS_int_data <- as.data.frame(WCBS_all_data)

#WCBS_int_data[,grepl("Comp.", names(WCBS_int_data)) & ! names(WCBS_int_data) %in%  names(Rich_wcbs_mod5$coefficients)] <- NA
#WCBS_int_data <- WCBS_int_data[,which(names(WCBS_int_data) %in%  c("Richness", "YEAR", names(Rich_wcbs_mod5$coefficients)))]
WCBS_int_data$SURVEY <- "WCBS"

LM0465_int_data <- buttrich_env

#LM0465_int_data[,grepl("Comp.", names(LM0465_int_data)) & ! names(LM0465_int_data) %in%  names(buttrich_mod1$coefficients)] <- NA
#LM0465_int_data <- LM0465_int_data[,which(names(LM0465_int_data) %in%  c("Richness", "YEAR", "N_VISITS_MAYTOAUGUST", names(buttrich_mod1$coefficients)))]
LM0465_int_data$SURVEY <- "LandSpAES"

UKBMS_int_data <- as.data.frame(UKBMS_all_data)

#UKBMS_int_data[,grepl("Comp.", names(UKBMS_int_data)) & ! names(UKBMS_int_data) %in%  names(Rich_ukbms_mod5$coefficients)] <- NA
#UKBMS_int_data <- UKBMS_int_data[,which(names(UKBMS_int_data) %in%  c("Richness", "YEAR", names(Rich_ukbms_mod5$coefficients)))]
UKBMS_int_data$SURVEY <- "UKBMS"


#integrated model


int_dat <- smartbind(LM0465_int_data, WCBS_int_data, UKBMS_int_data)

int_dat$SURVEY <- as.factor(int_dat$SURVEY)
int_dat$YEAR <- as.factor(int_dat$YEAR)

int_dat$N_VISITS_MAYTOAUGUSTSC <- rescale(int_dat$N_VISITS_MAYTOAUGUST)

#use rescale instead of scale to tackle convergence issues
#int_dat[,2:29] <- scale(int_dat[,2:29])
for (i in 2:29){
  int_dat[i] <- rescale(int_dat[,i])}

#int_dat$AES1KM <- rescale(int_dat$AES1KM)
#int_dat$AES3KM <- rescale(int_dat$AES3KM)

#check nb for LM0465 (not chosen as overdipersion not evident) - seems fine
# buttrich_mod1nb <-  glm.nb(formula = Richness ~ AES1KM * AES3KM + N_VISITS_MAYTOAUGUST + 
#         YEAR + Comp.1 + Comp.2 + Comp.6 + Comp.11 + Comp.15 + Comp.21 + Comp.23, data = buttrich_env)
#   
# int_mod1 <- glmer.nb(Richness ~ AES1KM * AES3KM  + N_VISITS_MAYTOAUGUST + YEAR + Comp.1 + Comp.2 + Comp.6 + Comp.11 + Comp.15 + Comp.21 + Comp.23 + (AES1KM*AES3KM|SURVEY) , data = int_dat)
#   
# summary(int_mod1)

#look at model where only random intercept included

#int_mod2 <- glmer.nb(Richness ~ AES1KM * AES3KM  + N_VISITS_MAYTOAUGUST + YEAR + Comp.1 + Comp.2 + Comp.6 + Comp.11 + Comp.15 + Comp.21 + Comp.23 + (1|SURVEY) , data = int_dat)

#summary(int_mod2)

#poisson

buttrich_int_mod1 <- glmer(Richness ~ AES1KM * AES3KM  + N_VISITS_MAYTOAUGUSTSC + YEAR + Comp.1 + Comp.2 + Comp.6 + Comp.11 + Comp.15 + Comp.21 + Comp.23 + (1|SURVEY) , data = int_dat, family = "poisson")

summary(buttrich_int_mod1)

source("validation_function.R")

v_rich_int <- val_fun(buttrich_int_mod1, type = "glmer")

v_rich <- val_fun(buttrich_mod1, type = "glm")

# library(glmmTMB)
# 
# int_mod2bTMB <- glmmTMB(Richness ~ AES1KM * AES3KM  + N_VISITS_MAYTOAUGUST + YEAR + Comp.1 + Comp.2 + Comp.6 + Comp.11 + Comp.15 + Comp.21 + Comp.23 + (1|SURVEY) , data = int_dat, family = "poisson")
# 
# summary(int_mod2bTMB)

#try with just WCBS



int_dat <- smartbind(LM0465_int_data, WCBS_int_data)

int_dat$SURVEY <- as.factor(int_dat$SURVEY)
int_dat$YEAR <- as.factor(int_dat$YEAR)

int_dat$N_VISITS_MAYTOAUGUSTSC <- rescale(int_dat$N_VISITS_MAYTOAUGUST)

#use rescale instead of scale to tackle convergence issues
#int_dat[,2:29] <- scale(int_dat[,2:29])
for (i in 2:29){
  int_dat[i] <- rescale(int_dat[,i])}


buttrich_int_mod2 <- glmer(Richness ~ AES1KM * AES3KM  + N_VISITS_MAYTOAUGUSTSC + YEAR + Comp.1 + Comp.2 + Comp.6 + Comp.11 + Comp.15 + Comp.21 + Comp.23 + (1|SURVEY) , data = int_dat, family = "poisson")

summary(buttrich_int_mod2)

source("validation_function.R")

v_rich_int2 <- val_fun(buttrich_int_mod2, type = "glmer")

v_rich2 <- val_fun(buttrich_mod2, type = "glm")



#simplify

# int_mod2a <- glmer.nb(Richness ~ AES1KM*AES3KM  + (1|SURVEY) , data = int_dat)
# 
# summary(int_mod2a)


#int_mod1 <- glm.nb(Richness ~ AES1KM * AES3KM + N_VISITS_MAYTOAUGUST + YEAR + Comp.1 + Comp.2 + Comp.6 + Comp.15 + Comp.21 + Comp.23 + SURVEY, data = int_dat)

#Bayesian fitting using Stan and brms
# 
# int_mod3 <- brm(formula = Richness ~ AES1KM * AES3KM + N_VISITS_MAYTOAUGUST + 
#                   YEAR + Comp.1 + Comp.2 + Comp.6 + Comp.11 + Comp.15 + Comp.21 + Comp.23 + (AES1KM*AES3KM|SURVEY), data = int_dat, family = negbinomial(link = "log", link_shape = "log"), refresh = 0)
# 
# summary(int_mod3)
  


###Abundance - just WCBS

#for consistency rename LM0465 data columns
names(buttabund_env_vars)[c(1:7)] <- c("SURVEY_SQUARE", "YEAR", "NCA", "AES1KM", "AES3KM", "Abundance", "N_VISITS_MAYTOAUGUST")

#' For each scheme, set the value of the PCA axes for components not selected in the individual scheme models to NA. This *should* mean that axes relationships will be shared across schemes only if the same axes are selected. 
#' 
#' Checked individual schemes models that no strongly constrasting patterns in relation to PCA axes - overall where scores are shared the relationships are broadly consistent between schemes. Doesn't matter too much if there are slight differences as we are not looking to reduce uncertainty in these estimated relationships

WCBS_int_data <- as.data.frame(WCBS_all_data)

#WCBS_int_data[,grepl("Comp.", names(WCBS_int_data)) & ! names(WCBS_int_data) %in%  names(Rich_wcbs_mod5$coefficients)] <- NA
#WCBS_int_data <- WCBS_int_data[,which(names(WCBS_int_data) %in%  c("Richness", "YEAR", names(Rich_wcbs_mod5$coefficients)))]
WCBS_int_data$SURVEY <- "WCBS"

LM0465_int_data <- buttabund_env_vars

#LM0465_int_data[,grepl("Comp.", names(LM0465_int_data)) & ! names(LM0465_int_data) %in%  names(buttrich_mod1$coefficients)] <- NA
#LM0465_int_data <- LM0465_int_data[,which(names(LM0465_int_data) %in%  c("Richness", "YEAR", "N_VISITS_MAYTOAUGUST", names(buttrich_mod1$coefficients)))]
LM0465_int_data$SURVEY <- "LandSpAES"


#integrated model


int_dat <- smartbind(LM0465_int_data, WCBS_int_data)

int_dat$SURVEY <- as.factor(int_dat$SURVEY)
int_dat$YEAR <- as.factor(int_dat$YEAR)

int_dat$N_VISITS_MAYTOAUGUSTSC <- rescale(int_dat$N_VISITS_MAYTOAUGUST)

#use rescale instead of scale to tackle convergence issues
#int_dat[,2:29] <- scale(int_dat[,2:29])
for (i in 8:35){
 int_dat[i] <- rescale(int_dat[,i])}

int_dat$AES1KMSc <- rescale(int_dat$AES1KM)
int_dat$AES3KMSc <- rescale(int_dat$AES3KM)

# AES1KM3 <- (int_dat$AES1KM - min(int_dat$AES1KM))/(max(int_dat$AES1KM) - min(int_dat$AES1KM))
# 
# AES1KM4 <- AES1KM3*(max(int_dat$AES1KM) - min(int_dat$AES1KM)) + min(int_dat$AES1KM)

#' LandSpAES model

summary(buttabund_mod2)

#' WCBS model

summary(Abund_wcbs_mod5)



buttabund_int_mod1 <- glmer.nb(Abundance ~ AES1KM * AES3KM  + N_VISITS_MAYTOAUGUSTSC + YEAR + Comp.1 + Comp.2 + Comp.9 + Comp.19 + Comp.6 + Comp.15 + Comp.20 + (1|SURVEY) , data = int_dat)

summary(buttabund_int_mod1)

# buttabund_int_mod2 <- glmer.nb(Abundance ~ AES1KM * AES3KM  + N_VISITS_MAYTOAUGUST + YEAR + Comp.1 + Comp.2 + Comp.9 + Comp.19 + Comp.6 + Comp.15 + Comp.20 + (1|SURVEY) , data = int_dat)
# 
# summary(buttabund_int_mod2)
# 
# coef(buttabund_int_mod1)$SURVEY[1,2] * (max(int_dat$AES1KM) - min(int_dat$AES1KM)) + min(int_dat$AES1KM)

# buttabund_TMB <- glmmTMB(Abundance ~ AES1KM * AES3KM  + N_VISITS_MAYTOAUGUST + YEAR + Comp.1 + Comp.2 + Comp.9 + Comp.19 + Comp.6 + Comp.15 + Comp.20 + (1|SURVEY) , data = int_dat, family = "nbinom2")
# 
# summary(buttabund_TMB)


v_abund_int <- val_fun(buttabund_int_mod1, type = "glmer.nb")

v_abund <- val_fun(buttabund_mod2, type = "glm.nb")

v_abund_int
v_abund


###Diversity - all 3

#for consistency rename LM0465 data columns
names(buttdiv_env)[c(1:7)] <- c("SURVEY_SQUARE", "YEAR", "NCA", "AES1KM", "AES3KM", "Shannon_diversity", "N_VISITS_MAYTOAUGUST")

#' For each scheme, set the value of the PCA axes for components not selected in the individual scheme models to NA. This *should* mean that axes relationships will be shared across schemes only if the same axes are selected. 
#' 
#' Checked individual schemes models that no strongly constrasting patterns in relation to PCA axes - overall where scores are shared the relationships are broadly consistent between schemes. Doesn't matter too much if there are slight differences as we are not looking to reduce uncertainty in these estimated relationships

WCBS_int_data <- as.data.frame(WCBS_all_data)

#WCBS_int_data[,grepl("Comp.", names(WCBS_int_data)) & ! names(WCBS_int_data) %in%  names(Rich_wcbs_mod5$coefficients)] <- NA
#WCBS_int_data <- WCBS_int_data[,which(names(WCBS_int_data) %in%  c("Richness", "YEAR", names(Rich_wcbs_mod5$coefficients)))]
WCBS_int_data$SURVEY <- "WCBS"

LM0465_int_data <- buttdiv_env

#LM0465_int_data[,grepl("Comp.", names(LM0465_int_data)) & ! names(LM0465_int_data) %in%  names(buttrich_mod1$coefficients)] <- NA
#LM0465_int_data <- LM0465_int_data[,which(names(LM0465_int_data) %in%  c("Richness", "YEAR", "N_VISITS_MAYTOAUGUST", names(buttrich_mod1$coefficients)))]
LM0465_int_data$SURVEY <- "LandSpAES"

UKBMS_int_data <- as.data.frame(UKBMS_all_data)

#UKBMS_int_data[,grepl("Comp.", names(UKBMS_int_data)) & ! names(UKBMS_int_data) %in%  names(Rich_ukbms_mod5$coefficients)] <- NA
#UKBMS_int_data <- UKBMS_int_data[,which(names(UKBMS_int_data) %in%  c("Richness", "YEAR", names(Rich_ukbms_mod5$coefficients)))]
UKBMS_int_data$SURVEY <- "UKBMS"


#integrated model


int_dat <- smartbind(LM0465_int_data, WCBS_int_data, UKBMS_int_data)

int_dat$SURVEY <- as.factor(int_dat$SURVEY)
int_dat$YEAR <- as.factor(int_dat$YEAR)

int_dat$N_VISITS_MAYTOAUGUSTSC <- rescale(int_dat$N_VISITS_MAYTOAUGUST)

#use rescale instead of scale to tackle convergence issues
#int_dat[,2:29] <- scale(int_dat[,2:29])
#for (i in 2:29){
#  int_dat[i] <- rescale(int_dat[,i])}

#int_dat$AES1KM <- rescale(int_dat$AES1KM)
#int_dat$AES3KM <- rescale(int_dat$AES3KM)




#' LandSpAES model

summary(buttdiv_mod2)

#' WCBS model

summary(Div_wcbs_mod5)

#' UKBMS model

summary(Div_ukbms_mod5)


# 
# buttdiv_int_mod1 <- lme(exp(Shannon_diversity) ~ AES1KM * AES3KM  + N_VISITS_MAYTOAUGUSTSC + YEAR + Comp.1 + Comp.2 + Comp.13 + Comp.15 + Comp.7 + Comp.3 + Comp.21, random = ~1|SURVEY , data = int_dat, na.action = na.omit)
# 
# summary(buttdiv_int_mod1)

buttdiv_int_mod1 <- lme(exp(Shannon_diversity) ~ AES1KM * AES3KM  + N_VISITS_MAYTOAUGUSTSC + YEAR + Comp.1 + Comp.2 + Comp.13 + Comp.15 + Comp.7 + Comp.3 + Comp.21, random = ~AES1KM*AES3KM|SURVEY , data = int_dat, na.action = na.omit)

summary(buttdiv_int_mod1)


v_div_int <- val_fun(buttdiv_int_mod1, type = "lme")

v_div <- val_fun(buttdiv_mod2, type = "lm")

v_div_int
v_div


val_out <- rbind(v_rich, v_rich_int, v_rich2, v_rich_int2, v_abund, v_abund_int, v_div, v_div_int)
val_out$Response <- c("Richness", "Richness","RichnessWCBSonly", "RichnessWCBSonly", "Abundance", "Abundance", "Diversity", "Diversity")
val_out$Model <- c("LandSpAES", "Integrated","LandSpAES", "Integrated", "LandSpAES", "Integrated","LandSpAES", "Integrated")

write.csv(val_out, "Validation butterfly models.csv")

#### Plots


#check comparable visit nos

(4 - min(int_dat$N_VISITS_MAYTOAUGUST ))/(max(int_dat$N_VISITS_MAYTOAUGUST ) - min(int_dat$N_VISITS_MAYTOAUGUST ))


library(ggeffects)
library(ggplot2)
library(gridExtra)
library(patchwork)
theme_set(theme_classic())
my_col <- unname(palette.colors()[c(8,9)])


##Richness



p1 <- ggpredict(buttrich_mod1, "AES1KM[-0.48:5.02, by = 0.5]",
                condition = c("AES3KM" = "0.01",
                              "ROUND_NUMBER" = "4",
                              "SURVEY_YEAR" = "2018",
                              "Comp.1" = "0",
                              "Comp.2" = "0"))


p1_int <- ggpredict(buttrich_int_mod1, "AES1KM[-0.48:6.52, by = 0.5]",
                    condition = c("AES3KM" = "0.01",
                                  "N_VISITS_MAYTOAUGUSTSC" = "0.0461",
                                  "YEAR" = "2018"
                                  ), type = "fixed")

p1$group <- "LandSpAES"
p1_int$group <- "Integrated"

#get raw data for rug
p1_raw <- attr(p1, "rawdata") %>%
  mutate(x = (x*6966.606 + 3340.947)/1000,
         group = "LandSpAES")
p1_int_raw <- attr(p1_int, "rawdata") %>%
  mutate(x = (x*6966.606 + 3340.947)/1000,
         group = "Integrated")


p4 <- rbind(p1, p1_int) %>%
  mutate(x = (x*6966.606 + 3340.947)/1000,
         group = forcats::fct_inorder(as.factor(group))) # orders factor in order it appears


richint_1km <- ggplot(p4, aes(x = x, y = predicted, colour = group, fill = group)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.3, colour = NA) +
  geom_line() + 
  geom_rug(data = p1_raw, aes(x = x, y = response, colour = group),
           sides = "b") +
  geom_rug(data = p1_int_raw, aes(x = x, y = response, colour = group),
           sides = "t") +
  coord_cartesian(clip = "off") +
  scale_fill_manual(aesthetics = c("fill","colour"),
                    values = my_col,
                    name = "Survey") +
  scale_y_continuous(limits = c(5,30), expand = c(0,0))+
  labs(x = "AES 1km ('000s)", y = "Predicted Butterfly Richness") +
  NULL

richint_1km

ggsave("report/Butterfly richness 1km integrated model.png", height = 800, width = 1000, units = "mm", scale = 0.15)

p2_int <- ggpredict(buttrich_int_mod1, "AES3KM[-0.86:9.64, by = 0.5]",
                    condition = c("AES1KM" = "0",
                                  "N_VISITS_MAYTOAUGUSTSC" = "0.0461",
                                  "YEAR" = "2018"), type = "fixed")


p2 <- ggpredict(buttrich_mod1, "AES3KM[-0.86:5.14, by = 0.5]",
                condition = c("AES1KM" = "0",
                              "ROUND_NUMBER" = "4",
                              "SURVEY_YEAR" = "2018"))

p2$group <- "LandSpAES"
p2_int$group <- "Integrated"

#get raw data for rug
p2_raw <- attr(p2, "rawdata") %>%
  mutate(x = (x*3745.525 + 3245.172)/1000,
         group = "LandSpAES")
p2_int_raw <- attr(p2_int, "rawdata") %>%
  mutate(x = (x*3745.525 + 3245.172)/1000,
         group = "Integrated")



p4 <- rbind(p2, p2_int) %>%
  mutate(x = (x*3745.525 + 3245.172)/1000,
         group = forcats::fct_inorder(as.factor(group))) # orders factor in order it appears


richint_3km <- ggplot(p4, aes(x = x, y = predicted, colour = group, fill = group)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.3, colour = NA) +
  geom_line() + 
  geom_rug(data = p2_raw, aes(x = x, y = response, colour = group),
           sides = "b") +
  geom_rug(data = p2_int_raw, aes(x = x, y = response, colour = group),
           sides = "t") +
  coord_cartesian(clip = "off") +
  scale_fill_manual(aesthetics = c("fill","colour"),
                    values = my_col,
                    name = "Survey") +
  scale_y_continuous(limits = c(5,30), expand = c(0,0))+
  labs(x = "AES 3km ('000s)", y = "Predicted Butterfly Richness") +
  NULL

richint_3km

ggsave("report/Butterfly richness 3km integrated model.png", height = 800, width = 1000, units = "mm", scale = 0.15)

richint_1km + richint_3km + plot_layout(guides='collect') &
  theme(legend.position='bottom')

ggsave("report/Butterfly richness combined 1km and 3km integrated model.png", height = 800, width = 1200, units = "mm", scale = 0.15)

p3_int <- ggpredict(buttrich_int_mod1, c("AES1KM[-0.48:6.52, by = 0.5]",
                                         "AES3KM[-0.8,-0.199,1.803]"),
                    condition = c("N_VISITS_MAYTOAUGUSTSC" = "0.0461",
                                  "YEAR" = "2018"), type = "fixed")

p3 <- ggpredict(buttrich_mod1, c("AES1KM[-0.48:5.02, by = 0.5]",
                                         "AES3KM[-0.8,-0.199,1.803]"),
                    condition = c("ROUND_NUMBER" = "4",
                                  "SURVEY_YEAR" = "2018"), type = "fixed")

p3_int$AES3KM <- rep(c(-0.8,-0.199,1.803),length.out = nrow(p3_int))
p3$AES3KM <- rep(c(-0.8,-0.199,1.803),length.out = nrow(p3))

p3$group <- "LandSpAES"
p3_int$group <- "Integrated"

p4 <- rbind(p3, p3_int) %>%
  mutate(x = (x*6966.606 + 3340.947)/1000,
         AES3KM = as.factor(paste("AES 3km:",
                                  round(AES3KM*3745.525 + 3245.172, -1)))) %>%
  mutate(AES3KM = forcats::fct_inorder(AES3KM),
         group = forcats::fct_inorder(as.factor(group)))


ggplot(p4, aes(x = x, y = predicted, colour = group, fill = group)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.3, colour = NA) +
  geom_line() + 
  scale_fill_manual(aesthetics = c("fill","colour"),
                    values = my_col,
                    name = "Survey") +
  facet_wrap(~AES3KM)+
  scale_y_continuous(limits = c(0,40), expand = c(0,0))+
  labs(x = "AES 1km ('000s)", y = "Predicted Butterfly Richness") +
  NULL

ggsave("report/Butterfly richness 1km x 3km integrated model.png", height = 800, width = 1400, units = "mm", scale = 0.15)


##Richness - WCBS only



p1 <- ggpredict(buttrich_mod1, "AES1KM[-0.48:5.02, by = 0.5]",
                condition = c("AES3KM" = "0.01",
                              "ROUND_NUMBER" = "4",
                              "SURVEY_YEAR" = "2018",
                              "Comp.1" = "0",
                              "Comp.2" = "0"))


p1_int <- ggpredict(buttrich_int_mod2, "AES1KM[-0.48:5.52, by = 0.5]",
                    condition = c("AES3KM" = "0.01",
                                  "N_VISITS_MAYTOAUGUSTSC" = "0.4285714",
                                  "YEAR" = "2018"
                    ), type = "fixed")

p1$group <- "LandSpAES"
p1_int$group <- "Integrated"

#get raw data for rug
p1_raw <- attr(p1, "rawdata") %>%
  mutate(x = (x*6966.606 + 3340.947)/1000,
         group = "LandSpAES")
p1_int_raw <- attr(p1_int, "rawdata") %>%
  mutate(x = (x*6966.606 + 3340.947)/1000,
         group = "Integrated")


p4 <- rbind(p1, p1_int) %>%
  mutate(x = (x*6966.606 + 3340.947)/1000,
         group = forcats::fct_inorder(as.factor(group))) # orders factor in order it appears


richint_1km <- ggplot(p4, aes(x = x, y = predicted, colour = group, fill = group)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.3, colour = NA) +
  geom_line() + 
  geom_rug(data = p1_raw, aes(x = x, y = response, colour = group),
           sides = "b") +
  geom_rug(data = p1_int_raw, aes(x = x, y = response, colour = group),
           sides = "t") +
  coord_cartesian(clip = "off") +
  scale_fill_manual(aesthetics = c("fill","colour"),
                    values = my_col,
                    name = "Survey") +
  scale_y_continuous(limits = c(5,30), expand = c(0,0))+
  labs(x = "AES 1km ('000s)", y = "Predicted Butterfly Richness") +
  NULL

richint_1km

ggsave("report/Butterfly richness 1km integrated model - WCBS only.png", height = 800, width = 1000, units = "mm", scale = 0.15)

p2_int <- ggpredict(buttrich_int_mod2, "AES3KM[-0.86:9.64, by = 0.5]",
                    condition = c("AES1KM" = "0",
                                  "N_VISITS_MAYTOAUGUSTSC" = "0.4285714",
                                  "YEAR" = "2018"), type = "fixed")


p2 <- ggpredict(buttrich_mod1, "AES3KM[-0.86:5.14, by = 0.5]",
                condition = c("AES1KM" = "0",
                              "ROUND_NUMBER" = "4",
                              "SURVEY_YEAR" = "2018"))

p2$group <- "LandSpAES"
p2_int$group <- "Integrated"

#get raw data for rug
p2_raw <- attr(p2, "rawdata") %>%
  mutate(x = (x*3745.525 + 3245.172)/1000,
         group = "LandSpAES")
p2_int_raw <- attr(p2_int, "rawdata") %>%
  mutate(x = (x*3745.525 + 3245.172)/1000,
         group = "Integrated")



p4 <- rbind(p2, p2_int) %>%
  mutate(x = (x*3745.525 + 3245.172)/1000,
         group = forcats::fct_inorder(as.factor(group))) # orders factor in order it appears


richint_3km <- ggplot(p4, aes(x = x, y = predicted, colour = group, fill = group)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.3, colour = NA) +
  geom_line() + 
  geom_rug(data = p2_raw, aes(x = x, y = response, colour = group),
           sides = "b") +
  geom_rug(data = p2_int_raw, aes(x = x, y = response, colour = group),
           sides = "t") +
  coord_cartesian(clip = "off") +
  scale_fill_manual(aesthetics = c("fill","colour"),
                    values = my_col,
                    name = "Survey") +
  scale_y_continuous(limits = c(5,30), expand = c(0,0))+
  labs(x = "AES 3km ('000s)", y = "Predicted Butterfly Richness") +
  NULL

richint_3km

ggsave("report/Butterfly richness 3km integrated model - WCBS only.png", height = 800, width = 1000, units = "mm", scale = 0.15)

richint_1km + richint_3km + plot_layout(guides='collect') &
  theme(legend.position='bottom')

ggsave("report/Butterfly richness combined 1km and 3km integrated model - WCBS only.png", height = 800, width = 1200, units = "mm", scale = 0.15)

p3_int <- ggpredict(buttrich_int_mod2, c("AES1KM[-0.48:5.52, by = 0.5]",
                                         "AES3KM[-0.8,-0.199,1.803]"),
                    condition = c("N_VISITS_MAYTOAUGUSTSC" = "0.4285714",
                                  "YEAR" = "2018"), type = "fixed")

p3 <- ggpredict(buttrich_mod1, c("AES1KM[-0.48:5.02, by = 0.5]",
                                 "AES3KM[-0.8,-0.199,1.803]"),
                condition = c("ROUND_NUMBER" = "4",
                              "SURVEY_YEAR" = "2018"), type = "fixed")

p3_int$AES3KM <- rep(c(-0.8,-0.199,1.803),length.out = nrow(p3_int))
p3$AES3KM <- rep(c(-0.8,-0.199,1.803),length.out = nrow(p3))

p3$group <- "LandSpAES"
p3_int$group <- "Integrated"

p4 <- rbind(p3, p3_int) %>%
  mutate(x = (x*6966.606 + 3340.947)/1000,
         AES3KM = as.factor(paste("AES 3km:",
                                  round(AES3KM*3745.525 + 3245.172, -1)))) %>%
  mutate(AES3KM = forcats::fct_inorder(AES3KM),
         group = forcats::fct_inorder(as.factor(group)))


ggplot(p4, aes(x = x, y = predicted, colour = group, fill = group)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.3, colour = NA) +
  geom_line() + 
  scale_fill_manual(aesthetics = c("fill","colour"),
                    values = my_col,
                    name = "Survey") +
  facet_wrap(~AES3KM)+
  scale_y_continuous(limits = c(0,40), expand = c(0,0))+
  labs(x = "AES 1km ('000s)", y = "Predicted Butterfly Richness") +
  NULL

ggsave("report/Butterfly richness 1km x 3km integrated model - WCBS only.png", height = 800, width = 1400, units = "mm", scale = 0.15)



###Abundance


p1 <- ggpredict(buttabund_mod1, "AES1KM[-0.48:5.02, by = 0.5]",
                condition = c("AES3KM" = "0.01",
                              "ROUND_NUMBER" = "4",
                              "SURVEY_YEAR" = "2018"))


p1_int <- ggpredict(buttabund_int_mod1, "AES1KM[-0.48:5.52, by = 0.5]",
                    condition = c("AES3KM" = "0.01",
                                  "N_VISITS_MAYTOAUGUSTSC" = "0.4285714",
                                  "YEAR" = "2018"
                    ), type = "fixed")

p1$group <- "LandSpAES"
p1_int$group <- "Integrated"



#get raw data for rug
p1_raw <- attr(p1, "rawdata") %>%
  mutate(x = (x*6966.606 + 3340.947)/1000,
         group = "LandSpAES")
p1_int_raw <- attr(p1_int, "rawdata") %>%
  mutate(x = (x*6966.606 + 3340.947)/1000,
         group = "Integrated")



p4 <- rbind(p1, p1_int) %>%
  mutate(x = (x*6966.606 + 3340.947)/1000,
         group = forcats::fct_inorder(as.factor(group))) # orders factor in order it appears


abundint_1km <- ggplot(p4, aes(x = x, y = predicted, colour = group, fill = group)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.3, colour = NA) +
  geom_line() + 
  geom_rug(data = p1_raw, aes(x = x, y = response, colour = group),
           sides = "b") +
  geom_rug(data = p1_int_raw, aes(x = x, y = response, colour = group),
           sides = "t") +
  coord_cartesian(clip = "off") +
  scale_fill_manual(aesthetics = c("fill","colour"),
                    values = my_col,
                    name = "Survey") +
  scale_y_continuous(limits = c(5,1600), expand = c(0,0))+
  labs(x = "AES 1km ('000s)", y = "Predicted Butterfly Abundance") +
  NULL

abundint_1km 

ggsave("report/Butterfly abundance 1km integrated model.png", height = 800, width = 1000, units = "mm", scale = 0.15)

p2_int <- ggpredict(buttabund_int_mod1, "AES3KM[-0.86:9.64, by = 0.5]",
                    condition = c("AES1KM" = "0",
                                  "N_VISITS_MAYTOAUGUSTSC" = "0.4285714",
                                  "YEAR" = "2018"), type = "fixed")


p2 <- ggpredict(buttabund_mod1, "AES3KM[-0.86:5.14, by = 0.5]",
                condition = c("AES1KM" = "0",
                              "ROUND_NUMBER" = "4",
                              "SURVEY_YEAR" = "2018"))

p2$group <- "LandSpAES"
p2_int$group <- "Integrated"


#get raw data for rug
p2_raw <- attr(p2, "rawdata") %>%
  mutate(x = (x*3745.525 + 3245.172)/1000,
         group = "LandSpAES")
p2_int_raw <- attr(p2_int, "rawdata") %>%
  mutate(x = (x*3745.525 + 3245.172)/1000,
         group = "Integrated")



p4 <- rbind(p2, p2_int) %>%
  mutate(x = (x*3745.525 + 3245.172)/1000,
         group = forcats::fct_inorder(as.factor(group))) # orders factor in order it appears


abundint_3km <- ggplot(p4, aes(x = x, y = predicted, colour = group, fill = group)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.3, colour = NA) +
  geom_line() + 
  geom_rug(data = p2_raw, aes(x = x, y = response, colour = group),
           sides = "b") +
  geom_rug(data = p2_int_raw, aes(x = x, y = response, colour = group),
           sides = "t") +
  coord_cartesian(clip = "off") +
  scale_fill_manual(aesthetics = c("fill","colour"),
                    values = my_col,
                    name = "Survey") +
  scale_y_continuous(limits = c(5,1500), expand = c(0,0))+
  labs(x = "AES 3km ('000s)", y = "Predicted Butterfly Abundance") +
  NULL

abundint_3km

ggsave("report/Butterfly abundance 3km integrated model.png", height = 800, width = 1000, units = "mm", scale = 0.15)

abundint_1km + abundint_3km + plot_layout(guides='collect') &
  theme(legend.position='bottom')

ggsave("report/Butterfly abundance combined 1km and 3km integrated model.png", height = 800, width = 1200, units = "mm", scale = 0.15)


p3_int <- ggpredict(buttabund_int_mod1, c("AES1KM[-0.48:5.52, by = 0.5]",
                                         "AES3KM[-0.8,-0.199,1.803]"),
                    condition = c("N_VISITS_MAYTOAUGUSTSC" = "0.4285714",
                                  "YEAR" = "2018"), type = "fixed")

p3 <- ggpredict(buttabund_mod1, c("AES1KM[-0.48:5.02, by = 0.5]",
                                 "AES3KM[-0.8,-0.199,1.803]"),
                condition = c("ROUND_NUMBER" = "4",
                              "SURVEY_YEAR" = "2018"), type = "fixed")

p3_int$AES3KM <- rep(c(-0.8,-0.199,1.803),length.out = nrow(p3_int))
p3$AES3KM <- rep(c(-0.8,-0.199,1.803),length.out = nrow(p3))

p3$group <- "LandSpAES"
p3_int$group <- "Integrated"

p4 <- rbind(p3, p3_int) %>%
  mutate(x = x*6966.606 + 3340.947,
         AES3KM = as.factor(paste("AES 3km:",
                                  round(AES3KM*3745.525 + 3245.172, -1)))) %>%
  mutate(AES3KM = forcats::fct_inorder(AES3KM),
         group = forcats::fct_inorder(as.factor(group)))


ggplot(p4, aes(x = x, y = predicted, colour = group, fill = group)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.3, colour = NA) +
  geom_line() + 
  scale_fill_manual(aesthetics = c("fill","colour"),
                    values = my_col,
                    name = "Survey") +
  facet_wrap(~AES3KM)+
  scale_y_continuous(limits = c(0,3000), expand = c(0,0))+
  labs(x = "AES 1km", y = "Predicted Butterfly Abundance") +
  NULL

ggsave("report/Butterfly abundance 1km x 3km integrated model.png", height = 800, width = 1400, units = "mm", scale = 0.15)

###Diversity



p1 <- ggpredict(buttdiv_mod1, "AES1KM[-0.48:5.02, by = 0.5]",
                condition = c("AES3KM" = "0.01",
                              "ROUND_NUMBER" = "4",
                              "SURVEY_YEAR" = "2018",
                              "Comp.1" = "0",
                              "Comp.2" = "0"))


p1_int <- ggpredict(buttdiv_int_mod1, "AES1KM[-0.48:6.52, by = 0.5]",
                    condition = c("AES3KM" = "0.01",
                                  "N_VISITS_MAYTOAUGUSTSC" = "0.0461",
                                  "YEAR" = "2018"
                    ), type = "fixed")

p1$group <- "LandSpAES"
p1_int$group <- "Integrated"


#get raw data for rug
p1_raw <- attr(p1, "rawdata") %>%
  mutate(x = (x*6966.606 + 3340.947)/1000,
         group = "LandSpAES")
p1_int_raw <- attr(p1_int, "rawdata") %>%
  mutate(x = (x*6966.606 + 3340.947)/1000,
         group = "Integrated")


p4 <- rbind(p1, p1_int) %>%
  mutate(x = (x*6966.606 + 3340.947)/1000,
         group = forcats::fct_inorder(as.factor(group))) # orders factor in order it appears


divint_1km <- ggplot(p4, aes(x = x, y = predicted, colour = group, fill = group)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.3, colour = NA) +
  geom_line() + 
  geom_rug(data = p1_raw, aes(x = x, y = response, colour = group),
           sides = "b") +
  geom_rug(data = p1_int_raw, aes(x = x, y = response, colour = group),
           sides = "t") +
  coord_cartesian(clip = "off") +
  scale_fill_manual(aesthetics = c("fill","colour"),
                    values = my_col,
                    name = "Survey") +
  scale_y_continuous(limits = c(0,14), expand = c(0,0))+
  labs(x = "AES 1km ('000s)", y = "Predicted Butterfly Diversity") +
  NULL

divint_1km

ggsave("report/Butterfly diversity 1km integrated model random slope.png", height = 800, width = 1000, units = "mm", scale = 0.15)

p2_int <- ggpredict(buttdiv_int_mod1, "AES3KM[-0.86:9.64, by = 0.5]",
                    condition = c("AES1KM" = "0",
                                  "N_VISITS_MAYTOAUGUSTSC" = "0.0461",
                                  "YEAR" = "2018"), type = "fixed")


p2 <- ggpredict(buttdiv_mod1, "AES3KM[-0.86:5.14, by = 0.5]",
                condition = c("AES1KM" = "0",
                              "ROUND_NUMBER" = "4",
                              "SURVEY_YEAR" = "2018"))

p2$group <- "LandSpAES"
p2_int$group <- "Integrated"


#get raw data for rug
p2_raw <- attr(p2, "rawdata") %>%
  mutate(x = (x*3745.525 + 3245.172)/1000,
         group = "LandSpAES")
p2_int_raw <- attr(p2_int, "rawdata") %>%
  mutate(x = (x*3745.525 + 3245.172)/1000,
         group = "Integrated")



p4 <- rbind(p2, p2_int) %>%
  mutate(x = (x*3745.525 + 3245.172)/1000,
         group = forcats::fct_inorder(as.factor(group))) # orders factor in order it appears


divint_3km <- ggplot(p4, aes(x = x, y = predicted, colour = group, fill = group)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.3, colour = NA) +
  geom_line() + 
  geom_rug(data = p2_raw, aes(x = x, y = response, colour = group),
           sides = "b") +
  geom_rug(data = p2_int_raw, aes(x = x, y = response, colour = group),
           sides = "t") +
  coord_cartesian(clip = "off") +
  scale_fill_manual(aesthetics = c("fill","colour"),
                    values = my_col,
                    name = "Survey") +
  scale_y_continuous(limits = c(-1,12), expand = c(0,0))+
  labs(x = "AES 3km ('000s)", y = "Predicted Butterfly Diversity") +
  NULL

divint_3km

ggsave("report/Butterfly diversity 3km integrated model random slope.png", height = 800, width = 1000, units = "mm", scale = 0.15)

divint_1km + divint_3km + plot_layout(guides='collect') &
  theme(legend.position='bottom')

ggsave("report/Butterfly diversity combined 1km and 3km integrated model.png", height = 800, width = 1200, units = "mm", scale = 0.15)

p3_int <- ggpredict(buttdiv_int_mod1, c("AES1KM[-0.48:6.52, by = 0.5]",
                                         "AES3KM[-0.8,-0.199,1.803]"),
                    condition = c("N_VISITS_MAYTOAUGUSTSC" = "0.0461",
                                  "YEAR" = "2018"), type = "fixed")

p3 <- ggpredict(buttdiv_mod1, c("AES1KM[-0.48:5.02, by = 0.5]",
                                 "AES3KM[-0.8,-0.199,1.803]"),
                condition = c("ROUND_NUMBER" = "4",
                              "SURVEY_YEAR" = "2018"), type = "fixed")

p3_int$AES3KM <- rep(c(-0.8,-0.199,1.803),length.out = nrow(p3_int))
p3$AES3KM <- rep(c(-0.8,-0.199,1.803),length.out = nrow(p3))

p3$group <- "LandSpAES"
p3_int$group <- "Integrated"

p4 <- rbind(p3, p3_int) %>%
  mutate(x = (x*6966.606 + 3340.947)/1000,
         AES3KM = as.factor(paste("AES 3km:",
                                  round(AES3KM*3745.525 + 3245.172, -1)))) %>%
  mutate(AES3KM = forcats::fct_inorder(AES3KM),
         group = forcats::fct_inorder(as.factor(group)))


ggplot(p4, aes(x = x, y = predicted, colour = group, fill = group)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.3, colour = NA) +
  geom_line() + 
  scale_fill_manual(aesthetics = c("fill","colour"),
                    values = my_col,
                    name = "Survey") +
  facet_wrap(~AES3KM)+
  scale_y_continuous(limits = c(0,16), expand = c(0,0))+
  labs(x = "AES 1km ('000s)", y = "Predicted Butterfly Diversity") +
  NULL

ggsave("report/Butterfly diversity 1km x 3km integrated model random slope.png", height = 800, width = 1400, units = "mm", scale = 0.15)
