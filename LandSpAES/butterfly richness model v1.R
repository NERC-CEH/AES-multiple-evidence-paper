###Basic richness models

homedir <- getwd()

analysis_group <- "Headline analyses"
analysis_name <- "Butterfly species richness"
plot_name <- "Total butterfly species richness"


# folder setup for saving
dir <- config::get()
file.dir <- paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/")
dir.create(file.dir, recursive = TRUE)


# data collation and aggregation
source("collate butterfly data.R")

# functions
source("functions_script.R")

library(reshape)
library(reshape2)
library(vegan)
library(lme4)


##Now aggregate data for models



buttcount2$RICHNESS_ID <- buttcount2$BUTTERFLY_ID

#remove "other butterflies" id

buttcount2$RICHNESS_ID[buttcount2$RICHNESS_ID == 65] <- NA


length(buttcount2$RICHNESS_ID[is.na(buttcount2$RICHNESS_ID)])
# 1898/19093 no butterflies


##aggregate data

buttrich <- aggregate(RICHNESS_ID ~ NCA + SURVEY_SQUARE + SURVEY_YEAR + AES1KM + AES3KM, data = buttcount2, FUN = function(x) length(unique(x[!is.na(x)])), na.action = na.pass)
#should be 198 rows long = (54 squares*3 years) + 36 squares for year 1

buttrounds <- aggregate(ROUND_NUMBER ~ SURVEY_SQUARE + SURVEY_YEAR, data = buttcount2, FUN = function(x) length(unique(x)))

buttrich <- merge(buttrich, buttrounds, by.x = c("SURVEY_SQUARE", "SURVEY_YEAR"), by.y = c("SURVEY_SQUARE", "SURVEY_YEAR"))

#sunshine
buttsun <- aggregate(PERC_SUN ~ SURVEY_SQUARE + SURVEY_YEAR, data = buttvariables, FUN = function(x) mean(x, na.rm = TRUE))

buttrich <- merge(buttrich, buttsun, by.x = c("SURVEY_SQUARE", "SURVEY_YEAR"), by.y = c("SURVEY_SQUARE", "SURVEY_YEAR"))

#temperature
butttemp <- aggregate(SURVEY_TEMP_SHADE ~ SURVEY_SQUARE + SURVEY_YEAR, data = buttvisit, FUN = function(x) mean(x, na.rm = TRUE))

buttrich <- merge(buttrich, butttemp, by.x = c("SURVEY_SQUARE", "SURVEY_YEAR"), by.y = c("SURVEY_SQUARE", "SURVEY_YEAR"))


#' Plot histogram of data
#' 
fit_distrib(buttrich, 'RICHNESS_ID')

#not much evidence of over or underdispersion, use Poisson

library(lme4)

#scale predictors
buttrich$AES1KM <- scale(buttrich$AES1KM)
buttrich$AES3KM <- scale(buttrich$AES3KM)

buttrich$SURVEY_TEMP_SHADE <- scale(buttrich$SURVEY_TEMP_SHADE)
buttrich$PERC_SUN <- scale(buttrich$PERC_SUN)

#treat survey year as factor
buttrich$SURVEY_YEAR <- factor(buttrich$SURVEY_YEAR)


#go for Poisson


mod1 <- glmer(RICHNESS_ID ~ AES1KM*AES3KM + ROUND_NUMBER + PERC_SUN + SURVEY_TEMP_SHADE + SURVEY_YEAR + (1|NCA/SURVEY_SQUARE), data = buttrich, family = "poisson")

summary(mod1)


# Try year as a random effect 
mod1_yrre <- glmer.nb(RICHNESS_ID ~ AES1KM*AES3KM + ROUND_NUMBER + PERC_SUN + SURVEY_TEMP_SHADE + 
                        (1|NCA/SURVEY_SQUARE) + (1|SURVEY_YEAR), 
                      data = buttrich)
summary(mod1_yrre)



#outlier has big effect?

#outlier defined as AES1KM scores above 23061.94 - this removes square SP4469 in years 2017, 2018 and 2019 but not 2021

buttrich2 <- buttrich[buttrich$SURVEY_SQUARE != 'SP4469',]

mod2 <- glmer(RICHNESS_ID ~ AES1KM*AES3KM + ROUND_NUMBER + PERC_SUN + SURVEY_TEMP_SHADE + SURVEY_YEAR + (1|NCA/SURVEY_SQUARE), data = buttrich2, family = "poisson")

summary(mod2)

#' Removing outlier does not substantially affect conclusions or validation, probably keep it in


#survey square effect?

boxplot(resid(mod1) ~ buttrich$SURVEY_SQUARE)
#yes, but square level random effect won't fit

# ## Plots
#' 
#' # change factors for boxplot
buttrich$NCA <- fct_recode(buttrich$NCA, "High\nWeald" = "High Weald", "South Suffolk\nand North\nEssex Clayland" = "South Suffolk and North Essex Clayland", 
                           "Dunsmore\nand Feldon" = "Dunsmore and Feldon", "The\nFens" = "The Fens", 
                           "Dartmoor" = "Dartmoor", "Yorkshire\nDales" = "Yorkshire Dales")
buttrich$NCA <- factor(buttrich$NCA,
                       levels = c("High\nWeald", "South Suffolk\nand North\nEssex Clayland", 
                                  "Dunsmore\nand Feldon", "The\nFens", 
                                  "Dartmoor", "Yorkshire\nDales"))

pred.data1 = expand.grid(AES1KM=seq(min(buttrich$AES1KM), max(buttrich$AES1KM), 0.2), 
                         AES3KM = mean(buttrich$AES3KM), 
                         PERC_SUN = mean(buttrich$PERC_SUN),
                         SURVEY_TEMP_SHADE = mean(buttrich$SURVEY_TEMP_SHADE),
                         NCA = NA, 
                         ROUND_NUMBER = mean(buttrich$ROUND_NUMBER), 
                         SURVEY_YEAR = factor(2018, levels = levels(buttrich$SURVEY_YEAR)))

pred.data2 = expand.grid(AES1KM=mean(buttrich$AES1KM), 
                         AES3KM = seq(min(buttrich$AES3KM), max(buttrich$AES3KM), 0.2),
                         PERC_SUN = mean(buttrich$PERC_SUN),
                         SURVEY_TEMP_SHADE = mean(buttrich$SURVEY_TEMP_SHADE), 
                         NCA = NA, 
                         ROUND_NUMBER = mean(buttrich$ROUND_NUMBER), 
                         SURVEY_YEAR = factor(2017, levels = levels(buttrich$SURVEY_YEAR)))


## predict
pred.data1 <- pred_conf(mod = mod1, pred_data = pred.data1, col_id = 'RICHNESS_ID', re_name = 'NCA')
pred.data2 <- pred_conf(mod = mod1, pred_data = pred.data2, col_id = 'RICHNESS_ID', re_name = 'NCA')


## plotting
plot_results(predict_data = pred.data1,
             xvars = 'AES1KM', 
             x_var_names = 'Local level',
             yvar = 'RICHNESS_ID', 
             boxplot_x = 'NCA',
             original_data = buttrich,
             ylim = c(1,30),
             unscale = TRUE,
             write = TRUE,
             plot_dir = paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/"))

plot_results(predict_data = pred.data2,
             xvars = 'AES3KM', 
             x_var_names = 'Landscape level',
             yvar = 'RICHNESS_ID', 
             boxplot_x = 'NCA', 
             original_data = buttrich,
             ylim = c(1,30),
             unscale = TRUE,
             write = TRUE,
             plot_dir = paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/"))

#write out model results

results <- capture.output(summary(mod1))

cat(paste0(analysis_name, " model results"), results, 
    file=paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/", analysis_name," model results.txt"), sep="\n", append=TRUE)

results_outlier <- capture.output(summary(mod2))

cat(paste0(analysis_name, " model results - outlier removed"), results_outlier, 
    file=paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/", analysis_name," model results outlier removed.txt"), sep="\n", append=TRUE)



model_table(model = mod1,
            NCA = 'All NCA',
            outlier = 'All data')

model_table(model = mod2,
            NCA = 'All NCA',
            outlier = 'Outlier removed')


#######################################
#####     Lowland only models     #####
#######################################

## create lowland directory
dir.create(paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/lowland_only/"))

buttrich_lowland <- buttrich[!buttrich$NCA %in% c("Dartmoor", "Yorkshire\nDales"),]


#' Plot histogram of data
#' 
fit_distrib(buttrich_lowland, 'RICHNESS_ID')

#go for Poisson


mod1_lowland <- glmer(RICHNESS_ID ~ AES1KM*AES3KM + ROUND_NUMBER + PERC_SUN + SURVEY_TEMP_SHADE + SURVEY_YEAR + (1|NCA/SURVEY_SQUARE), data = buttrich_lowland, family = "poisson")

summary(mod1_lowland)
# model is singular - try without NCA

mod1_lowland2 <- glmer(RICHNESS_ID ~ AES1KM*AES3KM + ROUND_NUMBER + PERC_SUN + SURVEY_TEMP_SHADE + SURVEY_YEAR + (1|SURVEY_SQUARE), data = buttrich_lowland, family = "poisson")

summary(mod1_lowland2)
# failed to converge - try with different optimiser

# with NCA 
mod1_lowland3 <- glmer(RICHNESS_ID ~ AES1KM*AES3KM + ROUND_NUMBER + PERC_SUN + SURVEY_TEMP_SHADE + SURVEY_YEAR + 
                         (1|NCA/SURVEY_SQUARE), 
                       data = buttrich_lowland, 
                       family = "poisson", control=glmerControl(optimizer="bobyqa"))
# still singular - try without NCA

mod1_lowland4 <- glmer(RICHNESS_ID ~ AES1KM*AES3KM + ROUND_NUMBER + PERC_SUN + SURVEY_TEMP_SHADE + SURVEY_YEAR + 
                         (1|SURVEY_SQUARE), 
                       data = buttrich_lowland, 
                       family = "poisson", control=glmerControl(optimizer="bobyqa"))
summary(mod1_lowland4)
# no effect on AES variables compared with first model. Temp becomes significant.
# Use this for predictions...!


#########################
####    without outlier

buttrich2_lowland <- buttrich_lowland[buttrich_lowland$SURVEY_SQUARE != 'SP4469',]

mod2_lowland <- glmer(RICHNESS_ID ~ AES1KM*AES3KM + ROUND_NUMBER + PERC_SUN + SURVEY_TEMP_SHADE + SURVEY_YEAR + (1|NCA/SURVEY_SQUARE), data = buttrich2_lowland, family = "poisson")

summary(mod2_lowland)
# model is singular without outlier - try without NCA

mod2_lowland2 <- glmer(RICHNESS_ID ~ AES1KM*AES3KM + ROUND_NUMBER + PERC_SUN + SURVEY_TEMP_SHADE + SURVEY_YEAR + (1|SURVEY_SQUARE), data = buttrich2_lowland, family = "poisson")

summary(mod2_lowland2)
# model is singular without outlier - try same model with bobyqa optimiser, like model with outlier


mod2_lowland3 <- glmer(RICHNESS_ID ~ AES1KM*AES3KM + ROUND_NUMBER + PERC_SUN + SURVEY_TEMP_SHADE + SURVEY_YEAR + 
                         (1|SURVEY_SQUARE), 
                       data = buttrich2_lowland, family = "poisson", control=glmerControl(optimizer="bobyqa"))

summary(mod2_lowland3)
# all good, very little effect of outlier.

#survey square effect?

boxplot(resid(mod1_lowland4) ~ buttrich_lowland$SURVEY_SQUARE)
# some squares have negative residuals...

# ## Plots
#' 
#' # change factors for boxplot
buttrich_lowland$NCA <- factor(buttrich_lowland$NCA,
                               levels = c("High\nWeald", "South Suffolk\nand North\nEssex Clayland", 
                                          "Dunsmore\nand Feldon", "The\nFens", 
                                          "Dartmoor", "Yorkshire\nDales"))

pred.data1_lowland = expand.grid(AES1KM=seq(min(buttrich_lowland$AES1KM), max(buttrich_lowland$AES1KM), 0.2), 
                                 AES3KM = mean(buttrich_lowland$AES3KM), 
                                 PERC_SUN = mean(buttrich_lowland$PERC_SUN),
                                 SURVEY_TEMP_SHADE = mean(buttrich_lowland$SURVEY_TEMP_SHADE),
                                 NCA = NA, 
                                 ROUND_NUMBER = mean(buttrich_lowland$ROUND_NUMBER), 
                                 SURVEY_YEAR = factor(2018, levels = levels(buttrich_lowland$SURVEY_YEAR)))

pred.data2_lowland = expand.grid(AES1KM=mean(buttrich_lowland$AES1KM), 
                                 AES3KM = seq(min(buttrich_lowland$AES3KM), max(buttrich_lowland$AES3KM), 0.2),
                                 PERC_SUN = mean(buttrich_lowland$PERC_SUN),
                                 SURVEY_TEMP_SHADE = mean(buttrich_lowland$SURVEY_TEMP_SHADE), 
                                 NCA = NA, 
                                 ROUND_NUMBER = mean(buttrich_lowland$ROUND_NUMBER), 
                                 SURVEY_YEAR = factor(2017, levels = levels(buttrich_lowland$SURVEY_YEAR)))


## predict
pred.data1_lowland <- pred_conf(mod = mod1_lowland4, pred_data = pred.data1_lowland, col_id = 'RICHNESS_ID', re_name = 'SURVEY_SQUARE')
pred.data2_lowland <- pred_conf(mod = mod1_lowland4, pred_data = pred.data2_lowland, col_id = 'RICHNESS_ID', re_name = 'SURVEY_SQUARE')


## plotting
plot_results(predict_data = pred.data1_lowland,
             xvars = 'AES1KM', 
             x_var_names = 'Local level',
             yvar = 'RICHNESS_ID', 
             boxplot_x = 'NCA',
             original_data = buttrich_lowland,
             scaled_data = buttrich,
             ylim = c(1,30),
             unscale = TRUE,
             write = TRUE,
             lowland_only = TRUE,
             plot_dir = paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/lowland_only"))

plot_results(predict_data = pred.data2_lowland,
             xvars = 'AES3KM', 
             x_var_names = 'Landscape level',
             yvar = 'RICHNESS_ID', 
             boxplot_x = 'NCA', 
             original_data = buttrich_lowland,
             scaled_data = buttrich,
             ylim = c(1,30),
             unscale = TRUE,
             write = TRUE,
             lowland_only = TRUE,
             plot_dir = paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/lowland_only"))

#write out model results

results <- capture.output(summary(mod1_lowland4))

cat(paste0(analysis_name, " model results"), results, 
    file=paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/lowland_only/", analysis_name," model results.txt"), sep="\n", append=FALSE)

results_outlier <- capture.output(summary(mod2_lowland3))

cat(paste0(analysis_name, " model results - outlier removed"), results_outlier, 
    file=paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/lowland_only/", analysis_name," model results outlier removed.txt"), sep="\n", append=FALSE)


model_table(model = mod1_lowland4,
            NCA = 'Lowland only',
            outlier = 'All data')

model_table(model = mod2_lowland3,
            NCA = 'Lowland only',
            outlier = 'Outlier removed')


#######################################
#####        Habitat models       #####
#######################################

source("construct habitat variables.R")
source("collate botanical data.R")
source("collate floral data - transects.R")

#combine covariate datasets

covariates <- merge(hab_vars, floral_vars, by = c("SURVEY_SQUARE", "SURVEY_YEAR"))

covariates <- merge(covariates, botanical_vars, by = "SURVEY_SQUARE")

buttrich_habs <- merge(buttrich, covariates, by.x = c("SURVEY_SQUARE", "SURVEY_YEAR"), by.y = c("SURVEY_SQUARE", "SURVEY_YEAR"))

par(mfrow=c(2,2))
for (i in 11:22){
  plot(RICHNESS_ID ~ buttrich_habs[,i], data = buttrich_habs, xlab = names(buttrich_habs)[i])
}

mod3 <- glmer(RICHNESS_ID ~ AES1KM*AES3KM + ROUND_NUMBER + SURVEY_YEAR + PERC_SUN + SURVEY_TEMP_SHADE + Habitat.diversity + BL_wood + Woody.linear + Water.linear + (1|NCA/SURVEY_SQUARE), data = buttrich_habs, family = "poisson", control=glmerControl(optimizer="bobyqa"))

summary(mod3)

#remove square level
mod3 <- glmer(RICHNESS_ID ~ AES1KM*AES3KM + ROUND_NUMBER + SURVEY_YEAR + PERC_SUN + SURVEY_TEMP_SHADE + Habitat.diversity + BL_wood + Woody.linear + Water.linear + (1|NCA), data = buttrich_habs, family = "poisson", control=glmerControl(optimizer="bobyqa"))

summary(mod3)

results_habitats <- capture.output(summary(mod3))


cat(paste0(analysis_name, " model results habitat analyses"), results_habitats, file=paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/", analysis_name," model results habitat analyses.txt"), sep="\n", append=FALSE)


#' plant models
#' 

mod4 <- glmer(RICHNESS_ID ~ AES1KM*AES3KM + ROUND_NUMBER + SURVEY_YEAR + PERC_SUN + SURVEY_TEMP_SHADE + BOT_DIV + PERC_GRAM + FLORAL_RESOURCES + (1|NCA/SURVEY_SQUARE), data = buttrich_habs, family = "poisson", control=glmerControl(optimizer="bobyqa"))

summary(mod4)

#remove survey square
mod4 <- glmer(RICHNESS_ID ~ AES1KM*AES3KM + ROUND_NUMBER + SURVEY_YEAR + PERC_SUN + SURVEY_TEMP_SHADE + BOT_DIV + PERC_GRAM + FLORAL_RESOURCES + (1|NCA), data = buttrich_habs, family = "poisson", control=glmerControl(optimizer="bobyqa"))

summary(mod4)


#habitat diversity plot

pred.data1 = expand.grid(Habitat.diversity=seq(min(buttrich_habs$Habitat.diversity), max(buttrich_habs$Habitat.diversity), 0.05), 
                         AES3KM = mean(buttrich_habs$AES3KM), 
                         AES1KM = mean(buttrich_habs$AES1KM), 
                         PERC_SUN = mean(buttrich_habs$PERC_SUN),
                         SURVEY_TEMP_SHADE = mean(buttrich_habs$SURVEY_TEMP_SHADE),
                         BL_wood = mean(buttrich_habs$BL_wood),
                         Woody.linear = mean(buttrich_habs$Woody.linear),
                         Water.linear = mean(buttrich_habs$Water.linear),
                         NCA = NA, 
                         ROUND_NUMBER = mean(buttrich_habs$ROUND_NUMBER), 
                         SURVEY_YEAR = factor(2018, levels = levels(buttrich_habs$SURVEY_YEAR)))

pred.data1 <- pred_conf(mod = mod3, pred_data = pred.data1, col_id = 'RICHNESS_ID', re_name = 'NCA')


plot(exp(pred.data1$RICHNESS_ID) ~ pred.data1$Habitat.diversity)

plot_covariates(pred.data1, "Habitat.diversity", "RICHNESS_ID", original_data = buttrich_habs, unscale = FALSE, x_var_names = "Habitat diversity", write = TRUE, lowland_only = FALSE, plot_dir = paste0(dir$directories$outpath, analysis_group, "/", analysis_name,"/"))


#botanical richness plot

pred.data2 = expand.grid(BOT_DIV=seq(min(buttrich_habs$BOT_DIV), max(buttrich_habs$BOT_DIV), 0.02), 
                         AES3KM = mean(buttrich_habs$AES3KM), 
                         AES1KM = mean(buttrich_habs$AES1KM), 
                         PERC_SUN = mean(buttrich_habs$PERC_SUN),
                         SURVEY_TEMP_SHADE = mean(buttrich_habs$SURVEY_TEMP_SHADE),
                         PERC_GRAM = mean(buttrich_habs$PERC_GRAM),
                         FLORAL_RESOURCES = mean(buttrich_habs$FLORAL_RESOURCES),
                         NCA = NA, 
                         ROUND_NUMBER = mean(buttrich_habs$ROUND_NUMBER), 
                         SURVEY_YEAR = factor(2018, levels = levels(buttrich_habs$SURVEY_YEAR)))

pred.data2 <- pred_conf(mod = mod4, pred_data = pred.data2, col_id = 'RICHNESS_ID', re_name = 'NCA')


plot_covariates(pred.data2, "BOT_DIV", "RICHNESS_ID", original_data = buttrich_habs, unscale = FALSE, x_var_names = "Botanical diversity", write = TRUE, lowland_only = FALSE, plot_dir = paste0(dir$directories$outpath, analysis_group, "/", analysis_name,"/"))



results_plants <- capture.output(summary(mod4))


cat(paste0(analysis_name, " model results plant analyses"), results_plants, 
    file=paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/", analysis_name," model results plant analyses.txt"), sep="\n", append=FALSE)

model_table(model = mod3,
            NCA = 'All NCA',
            outlier = 'All data',
            saving_name = paste(analysis_name, "habitat model"))

model_table(model = mod4,
            NCA = 'All NCA',
            outlier = 'All data',
            saving_name = paste(analysis_name, "plant model"))



####################################################
#####        Lowland only habitat models       #####
####################################################

buttrich_habslowland <- buttrich_habs[!buttrich_habs$NCA %in% c("Dartmoor", "Yorkshire\nDales"),]


par(mfrow=c(2,2))
for (i in 11:22){
  plot(RICHNESS_ID ~ buttrich_habslowland[,i], data = buttrich_habslowland, xlab = names(buttrich_habslowland)[i])
}

mod5 <- glmer(RICHNESS_ID ~ AES1KM*AES3KM + ROUND_NUMBER + SURVEY_YEAR + PERC_SUN + SURVEY_TEMP_SHADE + Habitat.diversity + Water.linear + (1|NCA/SURVEY_SQUARE), data = buttrich_habslowland, family = "poisson", control=glmerControl(optimizer="bobyqa"))

summary(mod5)

#remove square level
mod5 <- glmer(RICHNESS_ID ~ AES1KM*AES3KM + ROUND_NUMBER + SURVEY_YEAR + PERC_SUN + SURVEY_TEMP_SHADE + Habitat.diversity + Water.linear + (1|NCA), data = buttrich_habslowland, family = "poisson", control=glmerControl(optimizer="bobyqa"))

summary(mod5)


results_habitats <- capture.output(summary(mod5))


cat(paste0(analysis_name, " model results habitat analyses"), results_habitats, file=paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/", "lowland_only/Lowland only ", analysis_name," model results habitat analyses.txt"), sep="\n", append=FALSE)


#' plant models
#' 

mod6 <- glmer(RICHNESS_ID ~ AES1KM*AES3KM + ROUND_NUMBER + SURVEY_YEAR + PERC_SUN + SURVEY_TEMP_SHADE + BOT_DIV + PERC_GRAM + FLORAL_RESOURCES + (1|NCA), data = buttrich_habslowland, family = "poisson", control=glmerControl(optimizer="bobyqa"))

summary(mod6)


results_plants <- capture.output(summary(mod6))


cat(paste0(analysis_name, " model results plant analyses"), results_plants, 
    file=paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/", "lowland_only/Lowland only ", analysis_name," model results plant analyses.txt"), sep="\n", append=FALSE)

model_table(model = mod5,
            NCA = 'Lowland only',
            outlier = 'All data',
            saving_name = paste(analysis_name, "habitat model"),
            save_dir = paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/lowland_only/"))

model_table(model = mod6,
            NCA = 'Lowland only',
            outlier = 'All data',
            saving_name = paste(analysis_name, "plant model"),
            save_dir = paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/lowland_only/"))



####################################################
#####        No square level random term       #####
####################################################


mod1_nosq <- update(mod1, . ~ . -(1|NCA/SURVEY_SQUARE) -(1|SURVEY_SQUARE) + (1|NCA))
mod2_nosq <- update(mod2, . ~ . -(1|NCA/SURVEY_SQUARE) -(1|SURVEY_SQUARE) + (1|NCA))
mod1_lowland_nosq <- update(mod1_lowland, . ~ . -(1|NCA/SURVEY_SQUARE) -(1|SURVEY_SQUARE) + (1|NCA))
mod2_lowland_nosq <- update(mod2_lowland, . ~ . -(1|NCA/SURVEY_SQUARE) -(1|SURVEY_SQUARE) + (1|NCA))
mod3_nosq <- update(mod3, . ~ . -(1|NCA/SURVEY_SQUARE) -(1|SURVEY_SQUARE) + (1|NCA))
mod4_nosq <- update(mod4, . ~ . -(1|NCA/SURVEY_SQUARE) -(1|SURVEY_SQUARE) + (1|NCA))
mod5_nosq <- update(mod5, . ~ . -(1|NCA/SURVEY_SQUARE) -(1|SURVEY_SQUARE) + (1|NCA))
mod6_nosq <- update(mod6, . ~ . -(1|NCA/SURVEY_SQUARE) -(1|SURVEY_SQUARE) + (1|NCA))

results_mod1_nosq <- capture.output(summary(mod1_nosq))
results_mod2_nosq <- capture.output(summary(mod2_nosq))
results_mod1_lowland_nosq <- capture.output(summary(mod1_lowland_nosq))
results_mod2_lowland_nosq <- capture.output(summary(mod2_lowland_nosq))
results_mod3_nosq <- capture.output(summary(mod3_nosq))
results_mod4_nosq <- capture.output(summary(mod4_nosq))
results_mod5_nosq <- capture.output(summary(mod5_nosq))
results_mod6_nosq <- capture.output(summary(mod6_nosq))


## create no square re directory
dir.create(paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/no_sq/"))

cat(paste0(analysis_name, " model results no square"), results_mod1_nosq, 
    file=paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/no_sq/No square ", analysis_name," model results.txt"), sep="\n", append=FALSE)
cat(paste0(analysis_name, " model results no outlier no square"), results_mod2_nosq, 
    file=paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/no_sq/No square ", analysis_name," model results no outlier.txt"), sep="\n", append=FALSE)
cat(paste0(analysis_name, " model results lowland only no square"), results_mod1_lowland_nosq, 
    file=paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/no_sq/No square ", analysis_name," model results lowland only.txt"), sep="\n", append=FALSE)
cat(paste0(analysis_name, " model results lowland only no outlier no square"), results_mod2_lowland_nosq, 
    file=paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/no_sq/No square ", analysis_name," model results lowland only no outlier.txt"), sep="\n", append=FALSE)
cat(paste0(analysis_name, " model results habitat analyses no square"), results_mod3_nosq, 
    file=paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/no_sq/No square ", analysis_name," model results habitat analyses.txt"), sep="\n", append=FALSE)
cat(paste0(analysis_name, " model results plant analyses no square"), results_mod4_nosq, 
    file=paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/no_sq/No square ", analysis_name," model results plant analyses.txt"), sep="\n", append=FALSE)
cat(paste0(analysis_name, " model results lowland only habitat analyses no square"), results_mod5_nosq, 
    file=paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/no_sq/No square ", analysis_name," model results lowland only habitat analyses.txt"), sep="\n", append=FALSE)
cat(paste0(analysis_name, " model results lowland only plant analyses no square"), results_mod6_nosq, 
    file=paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/no_sq/No square ", analysis_name," model results lowland only plant analyses.txt"), sep="\n", append=FALSE)

model_table(model = mod1_nosq,
            NCA = 'All NCA',
            outlier = 'All data',
            saving_name = paste("No square", analysis_name, "model"),
            save_dir = paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/no_sq/"))

model_table(model = mod2_nosq,
            NCA = 'All NCA',
            outlier = 'No outlier',
            saving_name = paste("No square", analysis_name, "model"),
            save_dir = paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/no_sq/"))

model_table(model = mod1_lowland_nosq,
            NCA = 'Lowland only',
            outlier = 'All data',
            saving_name = paste("No square", analysis_name, "model"),
            save_dir = paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/no_sq/"))

model_table(model = mod2_lowland_nosq,
            NCA = 'Lowland only',
            outlier = 'No outlier',
            saving_name = paste("No square", analysis_name, "model"),
            save_dir = paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/no_sq/"))

model_table(model = mod3_nosq,
            NCA = 'All NCA',
            outlier = 'All data',
            saving_name = paste("No square", analysis_name, "habitat model"),
            save_dir = paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/no_sq/"))

model_table(model = mod4_nosq,
            NCA = 'All NCA',
            outlier = 'All data',
            saving_name = paste("No square", analysis_name, "plant model"),
            save_dir = paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/no_sq/"))

model_table(model = mod5_nosq,
            NCA = 'Lowland only',
            outlier = 'All data',
            saving_name = paste("No square", analysis_name, "lowland only habitat model"),
            save_dir = paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/no_sq/"))

model_table(model = mod6_nosq,
            NCA = 'Lowland only',
            outlier = 'All data',
            saving_name = paste("No square", analysis_name, "lowland only plant model"),
            save_dir = paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/no_sq/"))



