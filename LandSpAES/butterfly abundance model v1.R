###Basic abundance models


analysis_group <- "Headline analyses"
analysis_name <- "Butterfly abundance"
plot_name <- "Total count of butterflies"


# folder setup for saving
dir <- config::get()
file.dir <- paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/")
dir.create(file.dir, recursive = TRUE)


# source data
source("collate butterfly data.R")

# source function
source("functions_script.R")


library(lme4)
library(effects)
library(tidyverse)


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


#' Look at distribution of abundance variable

fit_distrib(buttabund, 'BUTTERFLY_COUNT')

# fits negative binomial fairly well

#' Rescale predictors and fit model
#' 
buttabund$AES1KM <- scale(buttabund$AES1KM)
buttabund$AES3KM <- scale(buttabund$AES3KM)

buttabund$SURVEY_TEMP_SHADE <- scale(buttabund$SURVEY_TEMP_SHADE)
buttabund$PERC_SUN <- scale(buttabund$PERC_SUN)

#treat survye year as factor
buttabund$SURVEY_YEAR <- factor(buttabund$SURVEY_YEAR)


# models
mod1 <- glmer.nb(BUTTERFLY_COUNT ~ AES1KM*AES3KM + ROUND_NUMBER + PERC_SUN + SURVEY_TEMP_SHADE + SURVEY_YEAR + (1|NCA/SURVEY_SQUARE), data = buttabund)

summary(mod1)

# Try year as a random effect 
mod1_yrre <- glmer.nb(BUTTERFLY_COUNT ~ AES1KM*AES3KM + ROUND_NUMBER + PERC_SUN + SURVEY_TEMP_SHADE + 
                        (1|NCA/SURVEY_SQUARE) + (1|SURVEY_YEAR), control=glmerControl(optimizer="bobyqa"), 
                      data = buttabund)
summary(mod1_yrre)

#remove outlier
buttabund2 <- buttabund[buttabund$SURVEY_SQUARE != 'SP4469',]

mod2 <- glmer.nb(BUTTERFLY_COUNT ~ AES1KM*AES3KM + ROUND_NUMBER + PERC_SUN + SURVEY_TEMP_SHADE + SURVEY_YEAR + (1|NCA/SURVEY_SQUARE), data = buttabund2)

summary(mod2)


## check model
boxplot(resid(mod1) ~ buttabund$SURVEY_SQUARE)
#no pattern

#Can use mod1

# ## Plots
#' 
#' 
#' # change factors for boxplot
buttabund$NCA <- fct_recode(buttabund$NCA, "High\nWeald" = "High Weald", "South Suffolk\nand North\nEssex Clayland" = "South Suffolk and North Essex Clayland", 
                            "Dunsmore\nand Feldon" = "Dunsmore and Feldon", "The\nFens" = "The Fens", 
                            "Yorkshire\nDales" = "Yorkshire Dales")
buttabund$NCA <- factor(buttabund$NCA,
                        levels = c("High\nWeald", "South Suffolk\nand North\nEssex Clayland", 
                                   "Dunsmore\nand Feldon", "The\nFens", 
                                   "Dartmoor", "Yorkshire\nDales"))

# 1km
pred.data1 = expand.grid(AES1KM = seq(min(buttabund$AES1KM), max(buttabund$AES1KM), 0.2), 
                         AES3KM = mean(buttabund$AES3KM),
                         PERC_SUN = mean(buttabund$PERC_SUN),
                         SURVEY_TEMP_SHADE = mean(buttabund$SURVEY_TEMP_SHADE),
                         NCA = NA, 
                         ROUND_NUMBER = mean(buttabund$ROUND_NUMBER), 
                         SURVEY_YEAR = factor(2018, levels = levels(buttabund$SURVEY_YEAR)))

# 3km
pred.data2 = expand.grid(AES1KM=mean(buttabund$AES1KM), 
                         AES3KM = seq(min(buttabund$AES3KM), max(buttabund$AES3KM), 0.2), 
                         PERC_SUN = mean(buttabund$PERC_SUN),
                         SURVEY_TEMP_SHADE = mean(buttabund$SURVEY_TEMP_SHADE),
                         NCA = NA, 
                         ROUND_NUMBER = mean(buttabund$ROUND_NUMBER), 
                         SURVEY_YEAR = factor(2018, levels = levels(buttabund$SURVEY_YEAR)))

# interaction?


source("comparison_table.R")

comparison_table(model = mod1, pred.data = pred.data2, data = buttabund, term = "AES3KM", response = "BUTTERFLY_COUNT", low_pred = 250, high_pred = 10000)

## predict
pred.data1 <- pred_conf(mod = mod1, pred_data = pred.data1, col_id = 'BUTTERFLY_COUNT')
pred.data2 <- pred_conf(mod = mod1, pred_data = pred.data2, col_id = 'BUTTERFLY_COUNT')


## plotting (includes boxplot)
plot_results(predict_data = pred.data1,
             xvars = 'AES1KM', 
             x_var_names = 'Local level',
             yvar = 'BUTTERFLY_COUNT', 
             boxplot_x = 'NCA',
             original_data = buttabund,
             unscale = TRUE,
             write = TRUE,
             plot_dir = paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/"))

plot_results(predict_data = pred.data2,
             xvars = 'AES3KM', 
             x_var_names = 'Landscape level',
             yvar = 'BUTTERFLY_COUNT', 
             boxplot_x = 'NCA', 
             original_data = buttabund,
             unscale = TRUE,
             write = TRUE,
             plot_dir = paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/"))


#write out model results

results <- capture.output(summary(mod1))

cat(paste0(analysis_name, " model results"), results, file=paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/", analysis_name," model results.txt"), sep="\n", append=TRUE)


results_outlier <- capture.output(summary(mod2))

cat(paste0(analysis_name, " model results - outlier removed"), results_outlier, file=paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/", analysis_name," model results outlier removed.txt"), sep="\n", append=TRUE)



model_table(model = mod1,
            NCA = 'All NCA',
            outlier = 'All data')

model_table(model = mod2,
            NCA = 'All NCA',
            outlier = 'Outlier removed')


#######################################
#####     Lowland only models     #####
#######################################

buttabund_lowland <- buttabund[!buttabund$NCA %in% c("Dartmoor", "Yorkshire\nDales"),]

## create lowland directory
dir.create(paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/lowland_only/"))

#' Plot histogram of data
#' 
fit_distrib(buttabund_lowland, 'BUTTERFLY_COUNT')
# overdispersed - negative binomial fits best
# go for negbinom


# first model
mod1_lowland <- glmer.nb(BUTTERFLY_COUNT ~ AES1KM*AES3KM + ROUND_NUMBER + PERC_SUN + SURVEY_TEMP_SHADE + SURVEY_YEAR + (1|NCA/SURVEY_SQUARE), 
                         data = buttabund_lowland)
summary(mod1_lowland)
## is singular - maybe because of weird upland square/NCA combination? 

# try removing NCA...
mod1_lowland2 <- glmer.nb(BUTTERFLY_COUNT ~ AES1KM*AES3KM + ROUND_NUMBER + PERC_SUN + SURVEY_TEMP_SHADE + SURVEY_YEAR + (1|SURVEY_SQUARE), 
                          data = buttabund_lowland)
summary(mod1_lowland2)
# results identical


###remove outlier
buttabund_lowland2 <- buttabund_lowland[buttabund_lowland$SURVEY_SQUARE != 'SP4469',]

mod2_lowland <- glmer.nb(BUTTERFLY_COUNT ~ AES1KM*AES3KM + ROUND_NUMBER + PERC_SUN + SURVEY_TEMP_SHADE + SURVEY_YEAR + (1|NCA/SURVEY_SQUARE), 
                         data = buttabund_lowland2)
summary(mod2_lowland)
# is singular

# try removing NCA

mod2_lowland2 <- glmer.nb(BUTTERFLY_COUNT ~ AES1KM*AES3KM + ROUND_NUMBER + PERC_SUN + SURVEY_TEMP_SHADE + SURVEY_YEAR + (1|SURVEY_SQUARE), 
                          data = buttabund_lowland2)
summary(mod2_lowland2) # convergence warning

# try with bobyqa
mod2_lowland3 <- glmer.nb(BUTTERFLY_COUNT ~ AES1KM*AES3KM + ROUND_NUMBER + PERC_SUN + SURVEY_TEMP_SHADE + SURVEY_YEAR + 
                            (1|SURVEY_SQUARE), 
                          data = buttabund_lowland2, control=glmerControl(optimizer="bobyqa"))
summary(mod2_lowland3) ## all fine and very similar results to convergence model


## check model
boxplot(resid(mod1_lowland) ~ buttabund_lowland$SURVEY_SQUARE)


# 1km
pred.data_lowland1 <- expand.grid(AES1KM=seq(min(buttabund_lowland$AES1KM), max(buttabund_lowland$AES1KM), 0.2), 
                                  AES3KM = mean(buttabund_lowland$AES3KM), 
                                  PERC_SUN = mean(buttabund_lowland$PERC_SUN),
                                  SURVEY_TEMP_SHADE = mean(buttabund_lowland$SURVEY_TEMP_SHADE), 
                                  NCA = NA, 
                                  ROUND_NUMBER = mean(buttabund_lowland$ROUND_NUMBER), 
                                  SURVEY_YEAR = factor(2018, levels = levels(buttabund_lowland$SURVEY_YEAR)))

# 3km
pred.data_lowland2 <- expand.grid(AES1KM=mean(buttabund_lowland$AES1KM), 
                                  AES3KM = seq(min(buttabund_lowland$AES3KM), max(buttabund_lowland$AES3KM), 0.2), 
                                  PERC_SUN = mean(buttabund_lowland$PERC_SUN),
                                  SURVEY_TEMP_SHADE = mean(buttabund_lowland$SURVEY_TEMP_SHADE), 
                                  NCA = NA, 
                                  ROUND_NUMBER = mean(buttabund_lowland$ROUND_NUMBER), 
                                  SURVEY_YEAR = factor(2018, levels = levels(buttabund_lowland$SURVEY_YEAR)))


## predict
pred.data_lowland1 <- pred_conf(mod = mod1_lowland2, pred_data = pred.data_lowland1, 
                                col_id = 'BUTTERFLY_COUNT', re_name = 'SURVEY_SQUARE')
pred.data_lowland2 <- pred_conf(mod = mod1_lowland2, pred_data = pred.data_lowland2, 
                                col_id = 'BUTTERFLY_COUNT', re_name = 'SURVEY_SQUARE')


## plotting (includes boxplot)
plot_results(predict_data = pred.data_lowland1,
             xvars = 'AES1KM', 
             x_var_names = 'Local level',
             yvar = 'BUTTERFLY_COUNT', 
             boxplot_x = 'NCA',
             original_data = buttabund_lowland,
             scaled_data = buttabund,
             unscale = TRUE,
             write = TRUE,
             lowland_only = TRUE,
             plot_dir = paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/lowland_only"))

plot_results(predict_data = pred.data_lowland2,
             xvars = 'AES3KM', 
             x_var_names = 'Landscape level',
             yvar = 'BUTTERFLY_COUNT', 
             boxplot_x = 'NCA', 
             original_data = buttabund_lowland,
             scaled_data = buttabund,
             unscale = TRUE,
             write = TRUE,
             lowland_only = TRUE,
             plot_dir = paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/lowland_only"))


#write out model results

results_lowland <- capture.output(summary(mod1_lowland2))

cat(paste0(analysis_name, " model results"), results_lowland, 
    file=paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/", "lowland_only/Lowland only ", analysis_name," model results.txt"), sep="\n", append=TRUE)


results_outlier_lowland <- capture.output(summary(mod2_lowland3))

cat(paste0(analysis_name, " model results - outlier removed"), results_outlier_lowland, 
    file=paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/", "lowland_only/Lowland only ", analysis_name," model results outlier removed.txt"), sep="\n", append=TRUE)


model_table(model = mod1_lowland2,
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

buttabund_habs <- merge(buttabund, covariates, by.x = c("SURVEY_SQUARE", "SURVEY_YEAR"), by.y = c("SURVEY_SQUARE", "SURVEY_YEAR"))

par(mfrow=c(2,2))
for (i in 11:22){
  plot(BUTTERFLY_COUNT ~ buttabund_habs[,i], data = buttabund_habs, xlab = names(buttabund_habs)[i])
}

mod3 <- glmer.nb(BUTTERFLY_COUNT ~ AES1KM*AES3KM + ROUND_NUMBER + SURVEY_YEAR + PERC_SUN + SURVEY_TEMP_SHADE + Habitat.diversity + BL_wood + Woody.linear + Water.linear + (1|NCA/SURVEY_SQUARE), data = buttabund_habs, control=glmerControl(optimizer="bobyqa"))

summary(mod3)



results_habitats <- capture.output(summary(mod3))


cat(paste0(analysis_name, " model results habitat analyses"), results_habitats, file=paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/", analysis_name," model results habitat analyses.txt"), sep="\n", append=FALSE)


#' plant models
#' 

mod4 <- glmer.nb(BUTTERFLY_COUNT ~ AES1KM*AES3KM + ROUND_NUMBER + SURVEY_YEAR + PERC_SUN + SURVEY_TEMP_SHADE + BOT_DIV + PERC_GRAM + FLORAL_RESOURCES + (1|NCA/SURVEY_SQUARE), data = buttabund_habs,  control=glmerControl(optimizer="bobyqa"))

summary(mod4)


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

buttabund_habslowland <- buttabund_habs[!buttabund_habs$NCA %in% c("Dartmoor", "Yorkshire\nDales"),]


par(mfrow=c(2,2))
for (i in 11:22){
  plot(BUTTERFLY_COUNT ~ buttabund_habslowland[,i], data = buttabund_habslowland, xlab = names(buttabund_habslowland)[i])
}

mod5 <- glmer.nb(BUTTERFLY_COUNT ~ AES1KM*AES3KM + ROUND_NUMBER + SURVEY_YEAR + PERC_SUN + SURVEY_TEMP_SHADE + Habitat.diversity + Water.linear + (1|NCA/SURVEY_SQUARE), data = buttabund_habslowland)

summary(mod5)



results_habitats <- capture.output(summary(mod5))


cat(paste0(analysis_name, " model results habitat analyses"), results_habitats, file=paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/", "lowland_only/Lowland only ", analysis_name," model results habitat analyses.txt"), sep="\n", append=FALSE)


#' plant models
#' 

mod6 <- glmer.nb(BUTTERFLY_COUNT ~ AES1KM*AES3KM + ROUND_NUMBER + SURVEY_YEAR + PERC_SUN + SURVEY_TEMP_SHADE + BOT_DIV + PERC_GRAM + FLORAL_RESOURCES + (1|SURVEY_SQUARE), data = buttabund_habslowland)

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



