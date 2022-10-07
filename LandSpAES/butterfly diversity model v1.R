## Butterfly diversity models


analysis_group <- "Headline analyses"
analysis_name <- "Butterfly diversity"
plot_name <- "Butterfly Shannon diversity index"


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
library(nlme)
library(tidyverse)


##Now aggregate data for models



buttcount2$RICHNESS_ID <- buttcount2$BUTTERFLY_ID

#remove "other butterflies" id

buttcount2$RICHNESS_ID[buttcount2$RICHNESS_ID == 65] <- NA


length(buttcount2$RICHNESS_ID[is.na(buttcount2$RICHNESS_ID)])
# 1253/14266 no butterflies



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

buttrounds <- aggregate(ROUND_NUMBER ~ SURVEY_SQUARE + SURVEY_YEAR, data = buttcount2, FUN = function(x) length(unique(x)))

buttdiv <- merge(buttdiv, buttrounds, by.x = c("SURVEY_SQUARE", "SURVEY_YEAR"), by.y = c("SURVEY_SQUARE", "SURVEY_YEAR"))

#sunshine
buttsun <- aggregate(PERC_SUN ~ SURVEY_SQUARE + SURVEY_YEAR, data = buttvariables, FUN = function(x) mean(x, na.rm = TRUE))

buttdiv <- merge(buttdiv, buttsun, by.x = c("SURVEY_SQUARE", "SURVEY_YEAR"), by.y = c("SURVEY_SQUARE", "SURVEY_YEAR"))

#temperature
butttemp <- aggregate(SURVEY_TEMP_SHADE ~ SURVEY_SQUARE + SURVEY_YEAR, data = buttvisit, FUN = function(x) mean(x, na.rm = TRUE))

buttdiv <- merge(buttdiv, butttemp, by.x = c("SURVEY_SQUARE", "SURVEY_YEAR"), by.y = c("SURVEY_SQUARE", "SURVEY_YEAR"))



#check diversity is normally distributed
hist(buttdiv$SHANNON_DIV)
#normalish, slightly left skewed
#square transform
#buttdiv$SHANNON_DIV_2 <- (buttdiv$SHANNON_DIV)^2



#scale AES variables
buttdiv$AES1KM <- scale(buttdiv$AES1KM)
buttdiv$AES3KM <- scale(buttdiv$AES3KM)

# scale temp and sun
buttdiv$SURVEY_TEMP_SHADE <- scale(buttdiv$SURVEY_TEMP_SHADE)
buttdiv$PERC_SUN <- scale(buttdiv$PERC_SUN)

#treat survey year as factor
buttdiv$SURVEY_YEAR <- factor(buttdiv$SURVEY_YEAR)

#scatterplots

plot(buttdiv$SHANNON_DIV ~ buttdiv$AES1KM)
plot(buttdiv$SHANNON_DIV ~ buttdiv$AES3KM)

boxplot(buttdiv$SHANNON_DIV ~ buttdiv$ROUND_NUMBER)

#model with the nlme package not lme4 as the response is normally distributed
mod1 <- lme(SHANNON_DIV ~ AES1KM*AES3KM + SURVEY_TEMP_SHADE + PERC_SUN + ROUND_NUMBER + SURVEY_YEAR, 
            random = ~1|NCA/SURVEY_SQUARE, data = buttdiv)
summary(mod1)

plot(mod1)

# check model with year as random effect
mod1_yearre <- lmer(SHANNON_DIV ~ AES1KM*AES3KM + SURVEY_TEMP_SHADE + PERC_SUN + ROUND_NUMBER + 
                      (1|NCA/SURVEY_SQUARE) + (1|SURVEY_YEAR), data = buttdiv)
summary(mod1_yearre)


#check the effect of removing the outlier
buttdiv2 <- buttdiv[buttdiv$SURVEY_SQUARE != 'SP4469',]

mod2 <- lme(SHANNON_DIV ~ AES1KM*AES3KM + SURVEY_TEMP_SHADE + PERC_SUN + ROUND_NUMBER + SURVEY_YEAR, 
            random = ~1|NCA/SURVEY_SQUARE, data = buttdiv2)

summary(mod2)



#check year effect

boxplot(resid(mod1) ~ buttdiv$SURVEY_SQUARE)
# no problem with survey square

# ## Plots
#' 
#' # change factors for boxplot
buttdiv$NCA <- fct_recode(buttdiv$NCA, "High\nWeald" = "High Weald", "South Suffolk\nand North\nEssex Clayland" = "South Suffolk and North Essex Clayland", 
                          "Dunsmore\nand Feldon" = "Dunsmore and Feldon", "The\nFens" = "The Fens", 
                          "Yorkshire\nDales" = "Yorkshire Dales")
buttdiv$NCA <- factor(buttdiv$NCA,
                      levels = c("High\nWeald", "South Suffolk\nand North\nEssex Clayland", 
                                 "Dunsmore\nand Feldon", "The\nFens", 
                                 "Dartmoor", "Yorkshire\nDales"))



pred.data1 = expand.grid(AES1KM=seq(min(buttdiv$AES1KM), max(buttdiv$AES1KM), 0.2), 
                         AES3KM = mean(buttdiv$AES3KM), 
                         PERC_SUN = mean(buttdiv$PERC_SUN),
                         SURVEY_TEMP_SHADE = mean(buttdiv$SURVEY_TEMP_SHADE), 
                         NCA = NA, 
                         ROUND_NUMBER = mean(buttdiv$ROUND_NUMBER), 
                         SURVEY_YEAR = factor(2018, levels = levels(buttdiv$SURVEY_YEAR)))

pred.data2 = expand.grid(AES1KM=mean(buttdiv$AES1KM), 
                         AES3KM = seq(min(buttdiv$AES3KM), max(buttdiv$AES3KM), 0.2), 
                         PERC_SUN = mean(buttdiv$PERC_SUN),
                         SURVEY_TEMP_SHADE = mean(buttdiv$SURVEY_TEMP_SHADE), 
                         NCA = NA, 
                         ROUND_NUMBER = mean(buttdiv$ROUND_NUMBER), 
                         SURVEY_YEAR = factor(2018, levels = levels(buttdiv$SURVEY_YEAR)))


## predict
pred.data1 <- pred_conf(mod = mod1, pred_data = pred.data1, col_id = 'SHANNON_DIV', re_name = 'NCA')
pred.data2 <- pred_conf(mod = mod1, pred_data = pred.data2, col_id = 'SHANNON_DIV', re_name = 'NCA')


## plotting (includes boxplot)
plot_results(predict_data = pred.data1,
             xvars = 'AES1KM', 
             x_var_names = 'Local level',
             yvar = 'SHANNON_DIV', 
             boxplot_x = 'NCA',
             ylim = c(0, 3),
             original_data = buttdiv,
             unscale = TRUE,
             write = TRUE,
             plot_dir = paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/"))

plot_results(predict_data = pred.data2,
             xvars = 'AES3KM', 
             x_var_names = 'Landscape level',
             yvar = 'SHANNON_DIV', 
             boxplot_x = 'NCA',
             ylim = c(0, 3), 
             original_data = buttdiv,
             unscale = TRUE,
             write = TRUE,
             plot_dir = paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/"))



#write out model results

results <- capture.output(summary(mod1))

cat(paste0(analysis_name, " model results"), results, file=paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/", analysis_name," model results.txt"), sep="\n", append=TRUE)



results_outlier <- capture.output(summary(mod2))

cat(paste0(analysis_name, " model results outlier removed"), results, file=paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/", analysis_name," model results outlier removed.txt"), sep="\n", append=TRUE)



model_table(model = mod1,
            NCA = 'All NCA',
            outlier = 'All data')

model_table(model = mod2,
            NCA = 'All NCA',
            outlier = 'Outlier removed')

############################
####    Lowland only    ####
############################


## create lowland directory
dir.create(paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/lowland_only/"))

buttdiv_lowland <- buttdiv[!buttdiv$NCA %in% c("Dartmoor", "Yorkshire\nDales"),]





#scatterplots

plot(buttdiv_lowland$SHANNON_DIV ~ buttdiv_lowland$AES1KM)
plot(buttdiv_lowland$SHANNON_DIV ~ buttdiv_lowland$AES3KM)

boxplot(buttdiv_lowland$SHANNON_DIV ~ buttdiv_lowland$ROUND_NUMBER)

#model with the nlme package not lme4 as the response is normally distributed
mod1_lowland <- lme(SHANNON_DIV ~ AES1KM*AES3KM + SURVEY_TEMP_SHADE + PERC_SUN + ROUND_NUMBER + SURVEY_YEAR, 
                    random = ~1|NCA/SURVEY_SQUARE, data = buttdiv_lowland)
summary(mod1_lowland)

# all fine

plot(mod1_lowland)

#check the effect of removing the outlier
buttdiv2_lowland <- buttdiv_lowland[buttdiv_lowland$SURVEY_SQUARE != 'SP4469',]

mod2_lowland <- lme(SHANNON_DIV ~ AES1KM*AES3KM + SURVEY_TEMP_SHADE + PERC_SUN + ROUND_NUMBER + SURVEY_YEAR, 
                    random = ~1|NCA/SURVEY_SQUARE, data = buttdiv2_lowland)

summary(mod2_lowland)

# all fine
#check year effect

boxplot(resid(mod1_lowland) ~ buttdiv_lowland$SURVEY_SQUARE)


pred.data1_lowland = expand.grid(AES1KM=seq(min(buttdiv_lowland$AES1KM), max(buttdiv_lowland$AES1KM), 0.2), 
                                 AES3KM = mean(buttdiv_lowland$AES3KM), 
                                 PERC_SUN = mean(buttdiv_lowland$PERC_SUN),
                                 SURVEY_TEMP_SHADE = mean(buttdiv_lowland$SURVEY_TEMP_SHADE), 
                                 NCA = NA, 
                                 ROUND_NUMBER = mean(buttdiv_lowland$ROUND_NUMBER), 
                                 SURVEY_YEAR = factor(2018, levels = levels(buttdiv_lowland$SURVEY_YEAR)))

pred.data2_lowland = expand.grid(AES1KM=mean(buttdiv_lowland$AES1KM), 
                                 AES3KM = seq(min(buttdiv_lowland$AES3KM), max(buttdiv_lowland$AES3KM), 0.2), 
                                 PERC_SUN = mean(buttdiv_lowland$PERC_SUN),
                                 SURVEY_TEMP_SHADE = mean(buttdiv_lowland$SURVEY_TEMP_SHADE), 
                                 NCA = NA, 
                                 ROUND_NUMBER = mean(buttdiv_lowland$ROUND_NUMBER), 
                                 SURVEY_YEAR = factor(2018, levels = levels(buttdiv_lowland$SURVEY_YEAR)))


## predict
pred.data1_lowland <- pred_conf(mod = mod1_lowland, pred_data = pred.data1_lowland, col_id = 'SHANNON_DIV', re_name = 'NCA')
pred.data2_lowland <- pred_conf(mod = mod1_lowland, pred_data = pred.data2_lowland, col_id = 'SHANNON_DIV', re_name = 'NCA')


## plotting (includes boxplot)
plot_results(predict_data = pred.data1_lowland,
             xvars = 'AES1KM', 
             x_var_names = 'Local level',
             yvar = 'SHANNON_DIV', 
             boxplot_x = 'NCA',
             ylim = c(0, 3),
             original_data = buttdiv_lowland,
             scaled_data = buttdiv,
             lowland_only = TRUE,
             unscale = TRUE,
             write = TRUE,
             plot_dir = paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/lowland_only"))

plot_results(predict_data = pred.data2,
             xvars = 'AES3KM', 
             x_var_names = 'Landscape level',
             yvar = 'SHANNON_DIV', 
             boxplot_x = 'NCA',
             ylim = c(0, 3), 
             original_data = buttdiv_lowland,
             scaled_data = buttdiv,
             lowland_only = TRUE,
             unscale = TRUE,
             write = TRUE,
             plot_dir = paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/lowland_only"))



#write out model results

results <- capture.output(summary(mod1_lowland))

cat(paste0(analysis_name, " model results"), results, file=paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/lowland_only/", analysis_name," model results.txt"), sep="\n", append=TRUE)



results_outlier <- capture.output(summary(mod2_lowland))

cat(paste0(analysis_name, " model results outlier removed"), results, file=paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/lowland_only/", analysis_name," model results outlier removed.txt"), sep="\n", append=TRUE)


model_table(model = mod1_lowland,
            NCA = 'Lowland only',
            outlier = 'All data')

model_table(model = mod2_lowland,
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

buttdiv_habs <- merge(buttdiv, covariates, by.x = c("SURVEY_SQUARE", "SURVEY_YEAR"), by.y = c("SURVEY_SQUARE", "SURVEY_YEAR"))

par(mfrow=c(2,2))
for (i in 12:22){
  plot(SHANNON_DIV ~ buttdiv_habs[,i], data = buttdiv_habs, xlab = names(buttdiv_habs)[i])
}

mod3 <- lme(SHANNON_DIV ~ AES1KM*AES3KM + ROUND_NUMBER + SURVEY_YEAR + PERC_SUN + SURVEY_TEMP_SHADE + Habitat.diversity + BL_wood + Woody.linear + Water.linear, random = ~1|NCA/SURVEY_SQUARE, data = buttdiv_habs)

summary(mod3)


results_habitats <- capture.output(summary(mod3))


cat(paste0(analysis_name, " model results habitat analyses"), results_habitats, file=paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/", analysis_name," model results habitat analyses.txt"), sep="\n", append=FALSE)


#' plant models
#' 

mod4 <- lme(SHANNON_DIV ~ AES1KM*AES3KM + ROUND_NUMBER + SURVEY_YEAR + PERC_SUN + SURVEY_TEMP_SHADE + BOT_DIV + PERC_GRAM + FLORAL_RESOURCES, random = ~1|NCA/SURVEY_SQUARE, data = buttdiv_habs)

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

buttdiv_habslowland <- buttdiv_habs[!buttdiv_habs$NCA %in% c("Dartmoor", "Yorkshire\nDales"),]


par(mfrow=c(2,2))
for (i in 11:22){
  plot(SHANNON_DIV ~ buttdiv_habslowland[,i], data = buttdiv_habslowland, xlab = names(buttdiv_habslowland)[i])
}

mod5 <- lme(SHANNON_DIV ~ AES1KM*AES3KM + ROUND_NUMBER + SURVEY_YEAR + PERC_SUN + SURVEY_TEMP_SHADE + Habitat.diversity + Water.linear, random = ~1|NCA/SURVEY_SQUARE, data = buttdiv_habslowland)


summary(mod5)

results_habitats <- capture.output(summary(mod5))


cat(paste0(analysis_name, " model results habitat analyses"), results_habitats, file=paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/", "lowland_only/Lowland only ", analysis_name," model results habitat analyses.txt"), sep="\n", append=FALSE)


#' plant models
#' 

mod6 <- lme(SHANNON_DIV ~ AES1KM*AES3KM + ROUND_NUMBER + SURVEY_YEAR + PERC_SUN + SURVEY_TEMP_SHADE + BOT_DIV + PERC_GRAM + FLORAL_RESOURCES, random = ~1|NCA/SURVEY_SQUARE, data = buttdiv_habslowland)

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


mod1_nosq <- update(mod1, . ~ . , random = ~1|NCA)
mod2_nosq <- update(mod2, . ~ . , random = ~1|NCA)
mod1_lowland_nosq <- update(mod1_lowland, . ~ . , random = ~1|NCA)
mod2_lowland_nosq <- update(mod2_lowland, . ~ . , random = ~1|NCA)
mod3_nosq <- update(mod3, . ~ . , random = ~1|NCA)
mod4_nosq <- update(mod4, . ~ . , random = ~1|NCA)
mod5_nosq <- update(mod5, . ~ . , random = ~1|NCA)
mod6_nosq <- update(mod6, . ~ . , random = ~1|NCA)

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



