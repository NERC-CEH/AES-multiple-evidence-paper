##butterfly trait abundance models


homedir <- getwd()

analysis_group <- "Butterfly trait analyses"
analysis_name <- "Butterfly abundance - high mobility"
plot_name <- "Total count of high mobility butterflies"


# folder setup for saving
dir <- config::get()
file.dir <- paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/")
dir.create(file.dir, recursive = TRUE)


###Basic richness models
source("collate butterfly data.R")

# source function
source("functions_script.R")



##Now aggregate data for models

buttcount2$RICHNESS_ID <- buttcount2$BUTTERFLY_ID

#remove "other butterflies" id

buttcount2$RICHNESS_ID[buttcount2$RICHNESS_ID == 65] <- NA


length(buttcount2$RICHNESS_ID[is.na(buttcount2$RICHNESS_ID)])
# 1253/14266 no butterflies




##Select trait grouping of interest

db <- config::get("AES")

con <- odbcConnect(db$DSN, db$uid, db$pwd, believeNRows = FALSE)

butttraits <- sqlFetch(con, "TBL_BUTTERFLY_SPECIES")

close(con)


buttcount2$TRAIT <- butttraits$WINGSPAN_CATEGORY[match(buttcount2$BUTTERFLY_ID, butttraits$BUTTERFLY_ID)]

buttcount2$BUTTERFLY_COUNT <- ifelse(buttcount2$TRAIT == 3, buttcount2$BUTTERFLY_COUNT, NA)



#' **Calculate BUTTERFLY abundance - aggregated** - remove crop pests

buttabund <- aggregate(BUTTERFLY_COUNT ~ NCA + SURVEY_SQUARE + SURVEY_YEAR + AES1KM + AES3KM, data = buttcount2, FUN = function(x) sum(x))
#36 observations

buttrounds <- aggregate(ROUND_NUMBER ~ NCA + SURVEY_SQUARE + SURVEY_YEAR + AES1KM + AES3KM, data = buttcount2, FUN = function(x) length(unique(x)))

buttabund <- merge(buttabund, buttrounds, by.x = c("NCA", "SURVEY_SQUARE", "SURVEY_YEAR", "AES1KM", "AES3KM"), by.y = c("NCA", "SURVEY_SQUARE", "SURVEY_YEAR", "AES1KM", "AES3KM"), all = TRUE)

buttabund$BUTTERFLY_COUNT[is.na(buttabund$BUTTERFLY_COUNT)] <- 0

#sunshine
buttsun <- aggregate(PERC_SUN ~ SURVEY_SQUARE + SURVEY_YEAR, data = buttvariables, FUN = function(x) mean(x, na.rm = TRUE))

buttabund <- merge(buttabund, buttsun, by.x = c("SURVEY_SQUARE", "SURVEY_YEAR"), by.y = c("SURVEY_SQUARE", "SURVEY_YEAR"))

#temperature
butttemp <- aggregate(SURVEY_TEMP_SHADE ~ SURVEY_SQUARE + SURVEY_YEAR, data = buttvisit, FUN = function(x) mean(x, na.rm = TRUE))

buttabund <- merge(buttabund, butttemp, by.x = c("SURVEY_SQUARE", "SURVEY_YEAR"), by.y = c("SURVEY_SQUARE", "SURVEY_YEAR"))


#' Look at distribution of abundance variable
fit_distrib(buttabund, 'BUTTERFLY_COUNT')
## definitely negative binomial

#boxplot with year

boxplot(buttabund$BUTTERFLY_COUNT ~ buttabund$SURVEY_YEAR)

#' Rescale predictors and fit model
#' 
buttabund$AES1KM <- scale(buttabund$AES1KM)
buttabund$AES3KM <- scale(buttabund$AES3KM)

buttabund$SURVEY_TEMP_SHADE <- scale(buttabund$SURVEY_TEMP_SHADE)
buttabund$PERC_SUN <- scale(buttabund$PERC_SUN)

#treat survey year as factor
buttabund$SURVEY_YEAR <- factor(buttabund$SURVEY_YEAR)


library(lme4)

# models
mod1 <- glmer.nb(BUTTERFLY_COUNT ~ AES1KM*AES3KM + ROUND_NUMBER + PERC_SUN + SURVEY_TEMP_SHADE + SURVEY_YEAR + 
                   (1|NCA/SURVEY_SQUARE), data = buttabund)

summary(mod1)



#remove outlier
buttabund2 <- buttabund[buttabund$SURVEY_SQUARE != 'SP4469',]

mod2 <- glmer.nb(BUTTERFLY_COUNT ~ AES1KM*AES3KM + ROUND_NUMBER  + PERC_SUN + SURVEY_TEMP_SHADE + SURVEY_YEAR +  
                   (1|NCA/SURVEY_SQUARE), data = buttabund2)

summary(mod2)



# ## Plots
#' 
#' # change factors for boxplot
#' 

library(tidyverse)

buttabund$NCA <- fct_recode(buttabund$NCA, "High\nWeald" = "High Weald", "South Suffolk\nand North\nEssex Clayland" = "South Suffolk and North Essex Clayland", 
                            "Dunsmore\nand Feldon" = "Dunsmore and Feldon", "The\nFens" = "The Fens", 
                            "Yorkshire\nDales" = "Yorkshire Dales")
buttabund$NCA <- factor(buttabund$NCA,
                        levels = c("High\nWeald", "South Suffolk\nand North\nEssex Clayland", 
                                   "Dunsmore\nand Feldon", "The\nFens", 
                                   "Dartmoor", "Yorkshire\nDales"))


pred.data1 = expand.grid(AES1KM=seq(min(buttabund$AES1KM), max(buttabund$AES1KM), 0.2), 
                         AES3KM = mean(buttabund$AES3KM), 
                         PERC_SUN = mean(buttabund$PERC_SUN),
                         SURVEY_TEMP_SHADE = mean(buttabund$SURVEY_TEMP_SHADE),
                         NCA = NA, 
                         ROUND_NUMBER = mean(buttabund$ROUND_NUMBER), 
                         SURVEY_YEAR = factor(2018, levels = levels(buttabund$SURVEY_YEAR)))
pred.data2 = expand.grid(AES1KM=mean(buttabund$AES1KM), 
                         AES3KM = seq(min(buttabund$AES3KM), max(buttabund$AES3KM), 0.2),
                         PERC_SUN = mean(buttabund$PERC_SUN),
                         SURVEY_TEMP_SHADE = mean(buttabund$SURVEY_TEMP_SHADE),
                         NCA = NA, ROUND_NUMBER = mean(buttabund$ROUND_NUMBER), 
                         SURVEY_YEAR = factor(2018, levels = levels(buttabund$SURVEY_YEAR)))


## predict
pred.data1 <- pred_conf(mod = mod1, pred_data = pred.data1, col_id = 'BUTTERFLY_COUNT', re_name = 'NCA')
pred.data2 <- pred_conf(mod = mod1, pred_data = pred.data2, col_id = 'BUTTERFLY_COUNT', re_name = 'NCA')


## plotting (includes boxplot)
plot_results(predict_data = pred.data1,
             xvars = 'AES1KM', 
             x_var_names = 'Local level',
             yvar = 'BUTTERFLY_COUNT', 
             boxplot_x = 'NCA',
             ylim = c(0, 1000),
             original_data = buttabund,
             unscale = TRUE,
             write = TRUE,
             plot_dir = paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/"))

plot_results(predict_data = pred.data2,
             xvars = 'AES3KM', 
             x_var_names = 'Landscape level',
             yvar = 'BUTTERFLY_COUNT', 
             boxplot_x = 'NCA',
             ylim = c(0, 1000), 
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

############################
####    Lowland only    ####
############################


## create lowland directory
dir.create(paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/lowland_only"))


buttabund_lowland <- buttabund[!buttabund$NCA %in% c("Dartmoor", "Yorkshire\nDales"),]


#' Plot histogram of data
#' 
fit_distrib(buttabund_lowland, 'BUTTERFLY_COUNT')


#model with the nlme package not lme4 as the response is normally distributed
mod1_lowland <- glmer.nb(BUTTERFLY_COUNT ~ AES1KM*AES3KM + ROUND_NUMBER + PERC_SUN + SURVEY_TEMP_SHADE + SURVEY_YEAR + 
                           (1|NCA/SURVEY_SQUARE), 
                         data = buttabund_lowland)

summary(mod1_lowland)



#remove outlier
buttabund2_lowland <- buttabund_lowland[buttabund_lowland$SURVEY_SQUARE != 'SP4469',]

mod2_lowland <- glmer.nb(BUTTERFLY_COUNT ~ AES1KM*AES3KM + ROUND_NUMBER  + PERC_SUN + SURVEY_TEMP_SHADE + SURVEY_YEAR +  
                           (1|NCA/SURVEY_SQUARE), data = buttabund2_lowland)

summary(mod2_lowland)





pred.data1_lowland = expand.grid(AES1KM=seq(min(buttabund_lowland$AES1KM), max(buttabund_lowland$AES1KM), 0.2), 
                                 AES3KM = mean(buttabund_lowland$AES3KM), 
                                 PERC_SUN = mean(buttabund_lowland$PERC_SUN),
                                 SURVEY_TEMP_SHADE = mean(buttabund_lowland$SURVEY_TEMP_SHADE),
                                 NCA = NA, 
                                 ROUND_NUMBER = mean(buttabund_lowland$ROUND_NUMBER), 
                                 SURVEY_YEAR = factor(2018, levels = levels(buttabund_lowland$SURVEY_YEAR)))

pred.data2_lowland = expand.grid(AES1KM=mean(buttabund_lowland$AES1KM), 
                                 AES3KM = seq(min(buttabund_lowland$AES3KM), max(buttabund_lowland$AES3KM), 0.2),
                                 PERC_SUN = mean(buttabund_lowland$PERC_SUN),
                                 SURVEY_TEMP_SHADE = mean(buttabund_lowland$SURVEY_TEMP_SHADE),
                                 NCA = NA, ROUND_NUMBER = mean(buttabund_lowland$ROUND_NUMBER), 
                                 SURVEY_YEAR = factor(2018, levels = levels(buttabund_lowland$SURVEY_YEAR)))


## predict
pred.data1_lowland <- pred_conf(mod = mod1_lowland, pred_data = pred.data1_lowland, col_id = 'BUTTERFLY_COUNT', re_name = 'NCA')
pred.data2_lowland <- pred_conf(mod = mod1_lowland, pred_data = pred.data2_lowland, col_id = 'BUTTERFLY_COUNT', re_name = 'NCA')


## plotting (includes boxplot)
### Won't save because of long file names

plot_results(predict_data = pred.data1_lowland,
             xvars = 'AES1KM', 
             x_var_names = 'Local level',
             yvar = 'BUTTERFLY_COUNT', 
             boxplot_x = 'NCA',
             ylim = c(0, 1000),
             original_data = buttabund_lowland,
             scaled_data = buttabund,
             lowland_only = TRUE,
             unscale = TRUE,
             write = FALSE,
             plot_dir = paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/lowland_only"))

plot_results(predict_data = pred.data2,
             xvars = 'AES3KM', 
             x_var_names = 'Landscape level',
             yvar = 'BUTTERFLY_COUNT', 
             boxplot_x = 'NCA',
             ylim = c(0, 1000), 
             original_data = buttabund_lowland,
             scaled_data = buttabund,
             lowland_only = TRUE,
             unscale = TRUE,
             write = FALSE,
             plot_dir = paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/lowland_only/"))


# save model outputs 

results_lowland <- capture.output(summary(mod1_lowland))

cat(paste0(analysis_name, " model results"), results_lowland, 
    file=paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/", "lowland_only/Lowland only ", analysis_name," model results.txt"), sep="\n", append=TRUE)


results_outlier_lowland <- capture.output(summary(mod2_lowland))

cat(paste0(analysis_name, " model results - outlier removed"), results_outlier_lowland, 
    file=paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/", "lowland_only/Lowland only ", analysis_name," model results outlier removed.txt"), sep="\n", append=TRUE)



model_table(model = mod1_lowland,
            NCA = 'Lowland only',
            outlier = 'All data')

model_table(model = mod2_lowland,
            NCA = 'Lowland only',
            outlier = 'Outlier removed')



####################################################
#####        No square level random term       #####
####################################################


mod1_nosq <- update(mod1, . ~ . -(1|NCA/SURVEY_SQUARE) -(1|SURVEY_SQUARE) + (1|NCA))
mod2_nosq <- update(mod2, . ~ . -(1|NCA/SURVEY_SQUARE)  -(1|SURVEY_SQUARE) + (1|NCA))
mod1_lowland_nosq <- update(mod1_lowland, . ~ . -(1|NCA/SURVEY_SQUARE)  -(1|SURVEY_SQUARE) + (1|NCA))
mod2_lowland_nosq <- update(mod2_lowland, . ~ . -(1|NCA/SURVEY_SQUARE)  -(1|SURVEY_SQUARE) + (1|NCA))

results_mod1_nosq <- capture.output(summary(mod1_nosq))
results_mod2_nosq <- capture.output(summary(mod2_nosq))
results_mod1_lowland_nosq <- capture.output(summary(mod1_lowland_nosq))
results_mod2_lowland_nosq <- capture.output(summary(mod2_lowland_nosq))



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