##model validation

library(brms)

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

#RMSE

#k_Rich_LS <- kfold(Rich_LS_mod, K = 5)

sqrt(mean((fitted(Rich_LS_mod)[,1] - Rich_LS_mod$data$RICHNESS_ID)^2))
sqrt(mean((fitted(Rich_WCBS_mod)[,1] - Rich_WCBS_mod$data$Richness)^2))
sqrt(mean((fitted(Rich_UKBMS_mod)[,1] - Rich_UKBMS_mod$data$Richness)^2))

sqrt(mean((fitted(Abund_LS_mod)[,1] - Abund_LS_mod$data$BUTTERFLY_COUNT)^2))
sqrt(mean((fitted(Abund_WCBS_mod)[,1] - Abund_WCBS_mod$data$Abundance)^2))
sqrt(mean((fitted(Abund_UKBMS_mod)[,1] - Abund_UKBMS_mod$data$Abundance)^2))

sqrt(mean((fitted(Div_LS_mod)[,1] - Div_LS_mod$data$expSHANNON_DIV)^2))
sqrt(mean((fitted(Div_WCBS_mod)[,1] - Div_WCBS_mod$data$expDiversity)^2))
sqrt(mean((fitted(Div_UKBMS_mod)[,1] - Div_UKBMS_mod$data$expDiversity)^2))

#MAE
median(abs(fitted(Rich_LS_mod)[,1] - Rich_LS_mod$data$RICHNESS_ID))
median(abs(fitted(Rich_WCBS_mod)[,1] - Rich_WCBS_mod$data$Richness))
median(abs(fitted(Rich_UKBMS_mod)[,1] - Rich_UKBMS_mod$data$Richness))

median(abs(fitted(Abund_LS_mod)[,1] - Abund_LS_mod$data$BUTTERFLY_COUNT))
median(abs(fitted(Abund_WCBS_mod)[,1] - Abund_WCBS_mod$data$Abundance))
median(abs(fitted(Abund_UKBMS_mod)[,1] - Abund_UKBMS_mod$data$Abundance))

median(abs(fitted(Div_LS_mod)[,1] - Div_LS_mod$data$expSHANNON_DIV))
median(abs(fitted(Div_WCBS_mod)[,1] - Div_WCBS_mod$data$expDiversity))
median(abs(fitted(Div_UKBMS_mod)[,1] - Div_UKBMS_mod$data$expDiversity))

#CV

sqrt(mean((fitted(Rich_LS_mod)[,1] - Rich_LS_mod$data$RICHNESS_ID)^2))/mean(Rich_LS_mod$data$RICHNESS_ID)
sqrt(mean((fitted(Rich_WCBS_mod)[,1] - Rich_WCBS_mod$data$Richness)^2))/mean(Rich_WCBS_mod$data$Richness)
sqrt(mean((fitted(Rich_UKBMS_mod)[,1] - Rich_UKBMS_mod$data$Richness)^2))/mean(Rich_UKBMS_mod$data$Richness)

sqrt(mean((fitted(Abund_LS_mod)[,1] - Abund_LS_mod$data$BUTTERFLY_COUNT)^2))/mean(Abund_LS_mod$data$BUTTERFLY_COUNT)
sqrt(mean((fitted(Abund_WCBS_mod)[,1] - Abund_WCBS_mod$data$Abundance)^2))/mean(Abund_WCBS_mod$data$Abundance)
sqrt(mean((fitted(Abund_UKBMS_mod)[,1] - Abund_UKBMS_mod$data$Abundance)^2))/mean(Abund_UKBMS_mod$data$Abundance)

sqrt(mean((fitted(Div_LS_mod)[,1] - Div_LS_mod$data$expSHANNON_DIV)^2))/mean(Div_LS_mod$data$expSHANNON_DIV)
sqrt(mean((fitted(Div_WCBS_mod)[,1] - Div_WCBS_mod$data$expDiversity)^2))/mean(Div_WCBS_mod$data$expDiversity)
sqrt(mean((fitted(Div_UKBMS_mod)[,1] - Div_UKBMS_mod$data$expDiversity)^2))/mean(Div_UKBMS_mod$data$expDiversity)


#plots
plot(fitted(Abund_LS_mod)[,1] , Abund_LS_mod$data$BUTTERFLY_COUNT, xlim = c(0,10000), ylim = c(0,10000))
plot(fitted(Abund_WCBS_mod)[,1] , Abund_WCBS_mod$data$Abundance, xlim = c(0,10000), ylim = c(0,10000))
plot(fitted(Abund_UKBMS_mod)[,1] , Abund_UKBMS_mod$data$Abundance, xlim = c(0,10000), ylim = c(0,10000))
