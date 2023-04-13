##model validation

library(brms)

# folder setup for saving
dir <- config::get()
modpath <- dir$directories$models

my_col <- unname(palette.colors()[c(8,3,4)])


Rich_LS_mod <- readRDS(paste0(modpath, "LandSpAES_Richness_brm.RDS"))
Rich_WCBS_mod <- readRDS(paste0(modpath, "WCBS_Richness_brm.RDS"))
Rich_UKBMS_mod <- readRDS(paste0(modpath, "UKBMS_Richness_brm.RDS"))
Rich_Int_mod <- readRDS("Integrated_Richness_model.RDS")
Rich_Int_mod2 <- readRDS("Integrated_Richness_model_noUKBMS.RDS")

Abund_LS_mod <- readRDS(paste0(modpath, "LandSpAES_Abundance_brm.RDS"))
Abund_WCBS_mod <- readRDS(paste0(modpath, "WCBS_Abundance_brm.RDS"))
Abund_UKBMS_mod <- readRDS(paste0(modpath, "UKBMS_Abundance_brm.RDS"))
Abund_Int_mod <- readRDS("Integrated_Abundance_model.RDS")
Abund_Int_mod2 <- readRDS("Integrated_Abundance_model_noUKBMS.RDS")

Div_LS_mod <- readRDS(paste0(modpath, "LandSpAES_Diversity_brm.RDS"))
Div_WCBS_mod <- readRDS(paste0(modpath, "WCBS_Diversity_brm.RDS"))
Div_UKBMS_mod <- readRDS(paste0(modpath, "UKBMS_Diversity_brm.RDS"))
Div_Int_mod <- readRDS("Integrated_Diversity_model.RDS")
Div_Int_mod2 <- readRDS("Integrated_Diversity_model_noUKBMS.RDS")
#RMSE

k_Rich_LS <- kfold(Rich_LS_mod, K = 5, save_fits= TRUE, chains = 1)
k_Abund_LS <- kfold(Abund_LS_mod, K = 5, save_fits= TRUE, chains = 1)
k_Div_LS <- kfold(Div_LS_mod, K = 5, save_fits= TRUE, chains = 1)

# define a loss function
rmse <- function(y, yrep) {
  yrep_mean <- colMeans(yrep)
  sqrt(mean((yrep_mean - y)^2))
}

mae <- function(y, yrep) {
  yrep_mean <- colMeans(yrep)
  median(abs(yrep_mean - y))
}

cv <- function(y, yrep) {
  yrep_mean <- colMeans(yrep)
  sqrt(mean((yrep_mean - y)^2))/mean(y)
}


# predict responses and evaluate the loss
kfp <- kfold_predict(k_Rich_LS)
rmse(y = kfp$y, yrep = kfp$yrep)
mae(y = kfp$y, yrep = kfp$yrep)
cv(y = kfp$y, yrep = kfp$yrep)

kfp <- kfold_predict(k_Abund_LS)
rmse(y = kfp$y, yrep = kfp$yrep)
mae(y = kfp$y, yrep = kfp$yrep)
cv(y = kfp$y, yrep = kfp$yrep)

kfp <- kfold_predict(k_Div_LS)
rmse(y = kfp$y, yrep = kfp$yrep)
mae(y = kfp$y, yrep = kfp$yrep)
cv(y = kfp$y, yrep = kfp$yrep)

val_fun <- function(x, tr = "none"){
  MAE <- vector();rmse <- vector(); mean_x <- vector()
  xdat <- x$frame; resp <- x$call$formula[[2]]
  k_nos <- rep(1:5, each= nrow(x$frame)/5)
  k <- sample(k_nos, length(k_nos), replace = FALSE)
  
  for (i in 1:5){
    train_dat <- xdat[k != i,]
    test_dat <- xdat[k == i,]
    model <- NULL
    model <- glmmTMB(formula = x$call$formula, family = x$call$family, data = train_dat)
    if(is.null(model)){next} else {
      preds <- as.numeric(predict(model, test_dat, na.action = na.pass))
      if(tr == "log"){
        preds <- exp(as.numeric(predict(model, test_dat, na.action = na.pass)))
      }
      MAE[i] <- median(abs(test_dat[,names(test_dat) == resp] -preds))
      rmse[i] <- sqrt(mean((test_dat[,names(test_dat) == resp] - preds)^2))
      mean_x[i] <- mean(test_dat[,names(test_dat) == resp])
    } }
  meanMAE <- mean(MAE, na.rm = TRUE)
  meanRMSE <- mean(rmse, na.rm = TRUE)
  CV <- mean(rmse/mean_x, na.rm=TRUE)
  return(data.frame(MAE = meanMAE, RMSE = meanRMSE, CV = CV))
}

v1 <- val_fun(Abund_Int_mod, tr = "log")
v2 <- val_fun(Div_Int_mod)
v3 <- val_fun(Rich_Int_mod, tr = "log")

v1 <- val_fun(Abund_Int_mod2, tr = "log")
v2 <- val_fun(Div_Int_mod2)
v3 <- val_fun(Rich_Int_mod2, tr = "log")

r1 <- sqrt(mean((fitted(Rich_LS_mod)[,1] - Rich_LS_mod$data$RICHNESS_ID)^2))
r2 <- sqrt(mean((fitted(Rich_WCBS_mod)[,1] - Rich_WCBS_mod$data$Richness)^2))
r3 <- sqrt(mean((fitted(Rich_UKBMS_mod)[,1] - Rich_UKBMS_mod$data$Richness)^2))

r4 <- sqrt(mean((fitted(Abund_LS_mod)[,1] - Abund_LS_mod$data$BUTTERFLY_COUNT)^2))
r5 <- sqrt(mean((fitted(Abund_WCBS_mod)[,1] - Abund_WCBS_mod$data$Abundance)^2))
r6 <- sqrt(mean((fitted(Abund_UKBMS_mod)[,1] - Abund_UKBMS_mod$data$Abundance)^2))

r7 <- sqrt(mean((fitted(Div_LS_mod)[,1] - Div_LS_mod$data$expSHANNON_DIV)^2))
r8 <- sqrt(mean((fitted(Div_WCBS_mod)[,1] - Div_WCBS_mod$data$expDiversity)^2))
r9 <- sqrt(mean((fitted(Div_UKBMS_mod)[,1] - Div_UKBMS_mod$data$expDiversity)^2))

r10 <- sqrt(mean((fitted(Rich_Int_mod) - Rich_Int_mod$frame$Richness)^2))
r11 <- sqrt(mean((fitted(Abund_Int_mod) - Abund_Int_mod$frame$Abundance)^2))
r12 <- sqrt(mean((fitted(Div_Int_mod) - Div_Int_mod$frame$expShannon_diversity)^2))

#MAE
m1 <- median(abs(fitted(Rich_LS_mod)[,1] - Rich_LS_mod$data$RICHNESS_ID))
m2 <- median(abs(fitted(Rich_WCBS_mod)[,1] - Rich_WCBS_mod$data$Richness))
m3 <- median(abs(fitted(Rich_UKBMS_mod)[,1] - Rich_UKBMS_mod$data$Richness))

m4 <- median(abs(fitted(Abund_LS_mod)[,1] - Abund_LS_mod$data$BUTTERFLY_COUNT))
m5 <- median(abs(fitted(Abund_WCBS_mod)[,1] - Abund_WCBS_mod$data$Abundance))
m6 <- median(abs(fitted(Abund_UKBMS_mod)[,1] - Abund_UKBMS_mod$data$Abundance))

m7 <- median(abs(fitted(Div_LS_mod)[,1] - Div_LS_mod$data$expSHANNON_DIV))
m8 <- median(abs(fitted(Div_WCBS_mod)[,1] - Div_WCBS_mod$data$expDiversity))
m9 <- median(abs(fitted(Div_UKBMS_mod)[,1] - Div_UKBMS_mod$data$expDiversity))

m10 <- median(abs(fitted(Rich_Int_mod) - Rich_Int_mod$frame$Richness))
m11 <- median(abs(fitted(Abund_Int_mod) - Abund_Int_mod$frame$Abundance))
m12 <- median(abs(fitted(Div_Int_mod) - Div_Int_mod$frame$expShannon_diversity))


#CV

c1 <- sqrt(mean((fitted(Rich_LS_mod)[,1] - Rich_LS_mod$data$RICHNESS_ID)^2))/mean(Rich_LS_mod$data$RICHNESS_ID)
c2 <- sqrt(mean((fitted(Rich_WCBS_mod)[,1] - Rich_WCBS_mod$data$Richness)^2))/mean(Rich_WCBS_mod$data$Richness)
c3 <- sqrt(mean((fitted(Rich_UKBMS_mod)[,1] - Rich_UKBMS_mod$data$Richness)^2))/mean(Rich_UKBMS_mod$data$Richness)

c4 <- sqrt(mean((fitted(Abund_LS_mod)[,1] - Abund_LS_mod$data$BUTTERFLY_COUNT)^2))/mean(Abund_LS_mod$data$BUTTERFLY_COUNT)
c5 <- sqrt(mean((fitted(Abund_WCBS_mod)[,1] - Abund_WCBS_mod$data$Abundance)^2))/mean(Abund_WCBS_mod$data$Abundance)
c6 <- sqrt(mean((fitted(Abund_UKBMS_mod)[,1] - Abund_UKBMS_mod$data$Abundance)^2))/mean(Abund_UKBMS_mod$data$Abundance)

c7 <- sqrt(mean((fitted(Div_LS_mod)[,1] - Div_LS_mod$data$expSHANNON_DIV)^2))/mean(Div_LS_mod$data$expSHANNON_DIV)
c8 <- sqrt(mean((fitted(Div_WCBS_mod)[,1] - Div_WCBS_mod$data$expDiversity)^2))/mean(Div_WCBS_mod$data$expDiversity)
c9 <- sqrt(mean((fitted(Div_UKBMS_mod)[,1] - Div_UKBMS_mod$data$expDiversity)^2))/mean(Div_UKBMS_mod$data$expDiversity)

c10 <- sqrt(mean((fitted(Rich_Int_mod) - Rich_Int_mod$frame$Richness)^2))/mean(Rich_Int_mod$frame$Richness)
c11 <- sqrt(mean((fitted(Abund_Int_mod) - Abund_Int_mod$frame$Abundance)^2))/mean(Abund_Int_mod$frame$Abundance)
c12 <- sqrt(mean((fitted(Div_Int_mod) - Div_Int_mod$frame$expShannon_diversity)^2))/mean(Div_Int_mod$frame$expShannon_diversity)

#plots
plot(fitted(Abund_LS_mod)[,1] , Abund_LS_mod$data$BUTTERFLY_COUNT, xlim = c(0,10000), ylim = c(0,10000))
plot(fitted(Abund_WCBS_mod)[,1] , Abund_WCBS_mod$data$Abundance, xlim = c(0,10000), ylim = c(0,10000))
plot(fitted(Abund_UKBMS_mod)[,1] , Abund_UKBMS_mod$data$Abundance, xlim = c(0,10000), ylim = c(0,10000))

out <- data.frame(Dataset = rep(c("LandSpAES", "Integrated"),3), Model = rep(c("Abundance", "Diversity", "Richness"), each = 2),
                  MAE = c(m4, m11, m7, m12, m1, m10), RMSE = c(r4,r11, r7, r12, r1, r10), CV = c(c4, c11, c7, c12, c1, c10))
