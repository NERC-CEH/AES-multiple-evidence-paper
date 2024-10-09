##HPD tests


library(HDInterval)

# folder setup for saving
modpath <- getwd()


Rich_LS_mod <- readRDS(paste0(modpath, "LandSpAES_Richness_brm.RDS"))
Rich_WCBS_mod <- readRDS(paste0(modpath, "WCBS_Richness_brm.RDS"))
Rich_UKBMS_mod <- readRDS(paste0(modpath, "UKBMS_Richness_brm.RDS"))

Abund_LS_mod <- readRDS(paste0(modpath, "LandSpAES_Abundance_brm.RDS"))
Abund_WCBS_mod <- readRDS(paste0(modpath, "WCBS_Abundance_brm.RDS"))
Abund_UKBMS_mod <- readRDS(paste0(modpath, "UKBMS_Abundance_brm.RDS"))

Div_LS_mod <- readRDS(paste0(modpath, "LandSpAES_Diversity_brm.RDS"))
Div_WCBS_mod <- readRDS(paste0(modpath, "WCBS_Diversity_brm.RDS"))
Div_UKBMS_mod <- readRDS(paste0(modpath, "UKBMS_Diversity_brm.RDS"))


HPD_test <- function(mod1, mod2, terms){
  out <- data.frame()
  for(i in terms){
  mat1 <- as.matrix(mod1, variable = paste0("b_",i))
  mat2 <- as.matrix(mod2, variable = paste0("b_",i))
  sample1 <- sample(mat1, 4000, replace = FALSE)
  sample2 <- sample(mat2, 4000, replace = FALSE)
  diffs <- sample1-sample2
  int <- hdi(diffs)
  zero_test <- int[1] < 0 & int[2] > 0
  out <- rbind(out, int)
  names(out) <- c("Low", "High")
  }
  return(round(out,3))
  }

t1 <- HPD_test(Rich_LS_mod, Rich_WCBS_mod, c("AES1KM", "AES3KM", "AES1KM:AES3KM"))
t2 <- HPD_test(Rich_LS_mod, Rich_UKBMS_mod, c("AES1KM", "AES3KM", "AES1KM:AES3KM"))
t3 <- HPD_test(Rich_UKBMS_mod, Rich_WCBS_mod, c("AES1KM", "AES3KM", "AES1KM:AES3KM"))
t4 <- as.data.frame(rbind(t1,t2,t3))
t4$comp <- rep(c("LS-WCBS", "LS-UKBMS", "UKBMS-WCBS"),each = 3)
t4$term <- rep(c("AES1KM", "AES3KM", "AES1KM:AES3KM"),3)
t4$Response <- "Richness"
Rich_HPD <- t4

t1 <- HPD_test(Abund_LS_mod, Abund_WCBS_mod, c("AES1KM", "AES3KM", "AES1KM:AES3KM"))
t2 <- HPD_test(Abund_LS_mod, Abund_UKBMS_mod, c("AES1KM", "AES3KM", "AES1KM:AES3KM"))
t3 <- HPD_test(Abund_UKBMS_mod, Abund_WCBS_mod, c("AES1KM", "AES3KM", "AES1KM:AES3KM"))
t4 <- as.data.frame(rbind(t1,t2,t3))
t4$comp <- rep(c("LS-WCBS", "LS-UKBMS", "UKBMS-WCBS"),each = 3)
t4$term <- rep(c("AES1KM", "AES3KM", "AES1KM:AES3KM"),3)
t4$Response <- "Abundance"
Abund_HPD <- t4

t1 <- HPD_test(Div_LS_mod, Div_WCBS_mod, c("AES1KM", "AES3KM", "AES1KM:AES3KM"))
t2 <- HPD_test(Div_LS_mod, Div_UKBMS_mod, c("AES1KM", "AES3KM", "AES1KM:AES3KM"))
t3 <- HPD_test(Div_UKBMS_mod, Div_WCBS_mod, c("AES1KM", "AES3KM", "AES1KM:AES3KM"))
t4 <- as.data.frame(rbind(t1,t2,t3))
t4$comp <- rep(c("LS-WCBS", "LS-UKBMS", "UKBMS-WCBS"),each = 3)
t4$term <- rep(c("AES1KM", "AES3KM", "AES1KM:AES3KM"),3)
t4$Response <- "Diversity"
Div_HPD <- t4


all_HPD <- rbind(Rich_HPD, Div_HPD, Abund_HPD)
write.csv(all_HPD, "HPDI checks - intervals.csv")
