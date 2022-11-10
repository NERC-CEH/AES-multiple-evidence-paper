#' # Ordination of NCAs
#' 
library(readxl)
library(ggcorrplot)

#' ## Collate data
#' 
#' set path to landscape and env data
#' 
dir <- config::get()
fpath <- dir$directories$envdata
#' 
#' 
#' ## 1km square level data:
#' 
#' Mean elevation
#' 
mean_elev <- read.csv(paste0(fpath, "IHDTM_1km/elevation_1km.csv"))
#' 
#' SD elevation
#' 
SD_elev <- read.csv(paste0(fpath, "IHDTM_1km/elevation_sd_1km.csv"))
#' 
#' Climate
#'
climate <- read.csv(paste0(fpath,"HADUK/HADUK_Rain_Temp_Summaries.csv"))

#need to add the 1km grid refs to climate data to allow matching to other datasets, have to do this via the ukbms site info
#also need the scores
fpath_ukbms <- dir$directories$ukbmsdata
wcbs_ukbms <- read.table(paste0(fpath_ukbms, "site_date_2017-2021.txt"), sep = "\t", header = TRUE)
fpath_scores <- dir$directories$scoredata
CSlocs <- read.csv(paste0(fpath_scores, "Butts_Gradient_Scores.csv"))

climate$SITENO <- wcbs_ukbms$SITENO[match(climate$GRIDREF, wcbs_ukbms$GRIDREF)]
climate$PLAN_NO <- CSlocs$CELLCODE[match(climate$SITENO, CSlocs$buttsurv.SITENO)]
climate$PLAN_NO[is.na(climate$PLAN_NO)] <- climate$GRIDREF[is.na(climate$PLAN_NO)]


#' 
#' Cover LCM (separate file per year)
#' 
LCM17 <- read.csv(paste0(fpath, "LCM/lcm17_1km.csv"))
LCM18 <- read.csv(paste0(fpath, "LCM/lcm18_1km.csv"))
LCM19 <- read.csv(paste0(fpath, "LCM/lcm19_1km.csv"))
LCM21 <- read.csv(paste0(fpath, "LCM/lcm21_1km.csv"))
#' 

#'
#' Length of hedgerows - note now does not include 0 lengths
#' 
#hedges_old <- read.csv(paste0(fpath, "Hedgerows/hedge_df_1km.csv"))
hedges <- read.csv(paste0(fpath, "Hedgerows/wlf_df_1km.csv"))
#' 
#' Slope
#' 
slope <- read.csv(paste0(fpath, "IHDTM_1km/slope_1km.csv"))
#' 
#' Aspect (southness)
#' 
aspect_south <- read.csv(paste0(fpath, "IHDTM_1km/southness_1km.csv"))
#' 
#' Aspect (eastness)
#' 
aspect_east <- read.csv(paste0(fpath, "IHDTM_1km/eastness_1km.csv"))
#' 
#' 
#' Parent material
#' 
soils <- read.csv(paste0(fpath, "SoilParentMaterial_V1_portal1km/soilparent_1km.csv"))
#' 
#' NCAs
#' 
#NCAs <- read.csv(paste0(fpath, "NCAs/nca_1km_dom.csv"))
#' 
#' 

#' 
#' 
#' 
#' ##Data to derive:
#' 
#' function to extract data from LCM sheets (for ordination use average across 3 years)

lcm_extract <- function(x, hab_name){
  x17 <- LCM17[LCM17$gridcode %in% x,]
  x18 <- LCM18[LCM18$gridcode %in% x,]
  x19 <- LCM19[LCM19$gridcode %in% x,]
  x21 <- LCM21[LCM21$gridcode %in% x,]
  if(length(x) > 1){
    x17 <- aggregate(cbind(area_msq, area_m_pc) ~ PLAN_NO, FUN = sum, data = x17)
    x18 <- aggregate(cbind(area_msq, area_m_pc) ~ PLAN_NO, FUN = sum, data = x18)
    x19 <- aggregate(cbind(area_msq, area_m_pc) ~ PLAN_NO, FUN = sum, data = x19)
    x21 <- aggregate(cbind(area_msq, area_pc) ~ PLAN_NO, FUN = sum, data = x21)
  }
  names(x17)[grep("area", names(x17))] <- paste0(names(x17)[grep("area", names(x17))], "_17")
  names(x18)[grep("area", names(x18))] <- paste0(names(x18)[grep("area", names(x18))], "_18")
  names(x19)[grep("area", names(x19))] <- paste0(names(x19)[grep("area", names(x19))], "_19")
  names(x21)[grep("area", names(x21))] <- paste0(names(x21)[grep("area", names(x21))], "_21")
  xall <- Reduce(function(...) merge(..., by = "PLAN_NO", all = TRUE), list(x17, x18, x19, x21))
  xall[is.na(xall)] <- 0
  xmean <- data.frame(PLAN_NO = xall$PLAN_NO, mean = rowMeans(xall[,grepl("area_msq", names(xall))]), area17 = xall$area_msq_17, area18 = xall$area_msq_18, area19 = xall$area_msq_19, area21 = xall$area_msq_21)
  names(xmean) <- c("PLAN_NO", hab_name, paste0(hab_name, "17"), paste0(hab_name, "18"), paste0(hab_name, "19"), paste0(hab_name, "21"))
  return(xmean)
}


#' Arable
#' 
arable <- lcm_extract(3, "Arable")
#'
#' Broadleaved woodland
#' 
broadleaved <- lcm_extract(1, "Broadleaved") 
#'
#' Coniferous woodland
#' 
coniferous <- lcm_extract(2, "Coniferous")
#' 
#' Improved grass
#' 
imp_grass <- lcm_extract(4, "ImpGrass")
#' 
#' Semi-natural grassland - removed heather and heather grassland
#' 
sn_grass <- lcm_extract(c(6,5,8), "SNGrass")
#'
#' Mountain bog heath
#' 
mbh <- lcm_extract(c(9,10,11,12), "MBH")
#' 
#' Coastal
#' 
coastal <- lcm_extract(c(15,16,17,18,19), "Coastal")
#' 
#' 
#' Calcium carbonate content
#' 
calc_carbonate_lookup <- data.frame(Cat = unique(soils$CARB_CNTNT), Conc = c(0,NA,0,1,0.5,0.5,1,1,1))

soils$calc_carbonate <- calc_carbonate_lookup$Conc[match(soils$CARB_CNTNT, calc_carbonate_lookup$Cat)]
#' 
#' Grain size - use raw median grain size
#' 
grain_lookup <- data.frame(Cat = unique(soils$PMM_GRAIN), Size = c(1, 2, 1, 0, 0.6, NA, 4, 3, 0.12, 0.03, 4, 0.6, 3, 0.03))

soils$grain_size <- grain_lookup$Size[match(soils$PMM_GRAIN, grain_lookup$Cat)]


soil_vars <- soils[,c(12:14)]

#join all datasets

all_vars <- Reduce(function(...) merge(..., by = "PLAN_NO"), 
                   list(NCAs, mean_elev[,2:3], SD_elev[,2:3], mean_rain[,2:3], SD_rain[,2:3], mean_temp[,2:3], 
                        SD_temp[,2:3],  slope[,2:3], aspect_south[,2:3], aspect_east[,2:3], soil_vars))

all_vars <- all_vars[complete.cases(all_vars),]

all_vars2 <- Reduce(function(...) merge(..., by = "PLAN_NO", all.x = TRUE), list(all_vars, SDA[,c(3,4)], hedges[,2:3], arable, broadleaved, coniferous, imp_grass, sn_grass, calcneut_grass, fenmarshswamp, heath, mbh, coastal, OFBP, spring_cereals, winter_cereals, maize, broadleaf_crops))

names(all_vars2)[1:18] <- c("PLAN_NO", "X", "NCA", "NCA_area", "NCA_perc", "Elevation", "SDElevation", "Rain", "SDRain", "Temp", "SDTemp", "Slope", "Southness", "Eastness", "CalcCarb", "GrainSize", "SDA", "Hedges")

all_vars2[is.na(all_vars2)] <- 0

# #write out table
# covariates <- all_vars2


#summarise to NCA level for ordination

NCA_vars <- aggregate(cbind(Elevation, SDElevation, Rain, SDRain, Temp, SDTemp, Slope, Southness, Eastness, Hedges, GrainSize, CalcCarb, SDA, Arable, Broadleaved, Coniferous, ImpGrass, SNGrass, CalcNeutGrass, HeathGrass, FenMarshSwamp, MBH, Coastal, OFBP, SpringCereals, WinterCereals, Maize, BroadleafCrops) ~ NCA, data = all_vars2, FUN = function(x) mean(x, na.rm = TRUE))

#add min and max elevation

NCA_vars2 <- aggregate(Elevation ~ NCA, data = all_vars2, FUN = function(x) min(x))
NCA_vars3 <- aggregate(Elevation ~ NCA, data = all_vars2, FUN = function(x) max(x))

names(NCA_vars2)[2] <- "MinElevation"
names(NCA_vars3)[2] <- "MaxElevation"

# #add habitat variables as proportion of NCA
# 
# NCA_vars4 <- aggregate(cbind() ~ NCA, data = all_vars2, FUN = function(x) sum(x, na.rm = TRUE))
# 
# NCA_squares <- aggregate(PLAN_NO ~ NCA, data = all_vars2, FUN = function(x) length(x))
# 
# NCA_vars4[,2:17] <- NCA_vars4[,2:18]/(NCA_squares$PLAN_NO*1000000)
# 

#combine
NCA_vars <- Reduce(function(...) merge(..., by = "NCA"), list(NCA_vars, NCA_vars2, NCA_vars3))

#rename NCAs to match NCA level variables

NCA_vars$NCA[c(2,8,90,133)] <- c("Avon Vale", "Blackmoor Vale and the Vale of Wardour", "North Yorkshire Moors and Cleveland Hills", "Breckland")


#match to NCA level variables

species_pools$NCA <- gsub("\\s*\\NCA\\s*\\([^\\)]+\\)", "", species_pools$NCA)

NCA_size$NCA <- gsub("\\s*\\NCA\\s*\\([^\\)]+\\)", "", NCA_size$M_TITLE)

NCA_vars$NCASize <- NCA_size$Area_km2[match(NCA_vars$NCA, NCA_size$NCA)]

NCA_vars$BirdPool <- species_pools$Mean_PriorityBirds[match(NCA_vars$NCA, species_pools$NCA)]

NCA_vars$ButterflyPool <- species_pools$Butterflies[match(NCA_vars$NCA, species_pools$NCA)]
  
pca1 <- princomp(scale(NCA_vars[,2:34]))


# Find max and min of first 3 PCs (for setting axis ranges)

png("Ordination with LM0465 NCAs v4.png")

plot(pca1$scores[,1], pca1$scores[,2], asp = 1, pch = 19, col = "grey", xlab= "PCA Axis 1", ylab = "PCA Axis 2", las = 1,
     xlim = c(-10,10), ylim = c(-5,10))

# Find our 6 NCAs (names, row indices and define abbreviations)
ncanames <- gsub(" NCA .*","", NCA_vars$NCA[grep("Dunsmore|Fens|x Clayland|High W|Dartm|Yorkshire D" ,  NCA_vars$NCA)])
ncaind2 <- grep(paste(ncanames, collapse = "|"), NCA_vars$NCA)
ncaadd2 <- c("DM","DF","HW","SS", "TF", "YD")

points(pca1$scores[ncaind2,1], pca1$scores[ncaind2,2], col = "blue", cex = 1.2, pch = 19)
text(pca1$scores[ncaind2,1], pca1$scores[ncaind2,2] , ncaadd2,  pos = c(1,3,1,3,1,3))

dev.off()

#biplot

#axes 1 and 2

png("Ordination biplot 1 v4 selected variables.png", height = 120, width = 120, units = "mm", res = 300, pointsize = 8)

plot(pca1$scores[,1], pca1$scores[,2], pch = 20, col = "grey65", xlab = "PCA Axis 1 (33%)", ylab = "PCA Axis 2 (12%)", xlim = c(-12,10), ylim = c(-10,12))

scale <- 1

lam <- (pca1$sdev[1:2]*sqrt(pca1$n.obs))^scale

len <- t(t(pca1$loadings[, 1:2]) * lam)*0.7


toplot <- c(22,13,1,2,3,4,17,15,10,33,5,14,26,28)
#toplot <- 1:nrow(len)

#row.names(len)[toplot] <- c("MeanAlt", "MeanPrec","MeanTemp","Moor and bog","Coastal","Birds","Broadleaf", "Arable","Area","Butterflies")

mapply(function(x,y) arrows(0, 0, x, y, col = "blue", length = .1),
       len[toplot,1], len[toplot,2])

textpos <- t(t(pca1$loadings[, 1:2]) * lam)*0.8

text(textpos[toplot,1], textpos[toplot,2], labels = row.names(len)[toplot], font = 2)

dev.off()

#axes 1 and 3

png("Ordination biplot 2 v4 selected variables.png", height = 120, width = 120, units = "mm", res = 300, pointsize = 8)
plot(pca1$scores[,1], pca1$scores[,3], pch = 20, col = "grey65", xlab = "PCA Axis 1 (33%)", ylab = "PCA Axis 3 (8%)", ylim = c(-10, 5), xlim = c(-10, 10))

scale <- 0.95

lam <- (pca1$sdev[c(1,3)]*sqrt(pca1$n.obs))^scale

len <- t(t(pca1$loadings[, c(1,3)]) * lam)*0.5


toplot <- c(23,32,5,14,26,12,10,29,1,2,3,18)
#toplot <- 1:nrow(len)

#row.names(len)[toplot] <- c("MeanAlt", "MeanPrec","MeanTemp","Moor and bog","Coastal","Birds","Broadleaf", "Arable","Area","Butterflies")

mapply(function(x,y) arrows(0, 0, x, y, col = "blue", length = .1),
       len[toplot,1], len[toplot,2])

textpos <- t(t(pca1$loadings[,c(1,3)]) * lam)*0.6

text(textpos[toplot,1], textpos[toplot,2], labels = row.names(len)[toplot], font = 2)

dev.off()


# Ordinations are not particularly useful - should also look at correlations at 1km level 

cor_mat <- round(cor(all_vars2[,c(6:19,23,27,31,35,39,43,47,51,55,59,63,67,71,75)], method = "spearman"),1)

ggcorrplot(cor_mat, type = "lower", hc.order = TRUE, lab = TRUE, lab_size = 3)

ggsave(height = 8, width = 8, units = "in", filename = "Correlation plot v4.png")

# NCA level correlations

cor_mat <- round(cor(NCA_vars[,c(2,3,4,6,7,11,9,13,14,15,16,18,19,20,22,23,24,26,28)], method = "spearman"),1)

ggcorrplot(cor_mat, type = "lower", hc.order = TRUE, lab = TRUE, lab_size = 3)

ggsave(height = 8, width = 8, units = "in", filename = "Correlation plot NCA v4.png")


# #ordination with subset of 19 selected variables
# 
# pca2 <- princomp(scale(NCA_vars[,c(2,3,4,6,7,11,9,13,14,15,16,18,19,20,22,23,24,26,28)]))
# 
# summary(pca2)
# 
# biplot(pca2, choices = c(1,3))
# 
# load_pca2 <- pca2$loadings[]
# 
# write.csv(load_pca2, "PCA 2 loadings.csv")
# 
# vars <- vector()
# for(i in 1:ncol(load_pca1)){
#   vars[i] <- row.names(load_pca1)[which.max(abs(load_pca1[,i]))]
# }
# 
# pca_importance <- function(x) {
#   vars <- x$sdev^2
#   vars <- vars/sum(vars)
#   rbind(`Standard deviation` = x$sdev, `Proportion of Variance` = vars, 
#         `Cumulative Proportion` = cumsum(vars))
# }
# 
# write.csv(pca_importance(pca2), "PCA 2 summary.csv")
# 
# #calculate ranks
# 
# rank_loads <- apply(load_pca2, 2, function(x) rank(abs(x)))
# 
# #weight by variable importance
# 
# comp_weights <- pca_importance(pca2)[2,]
# 
# weight_ranks <- rank_loads * rep(comp_weights, each= nrow(rank_loads))
# 
# rank_sum <- rowSums(weight_ranks)
# 
# sort(rank_sum, decreasing = TRUE)
# 
# #write.csv(sort(rank_sum, decreasing = TRUE), "Weighted rank sums loadings.csv")
# 
# 

## Ordination at square level to enable use of PCA scores in models ##

#subset variables to those that can be used in square level ordination
sq_vars <- c(1,6,7,8,9,10,11,12,13,14,15,16,17,18,19,23,27,31,35,39,43,47,51,55,59,63,67,71,75)

all_vars3 <- all_vars2[,sq_vars]

#scale
all_vars3sc <- scale(all_vars3[,2:29])

#fit PCA
pca1km <- princomp(all_vars3sc)

summary(pca1km)

#extract PCA scores per square
all_scores <- cbind(all_vars3$PLAN_NO,pca1km$scores)

save(all_scores, file = "Covariate selection/PCA scores per 1km square v2.Rdata")

#extract PCA loadings
pca1km_load <- pca1km$loadings[]

write.csv(pca1km_load, "Covariate selection/PCA 1km loadings v2.csv")

#extract axis variance
pca_importance <- function(x) {
  vars <- x$sdev^2
  vars <- vars/sum(vars)
  rbind(`Standard deviation` = x$sdev, `Proportion of Variance` = vars, 
        `Cumulative Proportion` = cumsum(vars))
}

write.csv(pca_importance(pca1km), "Covariate selection/PCA 1km summary v2.csv")


png("Ordination biplot 1km v2.png", height = 120, width = 120, units = "mm", res = 300, pointsize = 8)

plot(pca1km$scores[,1], pca1km$scores[,2], pch = 20, col = "grey65", xlab = "PCA Axis 1 (26%)", ylab = "PCA Axis 2 (12%)", xlim = c(-25,25), ylim = c(-25,25))

scale <- 1

lam <- (pca1km$sdev[1:2]*sqrt(pca1km$n.obs))^scale

len <- t(t(pca1km$loadings[, 1:2]) * lam)*0.05


toplot <- c(22,13,1,2,3,4,17,15,12,5,14,28,6)
#toplot <- 1:nrow(len)

#row.names(len)[toplot] <- c("MeanAlt", "MeanPrec","MeanTemp","Moor and bog","Coastal","Birds","Broadleaf", "Arable","Area","Butterflies")

mapply(function(x,y) arrows(0, 0, x, y, col = "blue", length = .1),
       len[toplot,1], len[toplot,2])

textpos <- t(t(pca1km$loadings[, 1:2]) * lam)*0.06

text(textpos[toplot,1], textpos[toplot,2], labels = row.names(len)[toplot], font = 2)

dev.off()


### correlations with AES gradients at national scale


fpath <- config_path

#' AES scores
AES <- read.csv(paste0(fpath, "Data/AES uptake/Outputs/All_1kmCells_Scored.csv"))

#' Get AES scores for WCBS
vars_aes <- AES %>%
  filter(CELLCODE %in% all_vars3$PLAN_NO) %>%
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
  tidyr::pivot_wider(names_from = c(scale, YEAR),
                     values_from = value) 

vars_aes$AES1KM <- rowMeans(vars_aes[,2:5])
vars_aes$AES3KM <- rowMeans(vars_aes[,6:9])

all_vars4 <- merge(all_vars3, vars_aes, by.x = "PLAN_NO", by.y = "CELLCODE")

cor(all_vars4[,c(13:29, 38, 39)], use = "pairwise.complete.obs", method = "spearman")
