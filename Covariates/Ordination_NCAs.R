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
  x17$YEAR <- 2017
  x18$YEAR <- 2018
  x19$YEAR <- 2019
  x21$YEAR <- 2021
  xall <- rbind(x17[,c("PLAN_NO", "YEAR", "area_msq")],x18[,c("PLAN_NO", "YEAR", "area_msq")],x19[,c("PLAN_NO", "YEAR", "area_msq")],x21[,c("PLAN_NO", "YEAR", "area_msq")])
  xall[is.na(xall)] <- 0
  #xmean <- data.frame(PLAN_NO = xall$PLAN_NO, mean = rowMeans(xall[,grepl("area_msq", names(xall))]), area17 = xall$area_msq_17, area18 = xall$area_msq_18, area19 = xall$area_msq_19, area21 = xall$area_msq_21)
  names(xall) <- c("PLAN_NO", "YEAR", hab_name)
  return(xall)
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

all_vars <- Reduce(function(...) merge(..., by = "PLAN_NO"), list(mean_elev[complete.cases(mean_elev),2:3], SD_elev[complete.cases(SD_elev),2:3],  slope[complete.cases(slope),2:3], aspect_south[complete.cases(aspect_south),2:3], aspect_east[complete.cases(aspect_east),2:3], soil_vars))

all_vars <- all_vars[complete.cases(all_vars),]

all_vars2 <- Reduce(function(...) merge(..., by = "PLAN_NO", all.x = TRUE), list(all_vars, hedges[,2:3])) 
  
#subset already to considered locations
all_locs <- merge(climate, all_vars2, all.x= TRUE, by = "PLAN_NO")

all_locs2 <- Reduce(function(...) merge(..., by = c("PLAN_NO", "YEAR"), all.x = TRUE), list(all_locs, arable, broadleaved, coniferous, imp_grass, sn_grass, mbh, coastal))        
             
#remove 2020 rows
all_locs2 <- all_locs2[all_locs2$YEAR != 2020,]
                                                                       

names(all_locs2)[11:18] <- c("Elevation", "SDElevation", "Slope", "Southness", "Eastness", "CalcCarb", "GrainSize", "Hedges")
all_locs2 <- all_locs2[,c(1,2,4:9,11:25)]

#replace NA with zero for hedge length and habitat coverage
all_locs2[,16:23][is.na(all_locs2[,16:23])] <- 0

#remove missing DTM and soil values
all_locs2 <- all_locs2[complete.cases(all_locs2),]

ord_locs <- cbind(data.frame(ID = paste(all_locs2$PLAN_NO, all_locs2$YEAR, sep = "_")),all_locs2[3:23])


## Ordination at square level to enable use of PCA scores in models ##

#extract axis variance
pca_importance <- function(x) {
  vars <- x$sdev^2
  vars <- vars/sum(vars)
  rbind(`Standard deviation` = x$sdev, `Proportion of Variance` = vars, 
        `Cumulative Proportion` = cumsum(vars))
}

#ordination 1 - climate

clim_vars <- ord_locs[,c(1:7)]

clim_pca <- prcomp(clim_vars[,2:7], scale = TRUE)

summary(clim_pca)
biplot(clim_pca)


#extract PCA scores per square
clim_scores <- cbind(clim_vars$ID, clim_pca$x)

#extract PCA loadings
clim_load <- clim_pca$rotation

pca_importance(clim_pca)


png("Climate PCA biplot.png", height = 120, width = 120, units = "mm", res = 300, pointsize = 8)

plot(clim_pca$x[,1], clim_pca$x[,2], pch = 20, col = "grey65", xlab = "PCA Axis 1 (61%)", ylab = "PCA Axis 2 (26%)", xlim = c(-15,10), ylim = c(-6,6))

scale <- 1
lam <- (clim_pca$sdev[1:2]*sqrt(nrow(clim_vars)))^scale
len <- t(t(clim_pca$rotation[, 1:2]) * lam)*0.05

toplot <- 1:nrow(len)

mapply(function(x,y) arrows(0, 0, x, y, col = "blue", length = .1),
       len[toplot,1], len[toplot,2])
textpos <- t(t(clim_pca$rotation[, 1:2]) * lam)*0.06
text(textpos[toplot,1], textpos[toplot,2], labels = row.names(len)[toplot], font = 2)

dev.off()


#ordination 2 - landscape

land_vars <- ord_locs[,c(1,8:14)]

land_pca <- prcomp(land_vars[,2:8], scale = TRUE)

summary(land_pca)
biplot(land_pca)


#extract PCA scores per square
land_scores <- cbind(land_vars$ID, land_pca$x)

#extract PCA loadings
land_load <- land_pca$rotation

pca_importance(land_pca)


png("Landscape PCA biplot.png", height = 120, width = 120, units = "mm", res = 300, pointsize = 8)

plot(land_pca$x[,1], land_pca$x[,2], pch = 20, col = "grey65", xlab = "PCA Axis 1 (30%)", ylab = "PCA Axis 2 (16%)", ylim = c(-4,6))

scale <- 1
lam <- (land_pca$sdev[1:2]*sqrt(nrow(land_vars)))^scale
len <- t(t(land_pca$rotation[, 1:2]) * lam)*0.05

toplot <- 1:nrow(len)

mapply(function(x,y) arrows(0, 0, x, y, col = "blue", length = .1),
       len[toplot,1], len[toplot,2])
textpos <- t(t(land_pca$rotation[, 1:2]) * lam)*0.06
text(textpos[toplot,1], textpos[toplot,2], labels = row.names(len)[toplot], font = 2)

dev.off()




#ordination 3 - habitat



######NEED TO UPDATE#######

###ADD COMPOSIITIONAL ANALYSIS####

hab_vars <- ord_locs[,c(1,8:14)]

hab_pca <- prcomp(hab_vars[,2:8], scale = TRUE)

summary(hab_pca)
biplot(hab_pca)


#extract PCA scores per square
hab_scores <- cbind(hab_vars$ID, hab_pca$x)

#extract PCA loadings
hab_load <- hab_pca$rotation

pca_importance(hab_pca)


png("Habitat PCA biplot.png", height = 120, width = 120, units = "mm", res = 300, pointsize = 8)

plot(hab_pca$x[,1], hab_pca$x[,2], pch = 20, col = "grey65", xlab = "PCA Axis 1 (30%)", ylab = "PCA Axis 2 (16%)", ylim = c(-4,6))

scale <- 1
lam <- (hab_pca$sdev[1:2]*sqrt(nrow(hab_vars)))^scale
len <- t(t(hab_pca$rotation[, 1:2]) * lam)*0.05

toplot <- 1:nrow(len)

mapply(function(x,y) arrows(0, 0, x, y, col = "blue", length = .1),
       len[toplot,1], len[toplot,2])
textpos <- t(t(hab_pca$rotation[, 1:2]) * lam)*0.06
text(textpos[toplot,1], textpos[toplot,2], labels = row.names(len)[toplot], font = 2)

dev.off()



### checks

#correlation of first 3 PCA axes from climate and land PCAs
cor(clim_pca$x[,1:3], land_pca$x[,1:3])
