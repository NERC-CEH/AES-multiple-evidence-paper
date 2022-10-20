#' # Habitat variable collation for scoping work
#' 
#' *v1 09/11/2020*
#' 
#' Susan Jarvis
#' 
#' For AES scoping work we require four habitat variables to be calculated from Land Cover Map for each year (2017-2019) for each 1km square:
#' 
#' 1. Area of arable
#' 
#' 2. Area of grassland
#' 
#' 3. Area of semi-natural habitat
#' 
#' 4. Diversity of aggregate classes
#' 
#'
#' ## Collate data
#' 
#' Load required packages
#' 
library(vegan)
library(reshape)
#' 
#' Set path to data directory
dir <- config::get()
fpath <- dir$directories$envdata
#' 
#' Read in LCM 1km summaries for each year
lcm17 <- read.csv(paste0(fpath,"LCM/lcm17_1km.csv"))

lcm18 <- read.csv(paste0(fpath,"LCM/lcm18_1km.csv"))

lcm19 <- read.csv(paste0(fpath,"LCM/lcm19_1km.csv"))

lcm21 <- read.csv(paste0(fpath,"LCM/lcm21_1km.csv"))

#' ## Extract each variable
#' 
#' Use area in msq rather than percentages as percentages are rounded slightly
#' 
#' ### 1 - Area of arable
#'

arable17 <- aggregate(area_msq ~ PLAN_NO, data = lcm17, subset = gridcode == 3, FUN = sum)
arable18 <- aggregate(area_msq ~ PLAN_NO, data = lcm18, subset = gridcode == 3, FUN = sum)
arable19 <- aggregate(area_msq ~ PLAN_NO, data = lcm19, subset = gridcode == 3, FUN = sum)
arable21 <- aggregate(area_msq ~ PLAN_NO, data = lcm21, subset = gridcode == 3, FUN = sum)

arable17$year <- 2017
arable18$year <- 2018
arable19$year <- 2019
arable21$year <- 2021

arableall <- rbind(arable17, arable18, arable19, arable21)


#' ### 2 - Area of grassland
#' 
#' Contains Improved grassland (4), neutral grassland (5), calcareous grassland (6), acid grassland (7) and fen, marsh, swamp (8). Excludes heather grassland (10)

grass17 <- aggregate(area_msq ~ PLAN_NO, data = lcm17, subset = gridcode %in% c(4,5,6,7,8), FUN = sum)
grass18 <- aggregate(area_msq ~ PLAN_NO, data = lcm18, subset = gridcode %in% c(4,5,6,7,8), FUN = sum)
grass19 <- aggregate(area_msq ~ PLAN_NO, data = lcm19, subset = gridcode %in% c(4,5,6,7,8), FUN = sum)
grass21 <- aggregate(area_msq ~ PLAN_NO, data = lcm21, subset = gridcode %in% c(4,5,6,7,8), FUN = sum)

grass17$year <- 2017
grass18$year <- 2018
grass19$year <- 2019
grass21$year <- 2021

grassall <- rbind(grass17, grass18, grass19, grass21)

#' ### 3 - Area of semi-natural habitat
#' 
#' There are three classifications of semi-natural habitat we will consider:
#' 
#' a) Total SNH - this includes acid grassland (7), bog (11), broadleaf woodland (1), calcareous grassland (6), coniferous woodland (2), fen marsh swamp (8), freshwater (14), heather (9), heather grassland (10), inland rock (12), neutral grassland (5) and saltmarsh (19)
#' 
snha17 <- aggregate(area_msq ~ PLAN_NO, data = lcm17, subset = gridcode %in% c(7,11,1,6,2,8,14,9,10,12,5,19), FUN = sum)
snha18 <- aggregate(area_msq ~ PLAN_NO, data = lcm18, subset = gridcode %in% c(7,11,1,6,2,8,14,9,10,12,5,19), FUN = sum)
snha19 <- aggregate(area_msq ~ PLAN_NO, data = lcm19, subset = gridcode %in% c(7,11,1,6,2,8,14,9,10,12,5,19), FUN = sum)
snha21 <- aggregate(area_msq ~ PLAN_NO, data = lcm21, subset = gridcode %in% c(7,11,1,6,2,8,14,9,10,12,5,19), FUN = sum)

snha17$year <- 2017
snha18$year <- 2018
snha19$year <- 2019
snha21$year <- 2021

snhaall <- rbind(snha17, snha18, snha19, snha21)

#' 
#' 
#' b) SNH minus acid grassland (7) - this includes all habitats in semi-natural grassland except acid grassland
#' 
snhb17 <- aggregate(area_msq ~ PLAN_NO, data = lcm17, subset = gridcode %in% c(11,1,6,2,8,14,9,10,12,5,19), FUN = sum)
snhb18 <- aggregate(area_msq ~ PLAN_NO, data = lcm18, subset = gridcode %in% c(11,1,6,2,8,14,9,10,12,5,19), FUN = sum)
snhb19 <- aggregate(area_msq ~ PLAN_NO, data = lcm19, subset = gridcode %in% c(11,1,6,2,8,14,9,10,12,5,19), FUN = sum)
snhb21 <- aggregate(area_msq ~ PLAN_NO, data = lcm21, subset = gridcode %in% c(11,1,6,2,8,14,9,10,12,5,19), FUN = sum)

snhb17$year <- 2017
snhb18$year <- 2018
snhb19$year <- 2019
snhb21$year <- 2021


snhball <- rbind(snhb17, snhb18, snhb19, snhb21)
#' 
#' c) semi-natural grassland - this includes neutral grassland (5), calcareous grassland (6), fen marsh swamp (8), heather grassland (10) and heather (9)
#' 
snhc17 <- aggregate(area_msq ~ PLAN_NO, data = lcm17, subset = gridcode %in% c(5,6,8,9,10), FUN = sum)
snhc18 <- aggregate(area_msq ~ PLAN_NO, data = lcm18, subset = gridcode %in% c(5,6,8,9,10), FUN = sum)
snhc19 <- aggregate(area_msq ~ PLAN_NO, data = lcm19, subset = gridcode %in% c(5,6,8,9,10), FUN = sum)
snhc21 <- aggregate(area_msq ~ PLAN_NO, data = lcm21, subset = gridcode %in% c(5,6,8,9,10), FUN = sum)

snhc17$year <- 2017
snhc18$year <- 2018
snhc19$year <- 2019
snhc21$year <- 2021

snhcall <- rbind(snhc17, snhc18, snhc19, snhc21)
#' 
#' 
#' ### 4 - Habitat diversity
#' 
#' Habitat diversity includes 9 of the 10 aggregate classes (minus saltwater (agg class 7))
#' 
#' First we need to create a lookup between target class and aggregate class
#' 
hablookup <- data.frame(target_class = 1:21, agg_class = c(1,2,3,4,5,5,5,5,6,6,6,6,7,8,8,9,9,9,9,10,10))

lcm17$agg_class <- hablookup$agg_class[match(lcm17$gridcode, hablookup$target_class)]
lcm18$agg_class <- hablookup$agg_class[match(lcm18$gridcode, hablookup$target_class)]
lcm19$agg_class <- hablookup$agg_class[match(lcm19$gridcode, hablookup$target_class)]
lcm21$agg_class <- hablookup$agg_class[match(lcm21$gridcode, hablookup$target_class)]

c17 <- cast(lcm17, PLAN_NO ~ agg_class, value = "area_msq", fun.aggregate = "sum", fill = 0, subset = lcm17$agg_class != 7)
c18 <- cast(lcm18, PLAN_NO ~ agg_class, value = "area_msq", fun.aggregate = "sum", fill = 0, subset = lcm18$agg_class != 7)
c19 <- cast(lcm19, PLAN_NO ~ agg_class, value = "area_msq", fun.aggregate = "sum", fill = 0, subset = lcm19$agg_class != 7)
c21 <- cast(lcm21, PLAN_NO ~ agg_class, value = "area_msq", fun.aggregate = "sum", fill = 0, subset = lcm21$agg_class != 7)


c17$div <- diversity(c17[,2:10], index = "shannon")
c18$div <- diversity(c18[,2:10], index = "shannon")
c19$div <- diversity(c19[,2:10], index = "shannon")
c21$div <- diversity(c21[,2:10], index = "shannon")

c17$year <- 2017
c18$year <- 2018
c19$year <- 2019
c21$year <- 2021

habdivall <- rbind(c17[,c(1,11,12)], c18[,c(1,11,12)], c19[,c(1,11,12)], c21[,c(1,11,12)])

#' Now we have all the variables we can combine all into a single dataframe with one entry per square/year combination
#' 

habitat_vars <- Reduce(function(...) merge(..., all = TRUE, by = c("PLAN_NO", "year")), list(arableall, grassall, snhaall, snhball, snhcall, habdivall))

names(habitat_vars) <- c("PLAN_NO", "year", "arable", "grass", "snh_a", "snh_b", "snh_c", "hab_div")

#' Check every 1km is included
#' 
length(unique(habitat_vars$PLAN_NO))
length(unique(lcm17$PLAN_NO))

#' Difference likely to be due to squares with all sea
#' 
nrow(lcm17[lcm17$agg_class == 7 & lcm17$area_m_pc > 99,])

#' accounts for about 1/3rd...
#' 
#' 
#' 

habitat_vars[,3:7] <- sapply(habitat_vars[,3:7], function(x) replace(x, is.na(x), 0))

#' 
save(habitat_vars, file = "Data/Habitat variables for scoping.RData")
