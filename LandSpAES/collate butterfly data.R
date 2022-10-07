#' exploratory analysis of 2018 models - butt transects
#' 
library(RODBC)

#try to connect to AES database

db <- config::get("AES")

con <- odbcConnect(db$DSN, db$uid, db$pwd, believeNRows = FALSE)

#extract relevant tables
#new views
buttdata <- sqlFetch(con, "VIEW_BUTTERFLY_DATA")
buttvariables <- sqlFetch(con, "VIEW_BUTTERFLY_VARIABLES_CT")
buttvisit <- sqlFetch(con, "VIEW_BUTTERFLY_VISITS")
buttvisit2 <- sqlFetch(con, "VIEW_BUTTERFLY_SEC_VISITS")
buttround <- sqlFetch(con, "TBL_ROUND")
butttransect <- sqlFetch(con, "TBL_TRANSECT_SECTION")
buttvisittbl <- sqlFetch(con, "TBL_BUTTERFLY_VISIT")
buttsquares <- sqlFetch(con, "TBL_SURVEY_SQUARE")
buttscores <- sqlFetch(con, "TBL_SQUARE_SCORES")
buttspecies <- sqlFetch(con, "TBL_BUTTERFLY_SPECIES")
close(con)


library(sqldf)

#' Now need to join butttransect to buttcount to add visits with no butts observed
#' 

buttvisit2$NEW_ID <- paste(buttvisit2$SURVEY_SQUARE, buttvisit2$SURVEY_YEAR, buttvisit2$ROUND_NUMBER, buttvisit2$TRANSECT_SECTION_NUMBER, sep = "_")

#buttdata$SURVEY_YEAR <- format(buttdata$SURVEY_DATE, "%Y")

buttdata$NEW_ID <- paste(buttdata$SURVEY_SQUARE, buttdata$SURVEY_YEAR, buttdata$ROUND_NUMBER, buttdata$TRANSECT_SECTION_NUMBER, sep = "_")


buttcount2 <- sqldf("select t1.NCA, t1.SURVEY_SQUARE, t1.SURVEY_YEAR, t1.ROUND_NUMBER, t1.TRANSECT_SECTION_NUMBER, t1.NEW_ID,t2.BUTTERFLY_NAME, t2.BUTTERFLY_SPECIES, t2.BUTTERFLY_COUNT, t2.CROP_PESTS, t2.BUTTERFLY_COUNT_ID, t2.BUTTERFLY_ID, t2.BUTTERFLY_VISIT_ID
                   from buttvisit2 as t1
                   left join buttdata as t2
                   on t1.NEW_ID = t2.NEW_ID")




#remove data from failed visits - visit ID 502 - clarify with Jo etc

failed_IDs <- buttvisittbl$BUTTERFLY_VISIT_ID[!is.na(buttvisittbl$FAILED_SURVEY)]
# 
buttcount2 <- buttcount2[!buttcount2$BUTTERFLY_VISIT_ID %in% failed_IDs,] 


#replace NA counts with 0
buttcount2$BUTTERFLY_COUNT[is.na(buttcount2$BUTTERFLY_COUNT_ID)] <- 0

#if butterfly species is entered but count is 0 change butterfly ID to NA
buttcount2$BUTTERFLY_ID[buttcount2$BUTTERFLY_COUNT == 0] <- NA


##add scores

buttscores$SURVEY_SQUARE <- buttsquares$SURVEY_SQUARE[match(buttscores$SURVEY_SQUARE_ID, buttsquares$SURVEY_SQUARE_ID)]

buttcount3 <- merge(buttcount2, buttscores, by.x = c("SURVEY_SQUARE", "SURVEY_YEAR"), by.y = c("SURVEY_SQUARE", "SURVEY_YEAR"))[,c(1:13,17,19)]

names(buttcount3)[14:15] <- c("AES1KM","AES3KM")


buttcount2 <- buttcount3


#' ##add floral resources
#' 
#' 
#' db <- config::get("AES")
#' con <- odbcConnect(db$DSN, db$uid, db$pwd, believeNRows = FALSE)
#' flowerindex <- sqlFetch(con, "TBL_FLOWER_INDEX") 
#' flowerspecies <- sqlFetch(con, "TBL_SPECIES")
#' flowerindexsurvey <- sqlFetch(con, "TBL_FLOWER_INDEX_SURVEY")
#' flowerdata <- sqlFetch(con, "VIEW_FLOWER_DATA")
#' flowervisits <- sqlFetch(con, "VIEW_FLOWER_VISITS")
#' close(con)
#' 
#' #' Need to create a flower index ID lookup with the midpoints of each class as suggested by Jo
#' #' Using the same range as Claire Carvell (although in Jo's email she suggests the ranges are different they seem to be the same in the protocol)
#' #'
#' flowerlookup <- data.frame(Index = c(0,1,2,3,4,5,6), Midpoint = c(0,3,15.5,113,600.5,3000.5,12000))
#' 
#' #' Add column to flower index with lookup midpoint value
#' #'
#' flowerdata$FLOWER_VALUE <- flowerlookup$Midpoint[match(flowerdata$FLOWER_INDEX,flowerlookup$Index)]
#' 
#' #' Can now sum by `NEW_ID` - this will sum across species and quadrats in one fell swoop
#' #'
#' 
#' flowerdata$SURVEY_YEAR <- format(flowerdata$SURVEY_DATE, "%Y")
#' 
#' flowerdata$NEW_ID <-  paste(flowerdata$SURVEY_SQUARE, flowerdata$SURVEY_YEAR, flowerdata$ROUND_NUMBER, flowerdata$TRANSECT_SECTION_NUMBER, sep = "_")
#' 
#' flowertransect <- aggregate(FLOWER_VALUE ~ NEW_ID, data = flowerdata, FUN = "sum")
#' 
#' #3161 obs
#' 
#' hist(flowertransect$FLOWER_VALUE)
#' 
#' 
#' ###join flower data based on NEW_ID
#' 
#' buttcount2$FLORAL_RESOURCES <- flowertransect$FLOWER_VALUE[match(buttcount2$NEW_ID, flowertransect$NEW_ID)]
#' 
#' ###177 transect sections have missing floral data - indicate no floral resources
#' 
#' 
#' #' Convert NA values to 0
#' #'
#' buttcount2$FLORAL_RESOURCES[is.na(buttcount2$FLORAL_RESOURCES)] <- 0
#' 
#' ##add floral richness
#' 
#' flowerrichtransect <- aggregate(SPECIES_ID ~ NEW_ID, data = flowerdata, FUN = function(x) length(unique(x[!is.na(x)])), na.action = "na.pass")
#' #1218 entries - same as flower transect
#' 
#' buttcount2$FLORAL_RESOURCES_RICH <- flowerrichtransect$SPECIES_ID[match(buttcount2$NEW_ID, flowerrichtransect$NEW_ID)]
#' 
#' buttcount2$FLORAL_RESOURCES_RICH[is.na(buttcount2$FLORAL_RESOURCES_RICH)] <- 0
#' 
#' # transect length
#' 
#' 
#' buttcount2$T_ID <- paste(buttcount2$SURVEY_SQUARE, buttcount2$TRANSECT_SECTION_NUMBER, sep = "_")
#' 
#' buttvariables$SURVEY_YEAR <- format(buttvariables$SURVEY_DATE, "%Y")
#' 
#' buttvariables$T_ID <- paste(buttvariables$SURVEY_SQUARE,  buttvariables$TRANSECT_SECTION_NUMBER, sep = "_")
#' 
#' buttlengths <- aggregate(SECTION_LENGTH ~ T_ID, data = buttvariables, FUN = function(x) max(x, na.rm=T))
#' 
#' buttcount2$TRANSECT_LENGTH <- buttlengths$SECTION_LENGTH[match(buttcount2$T_ID, buttlengths$T_ID)]

## allocate aggregates

## Correct issues with incorrect species name and aggregates

#allocate all to aggregate - change from 2018
buttcount2$BUTTERFLY_ID[buttcount2$BUTTERFLY_NAME == "Pearl-bordered Fritillary" & !is.na(buttcount2$BUTTERFLY_NAME)] <- 62

buttcount2$BUTTERFLY_ID[buttcount2$BUTTERFLY_NAME == "Small Pearl-bordered Fritillary" & !is.na(buttcount2$BUTTERFLY_NAME)] <- 62


#for small/essex skipper

sm_skip_sq <- aggregate(BUTTERFLY_COUNT ~ SURVEY_SQUARE + SURVEY_YEAR, data = buttcount2, FUN = function(x) sum(x), subset = BUTTERFLY_ID == 52)

es_skip_sq <- aggregate(BUTTERFLY_COUNT ~ SURVEY_SQUARE + SURVEY_YEAR, data = buttcount2, FUN = function(x) sum(x), subset = BUTTERFLY_ID == 14)

se_skip_sq <- aggregate(BUTTERFLY_COUNT ~ SURVEY_SQUARE + SURVEY_YEAR, data = buttcount2, FUN = function(x) sum(x), subset = BUTTERFLY_ID == 55)

se_skip_1 <- merge(sm_skip_sq, es_skip_sq, by.x = c("SURVEY_SQUARE", "SURVEY_YEAR"),by.y = c("SURVEY_SQUARE", "SURVEY_YEAR"), all = TRUE)

se_skip_2 <- merge(se_skip_1, se_skip_sq, by.x = c("SURVEY_SQUARE", "SURVEY_YEAR"),by.y = c("SURVEY_SQUARE", "SURVEY_YEAR"), all = TRUE)

names(se_skip_2) <- c("SURVEY_SQUARE", "SURVEY_YEAR", "SMALL_SKIP", "ESS_SKIP", "AGG")

#replace NA with 0
se_skip_2[,3:5][is.na(se_skip_2[,3:5])] <- 0

#add total recorded to spp
se_skip_2$tot <- rowSums(se_skip_2[,3:4])

#calc prop of each type
se_skip_2$P_sm <- se_skip_2$SMALL_SKIP/se_skip_2$tot
se_skip_2$P_es <- se_skip_2$ESS_SKIP/se_skip_2$tot




#for 3 squares with only aggregates recorded use average prop in NCA
se_skip_2$NCA <- buttcount2$NCA[match(se_skip_2$SURVEY_SQUARE, buttcount2$SURVEY_SQUARE)]

tapply(se_skip_2$P_sm, list(se_skip_2$NCA, se_skip_2$SURVEY_YEAR), FUN = function(x) mean(x, na.rm= T))

se_skip_2$P_sm[is.na(se_skip_2$P_sm) & se_skip_2$NCA  == "Dunsmore and Feldon" & se_skip_2$SURVEY_YEAR == 2017] <- 0.01388889

se_skip_2$P_sm[is.na(se_skip_2$P_sm) & se_skip_2$NCA  == "Dunsmore and Feldon" & se_skip_2$SURVEY_YEAR == 2019] <- 0.455803571

se_skip_2$P_sm[is.na(se_skip_2$P_sm) & se_skip_2$NCA == "High Weald" & se_skip_2$SURVEY_YEAR == 2017] <- 0.74390969

se_skip_2$P_es[is.na(se_skip_2$P_es)] <- 1-se_skip_2$P_sm[is.na(se_skip_2$P_es)] 


#allocate aggregate based on prop
se_skip_2$New_sm <- round(se_skip_2$AGG*se_skip_2$P_sm)
se_skip_2$New_ess <- round(se_skip_2$AGG*se_skip_2$P_es)

#change by hand
se_skip_2$New_ess[28] <- 5

#edit to ensure total new equal total aggregates
se_skip_2$tot_new <- rowSums(se_skip_2[,10:11])
se_skip_2$check <- se_skip_2$tot_new == se_skip_2$AGG


se_skip_3 <- se_skip_2[se_skip_2$tot_new > 0,c(1,2,5,10,11,12)]

##add back to buttcount2

aggrecs <- buttcount2[buttcount2$BUTTERFLY_ID == 55 & !is.na(buttcount2$BUTTERFLY_ID),c(2,1,8,9,12)]


aggrecs$BUTTERFLY_ID <- c(52,14,52,14,52,14,14,14,14,52,14,14,14,14,14,52,52,52,14,14,14,14,14,52,14,52,14,14,14,52,14,52,14,52,14,14,14,52,14,14,14,52,52)

buttcount2$BUTTERFLY_ID[buttcount2$BUTTERFLY_ID == 55 & !is.na(buttcount2$BUTTERFLY_ID)] <- aggrecs$BUTTERFLY_ID


###Green-veined/small whites

gvw_sq <- aggregate(BUTTERFLY_COUNT ~ SURVEY_SQUARE + SURVEY_YEAR + ROUND_NUMBER, data = buttcount2, FUN = function(x) sum(x), subset = BUTTERFLY_ID == 19)

sw_sq <- aggregate(BUTTERFLY_COUNT ~ SURVEY_SQUARE + SURVEY_YEAR + ROUND_NUMBER, data = buttcount2, FUN = function(x) sum(x), subset = BUTTERFLY_ID == 54)

gvs_sq <- aggregate(BUTTERFLY_COUNT ~ SURVEY_SQUARE + SURVEY_YEAR + ROUND_NUMBER, data = buttcount2, FUN = function(x) sum(x), subset = BUTTERFLY_ID == 64)

gvs_1 <- merge(gvw_sq, sw_sq, by.x = c("SURVEY_SQUARE", "SURVEY_YEAR", "ROUND_NUMBER"),by.y = c("SURVEY_SQUARE", "SURVEY_YEAR","ROUND_NUMBER"), all = TRUE)

gvs_2 <- merge(gvs_1, gvs_sq, by.x = c("SURVEY_SQUARE", "SURVEY_YEAR","ROUND_NUMBER"),by.y = c("SURVEY_SQUARE", "SURVEY_YEAR","ROUND_NUMBER"), all = TRUE)

names(gvs_2) <- c("SURVEY_SQUARE", "SURVEY_YEAR","ROUND_NUMBER", "GV_WHITE", "S_WHITE", "AGG")

#replace NA with 0
gvs_2[,4:6][is.na(gvs_2[,4:6])] <- 0

#add total recorded to spp
gvs_2$tot <- rowSums(gvs_2[,4:5])

#calc prop of each type
gvs_2$P_gv <- gvs_2$GV_WHITE/gvs_2$tot
gvs_2$P_sm <- gvs_2$S_WHITE/gvs_2$tot


#for squares with only aggregates recorded use average prop in NCA per year per round
gvs_2$NCA <- buttcount2$NCA[match(gvs_2$SURVEY_SQUARE, buttcount2$SURVEY_SQUARE)]

tapply(gvs_2$P_sm, list(gvs_2$NCA, gvs_2$ROUND_NUMBER, gvs_2$SURVEY_YEAR), FUN = function(x) mean(x, na.rm= T))

gvs_2$P_sm[is.na(gvs_2$P_sm) & gvs_2$NCA == "Yorkshire Dales" & gvs_2$SURVEY_YEAR == 2018 & gvs_2$ROUND_NUMBER == 1] <- 0.06666667
gvs_2$P_sm[is.na(gvs_2$P_sm) & gvs_2$NCA == "Yorkshire Dales" & gvs_2$SURVEY_YEAR == 2018 & gvs_2$ROUND_NUMBER == 2] <- 0.1000000
gvs_2$P_sm[is.na(gvs_2$P_sm) & gvs_2$NCA == "Yorkshire Dales" & gvs_2$SURVEY_YEAR == 2018 & gvs_2$ROUND_NUMBER == 3] <- 0.50416667
gvs_2$P_sm[is.na(gvs_2$P_sm) & gvs_2$NCA == "Yorkshire Dales" & gvs_2$SURVEY_YEAR == 2021 & gvs_2$ROUND_NUMBER == 3] <- 0.0000000
gvs_2$P_sm[is.na(gvs_2$P_sm) & gvs_2$NCA == "Dunsmore and Feldon" & gvs_2$SURVEY_YEAR == 2019 & gvs_2$ROUND_NUMBER == 1] <- 0.6363636
gvs_2$P_sm[is.na(gvs_2$P_sm) & gvs_2$NCA == "Dunsmore and Feldon" & gvs_2$SURVEY_YEAR == 2019 & gvs_2$ROUND_NUMBER == 2] <- 0.5862614
gvs_2$P_sm[is.na(gvs_2$P_sm) & gvs_2$NCA == "Dunsmore and Feldon" & gvs_2$SURVEY_YEAR == 2021 & gvs_2$ROUND_NUMBER == 3] <- 0.4600258
gvs_2$P_sm[is.na(gvs_2$P_sm) & gvs_2$NCA == "Dartmoor" & gvs_2$SURVEY_YEAR == 2019 & gvs_2$ROUND_NUMBER == 3] <- 0.0000000

gvs_2$P_gv[is.na(gvs_2$P_gv)] <- 1-gvs_2$P_sm[is.na(gvs_2$P_gv)] 


#allocate aggregate based on prop
gvs_2$New_gv <- round(gvs_2$AGG*gvs_2$P_gv)
gvs_2$New_sm <- round(gvs_2$AGG*gvs_2$P_sm)

#edit to ensure total new equal total aggregates
gvs_2$tot_new <- rowSums(gvs_2[,11:12])
gvs_2$check <- gvs_2$tot_new == gvs_2$AGG

#correct by hand
gvs_2$New_gv[c(3,100)] <- c(5,1)
gvs_2$New_sm[c(111,132)] <- c(1,3)
  
gvs_3 <- gvs_2[gvs_2$AGG > 0,c(1,2,3,10,11,12,13,14)]

##add back to buttcount2

aggrecs <- buttcount2[buttcount2$BUTTERFLY_ID == 64 & !is.na(buttcount2$BUTTERFLY_ID),c(2,1,4,7,9,12)]

for (i in 1:nrow(gvs_3)){
  rows <- aggrecs[aggrecs$SURVEY_YEAR == gvs_3[i,]$SURVEY_YEAR & aggrecs$SURVEY_SQUARE == gvs_3[i,]$SURVEY_SQUARE & aggrecs$ROUND_NUMBER == gvs_3[i,]$ROUND_NUMBER,]
  tot_count <- 0
  for(j in 1:nrow(rows)){
    if(gvs_3[i,"New_gv"]>tot_count){rows$BUTTERFLY_ID[j] <- 19} else {rows$BUTTERFLY_ID[j] <- 54}
    tot_count <- sum(rows[1:j, "BUTTERFLY_COUNT"])
  }
  aggrecs$BUTTERFLY_ID[aggrecs$SURVEY_YEAR == gvs_3[i,]$SURVEY_YEAR & aggrecs$SURVEY_SQUARE == gvs_3[i,]$SURVEY_SQUARE & aggrecs$ROUND_NUMBER == gvs_3[i,]$ROUND_NUMBER] <- rows$BUTTERFLY_ID
}

buttcount2$BUTTERFLY_ID[buttcount2$BUTTERFLY_ID == 64 & !is.na(buttcount2$BUTTERFLY_ID)] <- aggrecs$BUTTERFLY_ID
