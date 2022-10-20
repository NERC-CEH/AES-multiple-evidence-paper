#' # Citizen Science data scoping - UK Butterfly Monitoring Scheme
#' 
#' *v3 20/10/2022*
#' 
#' Susan Jarvis
#' 
#' 
#' There are six elements of the citizen science (CS) data scoping exercise for project 07678 Landscape AES National Extrapolation (LANE):
#' 
#' 1. Summarise quantity of data available
#' 
#' 2. Assess coverage of CS data along AES gradients
#' 
#' 3. Evaluate whether AES gradients are confounded with habitat variables
#' 
#' 4. Assess distribution in uplands vs lowlands
#' 
#' 5. Consider whether CS data are regionally biased
#' 
#' 6. Summarise differences in protocols
#' 
#' 
#' This script will cover elements 1-5, plus preliminary work on data collation and processing
#' 
#' 
#' ## Set up workspace
#' 
#' Load required packages
#' 
library(readxl)
library(sf)
library(reshape)
library(ggplot2)
library(dplyr)
theme_set(theme_bw())
#' 
#' ## Data collation
#' 
#' Firstly set path to project folder
#' 
dir <- config::get()
fpath <- dir$directories$ukbmsdata
scpath <- dir$directories$scoredata
envpath <- dir$directories$envdata

#' The datasets required for this scoping are as follows:
#' 
#' 1. UKBMS square locations
#' 
wcbs_ukbms <-  read.table(paste0(fpath, "site_date_2017-2021.txt"), sep = "\t", header = TRUE)
#'
#' Subset to UKBMS
#' 
ukbms <- wcbs_ukbms[wcbs_ukbms$SCHEME == "UKBMS",]
#'
#' About double the number of UKBMS locations than WCBS locations
#'
#' 2. Cleaned CS dataset grid references 
#' 
#' Need to set column types to correctly read in the site code column
#+ warning = FALSE, message = FALSE 
CSlocs <- read.csv(paste0(scpath, "Butts_Gradient_Scores.csv"))
#'
#' 3. AES scores
#' 
AES <- read.csv(paste0(scpath, "Butts_Gradient_Scores.csv"))
#'
#' 4. Habitat variables
#' 
load("Data/Habitat variables for scoping.Rdata")
#' 
#' 5. Severely Disadvantaged Areas
#' 
sda <- read.csv(paste0(envpath, "LFA/sda_1km.csv"))
#' 
#' 6. Map of England
#' 
ENG <- read_sf(paste0(envpath, "EnglandOutline/england.shp"))


#'
#' ## Data processing
#' 
#' Before we can start the scoping we need to make sure we can appropriately link up all the relevant datasets by their locations. To help with this John has produced a new worksheet with cleaned grid refs for each citizen science scheme. We will work with Eastings and Northings for our spatial reference system. The England shapefile is also in this system.
#' 
#' 1. Match UKBMS locations to the cleaned locations from John 
#' 
CSlocs <- CSlocs[!duplicated(CSlocs[,2:14]),2:14]#remove duplicate rows, some squares have more than one site number due to multiple transects in that location
#remove locations for masked squares
CSlocs <- CSlocs[CSlocs$MaskStatus == "OK",]#also includes squares with NA mask which are coastal

ukbms$GRIDREF_1KM <- CSlocs$buttsurv.GRIDREF_1km[match(ukbms$SITENO, CSlocs$buttsurv.SITENO)]

#' 2. Match UKBMS data to AES data - this requires information on both the location & year of survey
#' 
#' To achieve this we only need to restructure the AES data so that year is a variable rather than multiple columns
#' 

#' As suggested we will remove any squares from the AES data that have either >50% woodland or >30% urban or freshwater from both the all England set of squares and the scheme squares
#' 

AES <- AES[AES$MaskStatus == "OK" & !is.na(AES$MaskStatus),]
AES <- AES[!duplicated(AES[,2:14]),2:14]#remove duplicate rows, some squares have more than one site number due to multiple transects in that location

AES_melt_1km <- melt(AES, id.vars = c("CELLCODE", "buttsurv.SITENO", "buttsurv.GRIDREF_1km"), measure.vars = c("Score1km_2017", "Score1km_2018", "Score1km_2019", "Score1km_2021"))

AES_melt_1km$year <- as.numeric(substr(AES_melt_1km$variable, 10,13))
names(AES_melt_1km)[5] <- "AES1KM"

AES_melt_3km <- melt(AES, id.vars = c("CELLCODE", "buttsurv.SITENO", "buttsurv.GRIDREF_1km"), measure.vars = c("Score3km_2017", "Score3km_2018", "Score3km_2019", "Score3km_2021"))

AES_melt_3km$year <- as.numeric(substr(AES_melt_3km$variable, 10,13))
names(AES_melt_3km)[5] <- "AES3KM"

#remove 2020 for comparison

AES_melt_1km <- AES_melt_1km[AES_melt_1km$year %in% c(2017, 2018, 2019, 2021),]
AES_melt_3km <- AES_melt_3km[AES_melt_3km$year %in% c(2017, 2018, 2019, 2021),]


#' We can now merge the dataframes by grid reference and year to match the AES scores to the UKBMS data

ukbms_df <- merge(ukbms, AES_melt_1km, by.x = c("SITENO", "GRIDREF_1KM", "YEAR"), by.y = c("buttsurv.SITENO", "CELLCODE", "year"))

ukbms_df <- merge(ukbms_df, AES_melt_3km[,c(1,5,6)], by.x = c("GRIDREF_1KM", "YEAR"), by.y = c("CELLCODE", "year"))

#' 754 UKBMS squares don't match to AES scores.
#' 

ukbms_sqs <- unique(ukbms$SITENO[is.na(ukbms$GRIDREF_1KM)])
length(ukbms_sqs)



#' 
#' 3. Match the UKBMS data to the habitat data
#' 

ukbms_df2 <- merge(ukbms_df, habitat_vars, by.x = c("GRIDREF_1KM", "YEAR"), by.y = c("PLAN_NO", "year"))

#' All remaining squares have habitat information
#' 
#' 4. Match UKBMS data to SDA data
#' 
ukbms_df2 <- merge(ukbms_df2, sda[,c(3,4,5)],
                 by.x = "GRIDREF_1KM",
                 by.y = "PLAN_NO", 
                 all.x = TRUE)

ukbms_df2[is.na(ukbms_df2$area_msq),c(20,21)] <- 0
#'

#' ## Scoping
#' 
#' 1. Summarise quantity of data available
#' 
#' Total number of squares visited per year 

tapply(ukbms_df2$GRIDREF, ukbms_df2$YEAR, function(x) length(unique(x)))

#' Still have a reasonable amount of data after masking
#' 
#' Average number of repeat visits per year
#' 
tapply(ukbms_df2$N_VISITS_MAYTOAUGUST, ukbms_df2$YEAR, function (x) summary(x))
#' 
#' On average squares have around 14-15 visits. There is a big range, some squares only have 1 visit in a year while others have up to 52. Need to think about how to deal with this
#' 
#' 2. Assess coverage of CS data along AES gradients
#'
#' Don't have the full England coverage of AES scores for 2021

#' 
#' 3. Evaluate whether AES gradients are confounded with habitat variables
#' 
#' Within the UKBMS data we might find that AES scores are highly correlated with other landscape or habitat variables. To assess this we will calculate Spearman rank correlations and plot the correlations
#' 
m1 <- round(cor(ukbms_df2[,c(12:19)], 
                use = "pairwise.complete.obs", 
                method = "spearman"), 2)
m1[upper.tri(m1)] <- ""
m1 <- as.data.frame(m1)
m1

#' We can also plot these correlations to check for any non-linearities
#+ fig.width = 8, fig.height = 8  
panel.cor <- function(x, y, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y, use = "pairwise.complete.obs", method = "spearman"), digits=2)
  txt <- paste0("R = ", r)
  cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt)
}

pairs(ukbms_df2[,c(12:19)], upper.panel = panel.cor,
      col = "black", pch = 20)


#' There is isn't much evidence of relationships between AES gradients and habitat variables, slight positive correlations with semi-natural habitat type c
#' 
#' 
#' 4. Assess distribution in uplands vs lowlands
#' 
#' We will classify 1km squares as upland if more than 50% of their area is classified as Severely Disadvantaged Land. 
#' 
ukbms_df2$SDA_cat <- cut(ukbms_df2$area_m_pc, 
                       c(-1,50,101), 
                       c("Lowland", "Upland"))

#need to calculate SDA scores for all squares for comparison to whole countryside

habitat_vars$area_m_pc <- sda$area_m_pc[match(habitat_vars$PLAN_NO, sda$PLAN_NO)]

habitat_vars$area_m_pc[is.na(habitat_vars$area_m_pc)] <- 0
habitat_vars$SDA_cat <- cut(habitat_vars$area_m_pc,
                       c(-1,50,101),
                       c("Lowland", "Upland"))

sda_dat <- data.frame(SDA_cat = unlist(list(habitat_vars$SDA_cat, ukbms_df2$SDA_cat)),
                      Dataset = c(rep("All 1km sqs", nrow(habitat_vars)), rep("UKBMS squares", nrow(ukbms_df2))))

ggplot(sda_dat, aes(fill = SDA_cat, x = Dataset))+
  geom_bar(position = "fill", stat = "count")+
  scale_fill_brewer(palette="Paired")
#'
#' UKBMS squares have a lower proportion of upland squares than the national average suggesting the uplands are quite poorly represented in the UKBMS sample.
#' 
#' 5. Consider whether CS data are regionally biased
#' 
#' To do this we can map the locations to identify potential areas with low coverage
#+ warning = FALSE

library(rnrfa)

ukbms_df2$Ecent <- osg_parse(ukbms_df2$GRIDREF_1KM)$easting + 500
ukbms_df2$Ncent <- osg_parse(ukbms_df2$GRIDREF_1KM)$northing + 500

locs <- st_as_sf(data.frame(ukbms_df2), 
                 coords = c("Ecent", "Ncent"))
plot(st_geometry(locs), main = "UKBMS square locations", 
     pch = 20, cex = 0.7)
plot(ENG,  add = TRUE)

#' There is evidence of a bias towards the south and away from the north of England. There are large spatial gaps in Cumbria, Northumbria, Lincolnshire and other parts of the Midlands and North. 
#'
#'#' 
#' We can also look at whether there are any obvious patterns in the AES scores across the country. A simple way to do this is to look for correlations between the AES gradients and Eastings and Northings. This won't pick up more subtle regional patterns but would pick up if e.g. all the squares with high local AES were in the south.
#' 
pairs(ukbms_df2[,c(23,24,12,13)], upper.panel = panel.cor, labels = c("Easting", "Northing", "AES 1km", "AES 3km"))

#' Some indication that squares with very high local AES tend to be in the south, but overall no strong correlations.


#' 
#' 6. Summarise differences in protocols
#'
#' REQUIRES INPUT FROM MARC