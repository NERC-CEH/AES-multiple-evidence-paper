#' # Citizen Science data scoping - UK Butterfly Monitoring Scheme
#' 
#' *v3 19/11/2020*
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
fpath <- config_path

#' The datasets required for this scoping are as follows:
#' 
#' 1. UKBMS square locations
#' 
wcbs_ukbms <- read.csv(paste0(fpath, "Citizen science datasets/UKBMS_WCBS/site_data_2017-19.txt"))
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
CSlocs <- read_excel(paste0(fpath, "Citizen science datasets/CitizenScienceGridRefs.xlsx"), col_types = c("text", "text", "text", "numeric", "numeric"))
#'
#' 3. AES scores
#' 
AES <- read.csv(paste0(fpath, "AES uptake/Outputs/All_1kmCells_Scored.csv"))
#'
#' 4. Habitat variables
#' 
load("Data/Habitat variables for scoping.Rdata")
#' 
#' 5. Severely Disadvantaged Areas
#' 
sda <- read.csv(paste0(fpath, "Landscape and Environmental Data/LFA/sda_1km.csv"))
#' 
#' 6. Map of England
#' 
ENG <- read_sf("Data/england.shp")


#'
#' ## Data processing
#' 
#' Before we can start the scoping we need to make sure we can appropriately link up all the relevant datasets by their locations. To help with this John has produced a new worksheet with cleaned grid refs for each citizen science scheme. We will work with Eastings and Northings for our spatial reference system. The England shapefile is also in this system.
#' 
#' 1. Match UKBMS locations to the cleaned locations from John 
#' 
ukbms$Ecent <- CSlocs$Ecent[match(ukbms$SITENO, CSlocs$`SITE CODE`)]
ukbms$Ncent <- CSlocs$Ncent[match(ukbms$SITENO, CSlocs$`SITE CODE`)]
ukbms$GRIDREF_1KM <- CSlocs$GRIDREF_1KM[match(ukbms$SITENO, CSlocs$`SITE CODE`)]

#' 2. Match UKBMS data to AES data - this requires information on both the location & year of survey
#' 
#' To achieve this we only need to restructure the AES data so that year is a variable rather than multiple columns
#' 

#' As suggested we will remove any squares from the AES data that have either >50% woodland or >30% urban or freshwater from both the all England set of squares and the scheme squares
#' 

AES <- AES[AES$Masked == "",]
#removes about 20,000 squares

AES_melt_1km <- melt(AES, id.vars = c("CELLCODE", "Ecent", "Ncent", "M_TITLE", "nca_num", "NBuff", "NonMasked"), measure.vars = c("Sc17.1km", "Sc18.1km", "Sc19.1km", "Sc20.1km"))

AES_melt_1km$year <- as.numeric(paste0("20", substr(AES_melt_1km$variable, 3,4)))
names(AES_melt_1km)[9] <- "AES1KM"

AES_melt_3km <- melt(AES, id.vars = c("CELLCODE", "Ecent", "Ncent", "M_TITLE", "nca_num", "NBuff", "NonMasked"), measure.vars = c("Sc17.3km", "Sc18.3km", "Sc19.3km", "Sc20.3km"))

AES_melt_3km$year <- as.numeric(paste0("20", substr(AES_melt_3km$variable, 3,4)))
names(AES_melt_3km)[9] <- "AES3KM"

#remove 2020 for comparison

AES_melt_1km <- AES_melt_1km[AES_melt_1km$year %in% c(2017, 2018, 2019),]
AES_melt_3km <- AES_melt_3km[AES_melt_3km$year %in% c(2017, 2018, 2019),]


#' We can now merge the dataframes by grid reference and year to match the AES scores to the UKBMS data

ukbms_df <- merge(ukbms, AES_melt_1km, by.x = c("GRIDREF_1KM", "YEAR", "Ecent", "Ncent"), by.y = c("CELLCODE", "year", "Ecent", "Ncent"))

ukbms_df <- merge(ukbms_df, AES_melt_3km[,c(1,10,9)], by.x = c("GRIDREF_1KM", "YEAR"), by.y = c("CELLCODE", "year"))

#' 595 UKBMS squares don't match to AES scores.
#' 

ukbms_sqs <- unique(ukbms$GRIDREF_1KM[!ukbms$GRIDREF_1KM %in% ukbms_df$GRIDREF_1KM])
length(ukbms_sqs)

#' 38 of these are coastal squares. This was also the cause of missing AES scores for WCBS squares. In addition one square has NA values for 3km AES (SZ1790, coastal square near Bournemouth)
#' 
#' The remainder of squares are missing AES scores because they are in masked areas (i.e. > 50% woodland or > 30% urban or freshwater). There 1030 squares remaining so masking removed ~ 1/3 squares. This is a larger proportion than for WCBS.

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

ukbms_df2[is.na(ukbms_df2$area_msq),c(25,26)] <- 0
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
#' Firstly look at individual gradients
#' 
#' Plot the distribution of AES scores across England. Very time consuming to plot all 500,000 so take random sample of same size as CS data.  Graphs are plotted on a log scale because the distribution of AES scores is very skewed with many low values and few high values. For ease of interpretation vertical lines are added at scores of 500 and 5000 (used in design of LM0465 to differentiate Low, Medium and High AES squares).
#' 
#' 1km:
#' 
all_col <- rgb(1,0,0,0.2)
CS_col <- rgb(0,0,1,0.2)

plot_dat <- data.frame(AES1KM = c(AES_melt_1km$AES1KM[sample(1:nrow(AES_melt_1km), 
                                                             size = nrow(ukbms_df2), 
                                                             replace = FALSE)], 
                                  ukbms_df2$AES1KM), 
                       Source = c(rep("Whole countryside",nrow(ukbms_df2)),
                                  rep("UKBMS",nrow(ukbms_df2))))

ggplot(plot_dat) + 
  geom_histogram(aes(x = AES1KM, colour = Source, fill = Source), binwidth = 0.2) +
  facet_wrap(~Source, nrow = 2) +
  scale_color_manual(values = c(all_col, CS_col)) +
  scale_fill_manual(values = c(all_col, CS_col)) +
  labs(x = "1km AES score distribution") + 
  scale_x_log10() +
  scale_x_log10(oob=scales::oob_squish_infinite)+ 
  geom_vline(xintercept = 500) +
  geom_vline(xintercept = 5000) 

#' 3km:

plot_dat3 <- data.frame(AES3KM = c(AES_melt_3km$AES3KM[sample(1:nrow(AES_melt_3km), 
                                                              size = nrow(ukbms_df2), 
                                                              replace = FALSE)], 
                                   ukbms_df2$AES3KM), 
                        Source = c(rep("Whole countryside",nrow(ukbms_df2)),
                                   rep("UKBMS",nrow(ukbms_df2))))

ggplot(plot_dat3) + 
  geom_histogram(aes(x = AES3KM, colour = Source, fill = Source), binwidth = 0.2) +
  facet_wrap(~Source, nrow = 2) +
  scale_color_manual(values = c(all_col, CS_col)) +
  scale_fill_manual(values = c(all_col, CS_col)) +
  labs(x = "3km AES score distribution") + 
  scale_x_log10() +
  scale_x_log10(oob=scales::oob_squish_infinite)+ 
  geom_vline(xintercept = 500) +
  geom_vline(xintercept = 5000) 


#' Looks like UKBMS squares are slightly skewed to higher local AES now that masked squares have been removed, potential for these squares to be placed in "higher quality" areas?
#' 
#' Can also compare ranges
#' 
range(AES_melt_1km$AES1KM); range(ukbms_df2$AES1KM)
range(AES_melt_3km$AES3KM, na.rm = TRUE); range(ukbms_df2$AES3KM, na.rm=TRUE)

#' Ranges are ok, as expected UKBMS doesn't cover the very high ends of the local and landscape AES gradients
#' 
#' Also need to split into high, medium, low categories for both scores to compare coverage of both gradients.  These categories are defined as used in the design of LM0465. Tables show the proportion of squares in each of the 9 categories of local & landscape AES.
#' 
#' 
#' Low = 0 - 500
#' 
#' Medium = 501 - 5000
#' 
#' High = 5001 +
#' 
#need to create a combined all AES table

AES_all <- merge(AES_melt_1km[,c(1,9,10)], 
                 AES_melt_3km[,c(1,9,10)], 
                 by = c("CELLCODE", "year"))
AES_all$Cat_1KM <- cut(AES_all$AES1KM, 
                       c(-1,500,5000,1e9), 
                       c("Low","Medium","High"))
AES_all$Cat_3KM <- cut(AES_all$AES3KM, 
                       c(-1,500,5000,1e9), 
                       c("Low","Medium","High"))

#repeat for UKBMS

ukbms_df2$Cat_1KM <- cut(ukbms_df2$AES1KM, 
                        c(-1,500,5000,1e9), 
                        c("Low","Medium","High"))
ukbms_df2$Cat_3KM <- cut(ukbms_df2$AES3KM, 
                        c(-1,500,5000,1e9), 
                        c("Low","Medium","High"))

#table of proportions
Comb <- rbind(select(AES_all, AES1KM = Cat_1KM, AES3KM = Cat_3KM) %>%
                mutate(Survey = "Whole countryside", value = 1),
              select(ukbms_df2, AES1KM = Cat_1KM, AES3KM = Cat_3KM, value = 1) %>%
                mutate(Survey = "UKBMS", value = 1)) %>%
  count(Survey, AES1KM, AES3KM) %>%
  group_by(Survey) %>%
  mutate(prop = round(n/sum(n),3))

#+ fig.width = 8, fig.height = 4  
ggplot(na.omit(Comb), aes(x = AES1KM, y = AES3KM)) +
  facet_wrap(~Survey) +
  geom_tile(aes(fill = prop)) +
  geom_text(aes(label = prop)) +
  scale_fill_distiller(direction = 1, name = "proportion") +
  coord_fixed()



#' Overall UKBMS squares seem to cover all categories, and even have a higher proportion of High_Low squares than the national average. There is a higher proportion of squares in the high local categories than nationally.
#' 


#' 
#' 3. Evaluate whether AES gradients are confounded with habitat variables
#' 
#' Within the UKBMS data we might find that AES scores are highly correlated with other landscape or habitat variables. To assess this we will calculate Spearman rank correlations and plot the correlations
#' 
m1 <- round(cor(ukbms_df2[,c(17:24)], 
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

pairs(ukbms_df2[,c(17:24)], upper.panel = panel.cor,
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

AES_all <- merge(AES_all, sda[,c(3,4,5)],
                 by.x = "CELLCODE",
                 by.y = "PLAN_NO", 
                 all.x = TRUE)

AES_all[is.na(AES_all$area_msq),c(7,8)] <- 0

AES_all$SDA_cat <- cut(AES_all$area_m_pc,
                       c(-1,50,101),
                       c("Lowland", "Upland"))

sda_dat <- data.frame(SDA_cat = unlist(list(AES_all$SDA_cat, ukbms_df2$SDA_cat)),
                      Dataset = c(rep("All 1km sqs", nrow(AES_all)), rep("UKBMS squares", nrow(ukbms_df2))))

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
locs <- st_as_sf(data.frame(ukbms_df2[,c(1:5)]), 
                 coords = c("Ecent", "Ncent"))
plot(st_geometry(locs), main = "UKBMS square locations", 
     pch = 20, cex = 0.7)
plot(ENG,  add = TRUE)

#' There is evidence of a bias towards the south and away from the north of England. There are large spatial gaps in Cumbria, Northumbria, Lincolnshire and other parts of the Midlands and North. 
#'
#'#' 
#' We can also look at whether there are any obvious patterns in the AES scores across the country. A simple way to do this is to look for correlations between the AES gradients and Eastings and Northings. This won't pick up more subtle regional patterns but would pick up if e.g. all the squares with high local AES were in the south.
#' 
pairs(ukbms_df2[,c(3,4,17,18)], upper.panel = panel.cor, labels = c("Easting", "Northing", "AES 1km", "AES 3km"))

#' Some indication that squares with very high local AES tend to be in the south, but overall no strong correlations.


#' 
#' 6. Summarise differences in protocols
#'
#' REQUIRES INPUT FROM MARC