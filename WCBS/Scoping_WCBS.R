#' # Citizen Science data scoping - Wider Countryside Butterfly Survey
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
#' 1. WCBS square locations
#' 
ukbms_wcbs <- read.csv(paste0(fpath, "Citizen science datasets/UKBMS_WCBS/site_data_2017-19.txt"))
#'
#' Subset to WCBS
#' 
wcbs <- ukbms_wcbs[ukbms_wcbs$SCHEME == "WCBS",]
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
#' 1. Match WCBS locations to the cleaned locations from John 
#' 
wcbs$Ecent <- CSlocs$Ecent[match(wcbs$GRIDREF, CSlocs$GRIDREF_1KM)]
wcbs$Ncent <- CSlocs$Ncent[match(wcbs$GRIDREF, CSlocs$GRIDREF_1KM)]

#' 2. Match WCBS data to AES data - this requires information on both the location & year of survey
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


#' We can now merge the dataframes by grid reference and year to match the AES scores to the WCBS data

wcbs_df <- merge(wcbs, AES_melt_1km, by.x = c("GRIDREF", "YEAR", "Ecent", "Ncent"), by.y = c("CELLCODE", "year", "Ecent", "Ncent"))

wcbs_df <- merge(wcbs_df, AES_melt_3km[,c(1,10,9)], by.x = c("GRIDREF", "YEAR"), by.y = c("CELLCODE", "year"))

#' Now that we have masked out squares with > 50% woodland or > 30% urban or freshwater we have fewer WCBS squares remaining
#' 

wcbs_sqs <- unique(wcbs$GRIDREF[!wcbs$GRIDREF %in% wcbs_df$GRIDREF])
wcbs_sqs

#' Without masking we had only 8 squares which didn't have associated AES scores (because of their coastal locations). After masking we now have 179 squares without AES scores and 723 remaining which do have scores.
#' 

#' 3. Match the WCBS data to the habitat data
#' 

wcbs_df2 <- merge(wcbs_df, habitat_vars, by.x = c("GRIDREF", "YEAR"), by.y = c("PLAN_NO", "year"))

#' All remaining squares have habitat information
#' 
#' 4. Match WCBS data to SDA data
#' 
wcbs_df2 <- merge(wcbs_df2, sda[,c(3,4,5)],
                 by.x = "GRIDREF",
                 by.y = "PLAN_NO", 
                 all.x = TRUE)

wcbs_df2[is.na(wcbs_df2$area_msq),c(24,25)] <- 0
#'

#' ## Scoping
#' 
#' 1. Summarise quantity of data available
#' 
#' Total number of squares visited per year 

tapply(wcbs_df2$GRIDREF, wcbs_df2$YEAR, function(x) length(unique(x)))

#' Still have a reasonable amount of data after masking
#' 
#' Average number of repeat visits per year
#' 
tapply(wcbs_df2$N_VISITS_MAYTOAUGUST, wcbs_df2$YEAR, function (x) summary(x))
#' 
#' On average over 2 visits per square per year, with some squares getting up to 8 visits. Most squares seem to get 2 visits
#' 
#' 2. Assess coverage of CS data along AES gradients
#'
#' Firstly look at individual gradients
#' 
#' Plot the distribution of AES scores across England. Very time consuming to plot all 500,000 so take random sample of same size as CS data. Graphs are plotted on a log scale because the distribution of AES scores is very skewed with many low values and few high values. For ease of interpretation vertical lines are added at scores of 500 and 5000 (used in design of LM0465 to differentiate Low, Medium and High AES squares)
#' 
#' 1km:
#' 
all_col <- rgb(1,0,0,0.2)
CS_col <- rgb(0,0,1,0.2)

plot_dat <- data.frame(AES1KM = c(AES_melt_1km$AES1KM[sample(1:nrow(AES_melt_1km), 
                              size = nrow(wcbs_df2), 
                              replace = FALSE)], 
                              wcbs_df2$AES1KM), 
                       Source = c(rep("Whole countryside",nrow(wcbs_df2)),
                                  rep("WCBS",nrow(wcbs_df2))))

ggplot(plot_dat) + 
  geom_histogram(aes(x = AES1KM, colour = Source, fill = Source), binwidth = 0.2) +
  facet_wrap(~Source, nrow = 2) +
  scale_color_manual(values = c(all_col, CS_col)) +
  scale_fill_manual(values = c(all_col, CS_col)) +
  labs(x = "1km AES score distribution")+ 
  scale_x_log10() +
  scale_x_log10(oob=scales::oob_squish_infinite)+ 
  geom_vline(xintercept = 500) +
  geom_vline(xintercept = 5000) 

#' 3km:

plot_dat3 <- data.frame(AES3KM = c(AES_melt_3km$AES3KM[sample(1:nrow(AES_melt_3km), 
                               size = nrow(wcbs_df2), 
                               replace = FALSE)], 
                               wcbs_df2$AES3KM), 
                        Source = c(rep("Whole countryside",nrow(wcbs_df2)),
                                   rep("WCBS",nrow(wcbs_df2))))

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


#' Good similarity between distribution of AES scores in WCBS squares and nationally
#' 
#' Can also compare ranges
#' 
range(AES_melt_1km$AES1KM); range(wcbs_df2$AES1KM)
range(AES_melt_3km$AES3KM, na.rm = TRUE); range(wcbs_df2$AES3KM)

#' Ranges are reasonably comparable, although WCBS squares do not include the very high local and landscape AES scores. This is expected as these scores are very rare
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

#repeat for WCBS

wcbs_df2$Cat_1KM <- cut(wcbs_df2$AES1KM, 
                       c(-1,500,5000,1e9), 
                       c("Low","Medium","High"))
wcbs_df2$Cat_3KM <- cut(wcbs_df2$AES3KM, 
                       c(-1,500,5000,1e9), 
                       c("Low","Medium","High"))

#table of proportions
Comb <- rbind(select(AES_all, AES1KM = Cat_1KM, AES3KM = Cat_3KM) %>%
                mutate(Survey = "Whole countryside", value = 1),
              select(wcbs_df2, AES1KM = Cat_1KM, AES3KM = Cat_3KM, value= 1) %>%
                mutate(Survey = "WCBS", value = 1)) %>%
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


#' WCBS squares are very similar to all squares in terms of coverage of the 1km x 3km gradient categories. All categories have some representation although only small numbers of squares exist in the Low_High and High_Low categories (only 8 square visits (0.4%) in High_Low vs 405 (24%) in Medium_Medium)
#' 


#' 
#' 3. Evaluate whether AES gradients are confounded with habitat variables
#' 
#' Within the WCBS data we might find that AES scores are highly correlated with other landscape or habitat variables. To assess this we will calculate Spearman rank correlations and plot the correlations
#' 
m1 <- round(cor(wcbs_df2[,c(16:23)], 
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

pairs(wcbs_df2[,c(16:23)], upper.panel = panel.cor,
      col = "black", pch = 20)


#' It looks like neither the 1km or 3km AES scores are strongly correlated with any of the habitat variables. There is evidence of a correlation between 1km and 3km AES scores, and between some of the habitat variables.
#' 
#' 
#' 4. Assess distribution in uplands vs lowlands
#' 
wcbs_df2$SDA_cat <- cut(wcbs_df2$area_m_pc, 
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

sda_dat <- data.frame(SDA_cat = unlist(list(AES_all$SDA_cat, wcbs_df2$SDA_cat)),
                      Dataset = c(rep("All 1km sqs", nrow(AES_all)), rep("WCBS squares", nrow(wcbs_df2))))

ggplot(sda_dat, aes(fill = SDA_cat, x = Dataset))+
  geom_bar(position = "fill", stat = "count")+
  scale_fill_brewer(palette="Paired")
#' 
#' WCBS squares under sample upland areas. Compared to the national coverage of uplands, WCBS only have about half as many upland squares as would be expected if the scheme was completely representative.
#' 
#' 
#' 5. Consider whether CS data are regionally biased
#' 
#' To do this we can map the locations to identify potential areas with low coverage
#+ warning = FALSE
locs <- st_as_sf(data.frame(wcbs_df2[,c(1:5)]), 
                 coords = c("Ecent", "Ncent"))
plot(st_geometry(locs), main = "WCBS square locations", 
     pch = 20, cex = 0.7)
plot(ENG,  add = TRUE)

#' There is evidence of a bias towards the south, with sparser coverage in the Midlands and north. Therefore WCBS may not be entirely representative of the country, although there are still a fair number of squares in the north.
#' 
#' We can also look at whether there are any obvious patterns in the AES scores across the country. A simple way to do this is to look for correlations between the AES gradients and Eastings and Northings. This won't pick up more subtle regional patterns but would pick up if e.g. all the squares with high local AES were in the south.
#' 
pairs(wcbs_df2[,c(3,4,16,17)], upper.panel = panel.cor, labels = c("Easting", "Northing", "AES 1km", "AES 3km"))

#' No correlations between AES gradients and Eastings or Northings
#' 

#' 
#' 6. Summarise differences in protocols
#'
#' REQUIRES INPUT FROM MARC