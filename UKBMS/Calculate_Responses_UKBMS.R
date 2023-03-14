# UKBMS
library(dplyr)
library(readxl)
library(RODBC)

# useful function
shandiv <- function(x){
  x2 <- x[x>0]
  if(length(x2)<1) return(0)
  x2 <- x2/sum(x2)
  div <- -sum(x2*log(x2))
  return(div)
}


# Get data
dir <- config::get()
fpath <- dir$directories$ukbmsdata

wcbs_ukbms <- read.table(paste0(fpath, "site_date_2017-2021.txt"), sep = "\t", header = TRUE)

# Subset to UKBMS
ukbms_sites <- wcbs_ukbms[wcbs_ukbms$SCHEME == "UKBMS",]

# Get species data
wcbs_ukbms_species <- read.table(paste0(fpath,"count_data_2017-2021_v2.txt"), sep = ",", header = TRUE)

# filter to UKBMS
ukbms_species <- filter(wcbs_ukbms_species, SCHEME == "UKBMS" & SITENO %in% ukbms_sites$SITENO)

# Get transect data
ukbms_transectlength <- read.table(paste0(fpath,"ONLINE_sitedata_2017-2021.txt"), sep = "\t", header = TRUE)

ukbms_species$TRANSECT_LENGTH_NEW <- ukbms_transectlength$LENGTH[match(ukbms_species$SITENO, ukbms_transectlength$SITE_NO)]

ukbms_species$TRANSECT_LENGTH_NEW[is.na(ukbms_species$TRANSECT_LENGTH_NEW)|ukbms_species$TRANSECT_LENGTH_NEW == 0] <- ukbms_species$TRANSECT_LENGTH[is.na(ukbms_species$TRANSECT_LENGTH_NEW)|ukbms_species$TRANSECT_LENGTH_NEW == 0]

transect_mismatch <- unique(ukbms_species[ukbms_species$TRANSECT_LENGTH_NEW != ukbms_species$TRANSECT_LENGTH|is.na(ukbms_species$TRANSECT_LENGTH_NEW)|is.na(ukbms_species$TRANSECT_LENGTH),c(1,2,11)])


#species - site match

unique(ukbms_sites$SITENO[!ukbms_sites$SITENO %in% ukbms_species$SITENO])

#visits
ukbms_visit <- read.table(paste0(fpath,"visit_data_2017-2021.txt"), sep = "\t", header = TRUE)


# Get location data
scpath <- dir$directories$scoredata
CSlocs <- read.csv(paste0(scpath, "Butts_Gradient_Scores.csv"))
CSlocs <- CSlocs[!duplicated(CSlocs[,2:14]),2:14]#remove duplicate rows, some squares have more than one site number due to multiple transects in that location
#remove locations for masked squares
CSlocs <- CSlocs[CSlocs$MaskStatus == "OK",]#also includes squares with NA mask which are coastal
ukbms_locs <- CSlocs[CSlocs$buttsurv.SITENO %in% ukbms_sites$SITENO,]



# get transect length and add to site data
# check if each site has always the same transect length
ukbms_species %>% group_by(SITENO) %>% 
  summarise(n = length(unique(TRANSECT_LENGTH_NEW))) %>% 
  janitor::tabyl(n)
# every site has maximum 1 transect length
ukbms_transect <- dplyr::select(ukbms_species, SITENO, TRANSECT_LENGTH_NEW) %>%
  unique()
janitor::get_dupes(ukbms_transect, SITENO)
# no duplicates

# all metadata
ukbms_sites <- full_join(ukbms_sites, ukbms_transect) %>%
  dplyr::select(SITENO, EAST, NORTH, YEAR, N_VISITS_MAYTOAUGUST, TRANSECT_LENGTH_NEW) %>%
  left_join(dplyr::select(ukbms_locs, SITENO = `buttsurv.SITENO`,
                   buttsurv.GRIDREF_1km) %>%
  mutate(SITENO = as.numeric(SITENO)))


# Deal with taxonomy issues ####
unique(grep("Fritillary",ukbms_species$COMMON_NAME, value = TRUE, ignore.case = TRUE))
# no aggregate here

unique(grep("white",ukbms_species$COMMON_NAME, value = TRUE, ignore.case = TRUE))
# no aggregate here

unique(grep("skipper",ukbms_species$COMMON_NAME, value = TRUE, ignore.case = TRUE))
# aggregate for Essex and Small Skipper

ukbms_species %>% filter(COMMON_NAME == "Painted Lady") %>%
  .$COUNT %>% summary()


db <- config::get("AES")

con <- odbcConnect(db$DSN, db$uid, db$pwd, believeNRows = FALSE)

butttraits <- sqlFetch(con, "TBL_BUTTERFLY_SPECIES")

close(con)

#ugly code here to add traits
ukbms_species$TRAIT <- butttraits$WINGSPAN_CATEGORY[match(ukbms_species$SCI_NAME, butttraits$BUTTERFLY_SPECIES)]

##assign aggregate to trait values
ukbms_species$TRAIT[ukbms_species$SCI_NAME == "Thymelicus lineola/sylvestris"] <- 1


UKBMS_RESPONSES <- ukbms_species %>%
  mutate(COMMON_NAME = recode(COMMON_NAME,
                              "Essex Skipper" = "Essex/Small Skipper",
                              "Small Skipper" = "Essex/Small Skipper")) %>%
  group_by(SITENO, YEAR, COMMON_NAME, TRAIT) %>%
  summarise(COUNT = sum(COUNT),
            LOW_COUNT = sum(COUNT[TRAIT == 1], na.rm = TRUE),
            MED_COUNT = sum(COUNT[TRAIT == 2], na.rm = TRUE),
            HIGH_COUNT = sum(COUNT[TRAIT == 3], na.rm = TRUE)) %>%
  group_by(SITENO, YEAR) %>%
  summarise(Abundance = sum(COUNT),
            Richness = sum(COUNT>0),
            Shannon_diversity = shandiv(COUNT),
            Low_mobility_abund = sum(LOW_COUNT),
            Med_mobility_abund = sum(MED_COUNT),
            High_mobility_abund = sum(HIGH_COUNT)) %>%
  right_join(ukbms_sites) %>%
  # change to right join so sites with no butterflies are still included
  # remove site with unfeasibly high abundances
  filter(SITENO != 1063) %>%
  # remove site with missing species data (4856) 
  filter(SITENO != 4856) %>%
  # remove single species transects
  filter(!(SITENO %in% c(2127,2363,3107,3119,3128,3321,3322,3325)))

#replace NA values in Richness, Abundance and Shannon_diversity with 0
UKBMS_RESPONSES[,c("Abundance", "Richness", "Shannon_diversity",
                   "Low_mobility_abund","Med_mobility_abund","High_mobility_abund")][is.na(UKBMS_RESPONSES[,c("Abundance", "Richness", "Shannon_diversity",
                                                                                                              "Low_mobility_abund","Med_mobility_abund","High_mobility_abund")])] <- 0

#remove transect lengths less than 50m (6 entries)
UKBMS_RESPONSES <- UKBMS_RESPONSES[UKBMS_RESPONSES$TRANSECT_LENGTH_NEW > 49,]
  
#remove masked grid references

UKBMS_RESPONSES <- UKBMS_RESPONSES[!is.na(UKBMS_RESPONSES$buttsurv.GRIDREF_1km),]

UKBMS_RESPONSES <- UKBMS_RESPONSES[!duplicated(UKBMS_RESPONSES),]

summary(UKBMS_RESPONSES)
psych::multi.hist(select_if(UKBMS_RESPONSES, is.numeric))

plot(Abundance ~ N_VISITS_MAYTOAUGUST, UKBMS_RESPONSES)
plot(Richness ~ N_VISITS_MAYTOAUGUST, UKBMS_RESPONSES)
