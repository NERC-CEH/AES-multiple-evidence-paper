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

# Subset to WCBS
wcbs_sites <- wcbs_ukbms[wcbs_ukbms$SCHEME == "WCBS",]

# Get species data
wcbs_ukbms_species <- read.table(paste0(fpath,"count_data_2017-2021_v2.txt"), sep = ",", header = TRUE)

# filter to UKBMS
wcbs_species <- filter(wcbs_ukbms_species, SCHEME == "WCBS" & SITENO %in% wcbs_sites$SITENO)



#species - site match

unique(wcbs_sites$SITENO[!wcbs_sites$SITENO %in% wcbs_species$SITENO])
#one site with no species data


# Get location data
scpath <- dir$directories$scoredata
CSlocs <- read.csv(paste0(scpath, "Butts_Gradient_Scores.csv"))
CSlocs <- CSlocs[!duplicated(CSlocs[,2:14]),2:14]#remove duplicate rows, some squares have more than one site number due to multiple transects in that location
#remove locations for masked squares
CSlocs <- CSlocs[CSlocs$MaskStatus == "OK",]#also includes squares with NA mask which are coastal
wcbs_locs <- CSlocs[CSlocs$buttsurv.SITENO %in% wcbs_sites$SITENO,]



# check if each site has always the same transect length
wcbs_species %>% group_by(SITENO) %>% 
  summarise(n = length(unique(TRANSECT_LENGTH))) %>% 
  janitor::tabyl(n)
# every site has maximum 1 transect length
wcbs_transect <- select(wcbs_species, SITENO, TRANSECT_LENGTH) %>%
  unique()
janitor::get_dupes(wcbs_transect, SITENO)
# no duplicates


# all metadata
wcbs_sites <- full_join(wcbs_sites, wcbs_transect) %>%
  select(SITENO, EAST, NORTH, YEAR, N_VISITS_MAYTOAUGUST, TRANSECT_LENGTH) %>%
  full_join(select(wcbs_locs, SITENO = `buttsurv.SITENO`,
                   buttsurv.GRIDREF_1km) %>%
              mutate(SITENO = as.numeric(SITENO)))
summary(wcbs_sites)

# Deal with taxonomy issues ####

# Deal with taxonomy issues ####
unique(grep("Fritillary",wcbs_species$COMMON_NAME, value = TRUE, ignore.case = TRUE))
# no aggregate here

unique(grep("white",wcbs_species$COMMON_NAME, value = TRUE, ignore.case = TRUE))
# no aggregate here

unique(grep("skipper",wcbs_species$COMMON_NAME, value = TRUE, ignore.case = TRUE))
# aggregate for Essex and Small Skipper

wcbs_species %>% filter(COMMON_NAME == "Painted Lady") %>%
  .$COUNT %>% summary()


db <- config::get("AES")

con <- odbcConnect(db$DSN, db$uid, db$pwd, believeNRows = FALSE)

butttraits <- sqlFetch(con, "TBL_BUTTERFLY_SPECIES")

close(con)

#ugly code here to add traits
wcbs_species$TRAIT <- butttraits$WINGSPAN_CATEGORY[match(wcbs_species$SCI_NAME, butttraits$BUTTERFLY_SPECIES)]


wcbs_responses <- wcbs_species %>%
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
  right_join(wcbs_sites)

#replace one missing transect length
wcbs_responses$TRANSECT_LENGTH <- 2000

#replace NA values in Richness, Abundance and Shannon_diversity with 0
wcbs_responses[,c("Abundance", "Richness", "Shannon_diversity",
                  "Low_mobility_abund","Med_mobility_abund","High_mobility_abund")][is.na(wcbs_responses[,c("Abundance", "Richness", "Shannon_diversity",
                                                                                                            "Low_mobility_abund","Med_mobility_abund","High_mobility_abund")])] <- 0

#remove masked grid references

wcbs_responses <- wcbs_responses[!is.na(wcbs_responses$buttsurv.GRIDREF_1km),]

wcbs_responses <- wcbs_responses[!duplicated(wcbs_responses),]

summary(wcbs_responses)
psych::multi.hist(select_if(wcbs_responses, is.numeric))

plot(Abundance ~ N_VISITS_MAYTOAUGUST, wcbs_responses)
plot(Richness ~ N_VISITS_MAYTOAUGUST, wcbs_responses)
