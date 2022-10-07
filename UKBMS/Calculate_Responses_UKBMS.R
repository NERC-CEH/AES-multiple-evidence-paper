# UKBMS
library(dplyr)
library(readxl)

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

wcbs_ukbms <- read.csv(paste0(fpath, "Citizen science datasets/UKBMS_WCBS/site_data_2017-19.txt"))

# Subset to UKBMS
ukbms_sites <- wcbs_ukbms[wcbs_ukbms$SCHEME == "UKBMS",]

# Get location data
CSlocs <- read_excel(paste0(fpath, "Citizen science datasets/CitizenScienceGridRefs.xlsx"), col_types = c("text", "text", "text", "numeric", "numeric"))
ukbms_locs <- filter(CSlocs, SCHEME == "UKBMS")

# Get species data
wcbs_ukbms_species <- read.csv(paste0(fpath,"Citizen science datasets/UKBMS_WCBS/count _data_2017-2019.txt"))

# filter to UKBMS
ukbms_species <- filter(wcbs_ukbms_species, SCHEME == "UKBMS" & SITENO %in% ukbms_sites$SITENO)

# get transect length and add to site data
# check if each site has always the same transect length
ukbms_species %>% group_by(SITENO) %>% 
  summarise(n = length(unique(TRANSECT_LENGTH))) %>% 
  janitor::tabyl(n)
# every site has maximum 1 transect length
ukbms_transect <- select(ukbms_species, SITENO, TRANSECT_LENGTH) %>%
  unique()
janitor::get_dupes(ukbms_transect, SITENO)
# no duplicates

# all metadata
ukbms_sites <- full_join(ukbms_sites, ukbms_transect) %>%
  select(SITENO, EAST, NORTH, YEAR, N_VISITS_MAYTOAUGUST, TRANSECT_LENGTH) %>%
  full_join(select(ukbms_locs, SITENO = `SITE CODE`,
                   GRIDREF_1KM) %>%
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

UKBMS_RESPONSES <- ukbms_species %>%
  mutate(COMMON_NAME = recode(COMMON_NAME,
                              "Essex Skipper" = "Essex/Small Skipper",
                              "Small Skipper" = "Essex/Small Skipper")) %>%
  group_by(SITENO, YEAR, COMMON_NAME) %>%
  summarise(COUNT = sum(COUNT)) %>%
  summarise(Abundance = sum(COUNT),
            Richness = sum(COUNT>0),
            Shannon_diversity = shandiv(COUNT)) %>%
  inner_join(ukbms_sites) %>%
  # remove site with unfeasibly high abundances
  filter(SITENO != 1063)

summary(UKBMS_RESPONSES)
psych::multi.hist(select_if(UKBMS_RESPONSES, is.numeric))

plot(Abundance ~ N_VISITS_MAYTOAUGUST, UKBMS_RESPONSES)
plot(Abundance ~ N_VISITS_MAYTOAUGUST, UKBMS_RESPONSES)
