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
fpath <- config_path

wcbs_ukbms <- read.csv(paste0(fpath, "Citizen science datasets/UKBMS_WCBS/site_data_2017-19.txt"))

# Subset to UKBMS
wcbs_sites <- wcbs_ukbms[wcbs_ukbms$SCHEME == "WCBS",]

# Get location data
CSlocs <- read_excel(paste0(fpath, "Citizen science datasets/CitizenScienceGridRefs.xlsx"), col_types = c("text", "text", "text", "numeric", "numeric"))
wcbs_locs <- filter(CSlocs, SCHEME == "WCBS")

# Get species data
wcbs_ukbms_species <- read.csv(paste0(fpath,"Citizen science datasets/UKBMS_WCBS/count _data_2017-2019.txt"))

# filter to UKBMS
wcbs_species <- filter(wcbs_ukbms_species, SCHEME == "WCBS" & SITENO %in% wcbs_sites$SITENO)

# get transect length and add to site data
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
  full_join(select(wcbs_locs, SITENO = `SITE CODE`,
                   GRIDREF_1KM) %>%
              mutate(SITENO = as.numeric(SITENO)))
summary(wcbs_sites)

# Deal with taxonomy issues ####
wcbs_species %>% filter(COMMON_NAME == "Painted Lady") %>%
  .$COUNT %>% summary()

wcbs_responses <- wcbs_species %>%
  mutate(COMMON_NAME = recode(COMMON_NAME,
                              "Essex Skipper" = "Essex/Small Skipper",
                              "Small Skipper" = "Essex/Small Skipper")) %>%
  group_by(SITENO, YEAR, COMMON_NAME) %>%
  summarise(COUNT = sum(COUNT)) %>%
  summarise(Abundance = sum(COUNT),
            Richness = sum(COUNT>0),
            Shannon_diversity = shandiv(COUNT)) %>%
  inner_join(wcbs_sites)

summary(wcbs_responses)
psych::multi.hist(select_if(wcbs_responses, is.numeric))

plot(Abundance ~ N_VISITS_MAYTOAUGUST, wcbs_responses)
plot(Abundance ~ N_VISITS_MAYTOAUGUST, wcbs_responses)
