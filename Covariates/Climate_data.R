# This script is for reading the rainfall and temperature data from the HadUK-Grid
# files provided by the MetOffice, extracting the relevant summaries for the sites of
# interest and calculating the 5-year average for each site from 2017-2021.

# The six climate variables we are extracting are:
# * Annual rainfall
# * Winter rainfall
# * Summer rainfall
# * Annual average temperature
# * Winter minimum temperature
# * Summer maximum temperature

# All files were downloaded from the Met Office on the 2nd November. 

# Citation: Met Office; Hollis, D.; McCarthy, M.; Kendon, M.; Legg, T. (2022):
# HadUK-Grid Gridded Climate Observations on a 1km grid over the UK, v1.1.0.0
# (1836-2021). NERC EDS Centre for Environmental Data Analysis, 26 May 2022.
# doi:10.5285/bbca3267dc7d4219af484976734c9527.
# http://dx.doi.org/10.5285/bbca3267dc7d4219af484976734c9527

# The five year average comprises the five years before survey (NOT inclusive of
# survey). Note that the seasonal files contains 4 columns, in order: winter, spring,
# summer, autumn. Winter is the winter that ends in that year - i.e. winter in 2017
# is the winter 16/17. So that will have to be accounted for - we want the winter
# variable to start from  the previous winter, not the winter two past

# Load libraries
library(raster, exclude = "select")
library(odbc)
library(DBI)
library(sf)
library(dplyr)

# Get locations of the rainfall files
dir <- config::get()
cpath <- dir$directories$climatedata
files <- list.files(cpath)
scpath <- dir$directories$scoredata

# Get the location data for all sites
# UKBMS and WCBS
fpath <- dir$directories$ukbmsdata
wcbs_ukbms <-  read.table(paste0(fpath, "site_date_2017-2021.txt"), sep = "\t", 
                          header = TRUE)

# LandSpAES
db <- config::get("AES")
con <- dbConnect(odbc(), db$DSN, uid = db$uid, pwd = db$pwd)
landspaes <- dbReadTable(con, "TBL_SURVEY_SQUARE") %>%
  mutate(SURVEY_SQUARE = sapply(strsplit(SURVEY_SQUARE, "/"),"[",1)) 
landspaes_loc <- mutate(landspaes,
                        EAST = rnrfa::osg_parse(landspaes$SURVEY_SQUARE)$easting,
                        NORTH = rnrfa::osg_parse(landspaes$SURVEY_SQUARE)$northing,
                        SCHEME = "LandSpAES") %>%
  select(SCHEME, GRIDREF = SURVEY_SQUARE, EAST, NORTH)
landspaes_year <- data.frame(GRIDREF = rep(landspaes_loc$GRIDREF, each = 5),
                             YEAR = rep(2017:2021,nrow(landspaes_loc)))

# Combine sites into one data frame
full_data_year <- landspaes_loc %>%
  full_join(landspaes_year, by="GRIDREF") %>%
  full_join(select(wcbs_ukbms,  GRIDREF, SCHEME, EAST, NORTH, YEAR))

# Turn into spatial object using sf
all_sites <- select(full_data_year, SCHEME, GRIDREF, EAST, NORTH) %>%
  distinct() %>%
  st_as_sf(coords = c("EAST","NORTH"), crs = 27700)

# extract climate data
# rainfall
# read in rainfall data for every site for 2012-2020 (or 2013-2021 for winter)
rain <- lapply(2012:2020, function(x){
  annual_file <- grep(paste0("rainfall_.*_ann_",x),files, value = TRUE)
  annual_rain <- brick(paste0(cpath, annual_file))
  annual_siterain <- extract(annual_rain, all_sites)
  rownames(annual_siterain) <- all_sites$GRIDREF
  colnames(annual_siterain) <- paste0("Annual_",x)
  
  if(anyNA(annual_siterain)){
    annual_rain_NAsitess <- all_sites[is.na(annual_siterain[,1]),]
    annual_siterainNA <- as.matrix(extract(annual_rain, annual_rain_NAsitess, fun = mean,
                                           buffer = 2000), ncol=  1)
    rownames(annual_siterainNA) <- annual_rain_NAsitess$GRIDREF
    colnames(annual_siterainNA) <- paste0("Annual_",x)
  } else{
    annual_siterainNA <- matrix(dimnames = list(c("NA",paste0("Annual_",x))))
  }
  
  annual_siterain <- rbind(
    as.data.frame(annual_siterain) %>%
      tibble::rownames_to_column("GRIDREF"),
    as.data.frame(annual_siterainNA) %>%
      tibble::rownames_to_column("GRIDREF")
  ) %>% na.omit()
  
  # Summer
  summer_file <- grep(paste0("rainfall_.*_seas_",x),files, value = TRUE)
  summer_rain <- brick(paste0(cpath, summer_file))
  summer_siterain <- extract(summer_rain, all_sites)
  rownames(summer_siterain) <- all_sites$GRIDREF
  summer_siterain <- summer_siterain[,3,drop=FALSE]
  colnames(summer_siterain) <- paste0("Summer_",x)
  
  if(anyNA(summer_siterain)){
    summer_rain_NAsitess <- all_sites[is.na(summer_siterain[,1]),]
    summer_siterainNA <- extract(summer_rain, annual_rain_NAsitess, fun = mean,
                                 buffer = 2000)
    rownames(summer_siterainNA) <- summer_rain_NAsitess$GRIDREF
    summer_siterainNA <- summer_siterainNA[,3,drop=FALSE]
    colnames(summer_siterainNA) <- paste0("Summer_",x)
  } else{
    summer_siterainNA <- matrix(dimnames = list(c("NA",paste0("Summer_",x))))
  }
  
  summer_siterain <- rbind(
    as.data.frame(summer_siterain) %>%
      tibble::rownames_to_column("GRIDREF"),
    as.data.frame(summer_siterainNA) %>%
      tibble::rownames_to_column("GRIDREF")
  ) %>% na.omit()
  
  # Winter
  winter_file <- grep(paste0("rainfall_.*_seas_",x+1),files, value = TRUE)
  winter_rain <- brick(paste0(cpath, winter_file))
  winter_siterain <- extract(winter_rain, all_sites)
  rownames(winter_siterain) <- all_sites$GRIDREF
  winter_siterain <- winter_siterain[,1,drop=FALSE]
  colnames(winter_siterain) <- paste0("Winter_",x)
  
  if(anyNA(winter_siterain)){
    winter_rain_NAsitess <- all_sites[is.na(winter_siterain[,1]),]
    winter_siterainNA <- extract(winter_rain, annual_rain_NAsitess, fun = mean,
                                 buffer = 2000)
    rownames(winter_siterainNA) <- winter_rain_NAsitess$GRIDREF
    winter_siterainNA <- winter_siterainNA[,1,drop=FALSE]
    colnames(winter_siterainNA) <- paste0("Winter_",x)
  } else{
    winter_siterainNA <- matrix(dimnames = list(c("NA",paste0("Winter_",x))))
  }
  
  winter_siterain <- rbind(
    as.data.frame(winter_siterain) %>%
      tibble::rownames_to_column("GRIDREF"),
    as.data.frame(winter_siterainNA) %>%
      tibble::rownames_to_column("GRIDREF")
  ) %>% na.omit()
  
  all_data <- full_join(annual_siterain, summer_siterain, by = "GRIDREF") %>%
    full_join(winter_siterain, by = "GRIDREF")
  
})
# Combine list of dataframes into one dataframe
rain_dat <- lapply(rain, tibble::column_to_rownames,"GRIDREF") %>%
  bind_cols() %>%
  tibble::rownames_to_column("GRIDREF") %>%
  tidyr::pivot_longer(-GRIDREF, names_sep = "_", names_to = c("Season","Year"),
                      values_to = "rainfall")

# Calculate 5 year average for every year of survey
rain_summary <- lapply(2017:2021, function(x){
  year_range <- (x-5):(x-1)
  rain_summ <- rain_dat %>%
    filter(Year %in% year_range) %>%
    group_by(GRIDREF, Season) %>%
    summarise(rainfall = mean(rainfall), .groups = "drop") 
  colnames(rain_summ) <- c("GRIDREF","Season",x)
  rain_summ
})
# Combine list into one dataframe
rain_summary <- full_join(rain_summary[[1]],rain_summary[[2]]) %>%
  full_join(rain_summary[[3]]) %>%
  full_join(rain_summary[[4]]) %>%
  full_join(rain_summary[[5]]) %>%
  tidyr::pivot_longer(contains("20"), names_to ="YEAR", values_to = "rainfall") %>%
  mutate(Season = paste(Season,"Rainfall",sep ="_")) %>%
  tidyr::pivot_wider(names_from = Season, values_from = "rainfall")


# temperature
# read in temperature data for every site for 2012-2020 (or 2013-2021 for winter).
# Note the reading in of the daily mean (tas) for annual, maximum for summer (tasmax)
# and minimum (tasmin) for winter
temp <- lapply(2012:2020, function(x){
  annual_file <- grep(paste0("tas_.*_ann_",x),files, value = TRUE)
  annual_temp <- brick(paste0(cpath, annual_file))
  annual_sitetemp <- extract(annual_temp, all_sites)
  rownames(annual_sitetemp) <- all_sites$GRIDREF
  colnames(annual_sitetemp) <- paste0("Annual_",x)
  
  if(anyNA(annual_sitetemp)){
    annual_temp_NAsitess <- all_sites[is.na(annual_sitetemp[,1]),]
    annual_sitetempNA <- as.matrix(extract(annual_temp, annual_temp_NAsitess, fun = mean,
                                           buffer = 2500), ncol=  1)
    rownames(annual_sitetempNA) <- annual_temp_NAsitess$GRIDREF
    colnames(annual_sitetempNA) <- paste0("Annual_",x)
  } else{
    annual_sitetempNA <- matrix(dimnames = list(c("NA",paste0("Annual_",x))))
  }
  
  annual_sitetemp <- rbind(
    as.data.frame(annual_sitetemp) %>%
      tibble::rownames_to_column("GRIDREF"),
    as.data.frame(annual_sitetempNA) %>%
      tibble::rownames_to_column("GRIDREF")
  ) %>% na.omit()
  
  # Summer
  summer_file <- grep(paste0("tasmax_.*_seas_",x),files, value = TRUE)
  summer_temp <- brick(paste0(cpath, summer_file))
  summer_sitetemp <- extract(summer_temp, all_sites)
  rownames(summer_sitetemp) <- all_sites$GRIDREF
  summer_sitetemp <- summer_sitetemp[,3,drop=FALSE]
  colnames(summer_sitetemp) <- paste0("Summer_",x)
  
  if(anyNA(summer_sitetemp)){
    summer_temp_NAsitess <- all_sites[is.na(summer_sitetemp[,1]),]
    summer_sitetempNA <- extract(summer_temp, annual_temp_NAsitess, fun = mean,
                                 buffer = 2000)
    rownames(summer_sitetempNA) <- summer_temp_NAsitess$GRIDREF
    summer_sitetempNA <- summer_sitetempNA[,3,drop=FALSE]
    colnames(summer_sitetempNA) <- paste0("Summer_",x)
  } else{
    summer_sitetempNA <- matrix(dimnames = list(c("NA",paste0("Summer_",x))))
  }
  
  summer_sitetemp <- rbind(
    as.data.frame(summer_sitetemp) %>%
      tibble::rownames_to_column("GRIDREF"),
    as.data.frame(summer_sitetempNA) %>%
      tibble::rownames_to_column("GRIDREF")
  ) %>% na.omit()
  
  # Winter
  winter_file <- grep(paste0("tasmin_.*_seas_",x+1),files, value = TRUE)
  winter_temp <- brick(paste0(cpath, winter_file))
  winter_sitetemp <- extract(winter_temp, all_sites)
  rownames(winter_sitetemp) <- all_sites$GRIDREF
  winter_sitetemp <- winter_sitetemp[,1,drop=FALSE]
  colnames(winter_sitetemp) <- paste0("Winter_",x)
  
  if(anyNA(winter_sitetemp)){
    winter_temp_NAsitess <- all_sites[is.na(winter_sitetemp[,1]),]
    winter_sitetempNA <- extract(winter_temp, annual_temp_NAsitess, fun = mean,
                                 buffer = 2000)
    rownames(winter_sitetempNA) <- winter_temp_NAsitess$GRIDREF
    winter_sitetempNA <- winter_sitetempNA[,1,drop=FALSE]
    colnames(winter_sitetempNA) <- paste0("Winter_",x)
  } else{
    winter_sitetempNA <- matrix(dimnames = list(c("NA",paste0("Winter_",x))))
  }
  
  winter_sitetemp <- rbind(
    as.data.frame(winter_sitetemp) %>%
      tibble::rownames_to_column("GRIDREF"),
    as.data.frame(winter_sitetempNA) %>%
      tibble::rownames_to_column("GRIDREF")
  ) %>% na.omit()
  
  all_data <- full_join(annual_sitetemp, summer_sitetemp, by = "GRIDREF") %>%
    full_join(winter_sitetemp, by = "GRIDREF")
  
})
# Combine list into one dataframe
temp_dat <- lapply(temp, tibble::column_to_rownames,"GRIDREF") %>%
  bind_cols() %>%
  tibble::rownames_to_column("GRIDREF") %>%
  tidyr::pivot_longer(-GRIDREF, names_sep = "_", names_to = c("Season","Year"),
                      values_to = "temp")

# Calculate 5-year average for every year of survey
temp_summary <- lapply(2017:2021, function(x){
  year_range <- (x-5):(x-1)
  temp_summ <- temp_dat %>%
    filter(Year %in% year_range) %>%
    group_by(GRIDREF, Season) %>%
    summarise(temp = mean(temp), .groups = "drop") 
  colnames(temp_summ) <- c("GRIDREF","Season",x)
  temp_summ
})
# Combine list of files into one dataframe
temp_summary <- full_join(temp_summary[[1]],temp_summary[[2]]) %>%
  full_join(temp_summary[[3]]) %>%
  full_join(temp_summary[[4]]) %>%
  full_join(temp_summary[[5]]) %>%
  tidyr::pivot_longer(contains("20"), names_to ="YEAR", values_to = "temp") %>%
  mutate(Season = paste(Season,"Temp",sep ="_")) %>%
  tidyr::pivot_wider(names_from = Season, values_from = "temp")

# Combine rainfall and temperature data into one climate summary file
climdata <- full_join(rain_summary,temp_summary)

write.csv(climdata,paste0(scpath,"HADUK_Rain_Temp_Summaries.csv"),
          row.names = FALSE)

# Checks:
# Are there any NAs?
anyNA(climdata)
# No

# Are the ranges of the variables reasonable?
summary(climdata)
# Yes

# are there five years for every grid reference?
climdata %>% ungroup%>% group_by(GRIDREF) %>% count() %>% ungroup %>% count(n)
# Yes

# Are there squares missing?
wcbs_ukbms$GRIDREF[!(wcbs_ukbms$GRIDREF %in% climdata$GRIDREF)]
# Yes - compare to AES scores and see if those squares are also missing there
AES <- read.csv(paste0(scpath, "Butts_Gradient_Scores.csv"))
filter(AES, GRIDREF %in% c("TA406116", "TA399109"))
# No AES scores for these squares, and looking at them on the map indicates they are
# outside of England as on thin peninsula - assume acceptable for these to be missing
# from climate data
