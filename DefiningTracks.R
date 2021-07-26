#################################################################
# Script developed by Jen Cruz to clean and format location data #
# We also convert cleaned location data to tracks using atm package   # 
# We rely heavily on amt getting started vignette here:       #
# https://cran.r-project.org/web/packages/amt/vignettes/p1_getting_started.html#
#                                                               #
# Data are prairie falcon locations collected during Spring/Summer #
# of 2021 at Morley Nelson Birds of Prey NCA.                      #
# Data were collected for multiple individuals and at #
# different frequencies including 2 sec intervals when the individuals#
# were moving (every 2-3 days), and 4hr fixes otherwise to define #
# breeding season home range. # Movements also shifted to migration#
# collection after individuals left their breeding grounds. #
#################################################################

################## prep workspace ###############################

# Install new packages from "CRAN" repository if you don't have them. # 
install.packages( "tidyverse" ) #actually a collection of packages 
install.packages( "sp" )
install.packages( "amt" )

# load packages relevant to this script:
library( tidyverse ) #easy data manipulation
library( amt )
library( sp )

#####################################################################
## end of package load ###############

###################################################################
#### Load or create data -----------------------------------------
# Clean your workspace to reset your R environment. #
rm( list = ls() )

# Set working directory. This is the path to your Rstudio folder for this 
# project. If you are in your correct Rstudio project then it should be:
getwd()
# if so then:
workdir <- getwd()

# set path to where you can access your data #
# Note that the path will be different in yours than mine.#
datapath <- "Z:/Common/PrairieFalcons/"

#import GPS data
dataraw <- read.csv( file = paste0( datapath,"ctt_data_export_20210721152221.csv" ),
                     header = TRUE,
                     colClasses = c("serial" = "character"))
#view
head( dataraw ); dim( dataraw )
##############

#######################################################################
######## cleaning data ###############################################
# Data cleaning is crutial to accurate analysis # 

# GPS units often provide information on the quality of the fixes they #
# obtained.#
# These units from Cellular track technologies give us HDOP, VDOP and #
# time to fix information # 
# Start by viewing what those look like in our dataset #

hist( dataraw$vdop, breaks = 50 )
hist( dataraw$hdop, breaks = 50 )
hist( dataraw$time_to_fix )

# Remove 2D fixes and fixes where HDOP or VDOP ≥10 #
# (D’eon and Delparte, 2005; Poessel et al., 2016).#
# Also those where time to fix > 20min or with 0 satellites
unique(dataraw$cog )
#start by creating a new dataframe 
datadf <- dataraw 
#which columns do we have
colnames( datadf )
# Filter to remove inaccurate locations #remove some superfluous columns
datadf <- datadf %>% dplyr::filter( hdop < 10 ) %>%
  dplyr::filter( vdop < 10 ) %>%
  dplyr::filter( time_to_fix <= 20 ) %>% 
  dplyr::filter( nsats > 0 ) %>%
  dplyr::select( -inactivity, -geo, -data_voltage, -solar_current, -solar_charge )

#view
head( datadf ); dim( datadf )
#How many rows did we remove?
# Answer: 
dim( dataraw ) - dim( datadf )
# What % of data did we loose?
# Answer:
# 
# We also need to set a time column containing date and time information #
# in POSIX format
# We rely on lubridate for this. If you haven't use lubridate before #
# go here: https://cran.r-project.org/web/packages/lubridate/vignettes/lubridate.html
# to learn more about how to easily manipulate time and dates in R #
# Data are stored in year,month,day,hour,minute, second format. 
# We define correct format with lubridate and convert it to posixct
datadf$date <- lubridate::ymd_hms( datadf$GPS_YYYY.MM.DD_HH.MM.SS,
              tz = "MST" )
datadf$ts <- as.POSIXct( datadf$date )
#view
head( datadf ); dim( datadf )

# check if any data are missing
all( complete.cases( datadf ) )
# none so we can move on
# Our individual IDs are cumbersome so we create a new id column
datadf$id <- group_indices( datadf, serial )
# we also add month and day of year information using lubridate
datadf <- datadf %>% 
  mutate( mth = lubridate::month(date),
          jday = lubridate::yday(date) )

### Defining coordinate system and projection for the data ######
# location data were recorded using WGS 84  in lat long #
# We use the epsg code to define coordinate system for our sites #
# How? Google epsg WGS 84 # First result should  take you here: #
# https://spatialreference.org/ref/epsg/wgs-84/ 
# where you can find that epgs = 4326 for this coordinate system #
# If you are not familiar with geographic data, vector, rasters & #
# coordinate systems go here: 
# https://bookdown.org/robinlovelace/geocompr/spatial-class.html #
# to learn more. #

# For amt, crs need to be provided using sp package so:
crsdata <- sp::CRS( "+init=epsg:4326" )
# We also want to transform the lat longs to easting and northings #
# using UTM. For this we need to know what zone we are in. Go: #
# http://www.dmap.co.uk/utmworld.htm
crstracks <- sp::CRS( "+proj=utm +zone=12" )
#
# We are now ready to make tracks using atm package. We start with#
# one individual only. But which one? We first check sample size #
table( datadf$id )
# We can also get an idea of the data collection for each individual
# by plotting histograms
ggplot( datadf, aes( x = jday, group = id ) ) +
  theme_classic( base_size = 15 ) +
  geom_histogram( ) +
  facet_wrap( ~ id )
# What do the histograms tell you about the nature of the data #
# Sample size, intensity for different individuals? #
# Answer:
#


# Note that we only have 6 records for the 1st individual. #
# Individual 3 appears to have few records over multiple days #
# while others have high records over few days (probably collecting) #
# data in flight mode #

# We start with indiv 3:
tr.idv <- datadf %>% 
  #select data for the first individual only
  dplyr::filter( id == 3 ) %>% 
  # remove duplicate times
  dplyr::filter( !duplicated( ts ) ) %>% 
  #make track. Note you can add additional columns to it
  amt::make_track(.y = lat, .x = lon, .t = ts, id = id, mth = mth, jday = jday,
              speed = speed, alt = alt, 
              #note that we give it the data CRS to start
              crs = crsdata )
#we check it
class( tr.idv ); head( tr.idv ); dim( tr.idv )
# #we project object to UTM:
tr.idv <- amt::transform_coords( tr.idv, sp::CRS( "+proj=utm +zone=11" ) )
head( tr.idv)
# Note, a common mistake when working with spatial data is forgetting to match their projections.#
# This can introduce significant errors in your analysis. It is particularly common when you #
# are sourcing your data from multiple locations. #

#### obtaining basic data summaries ###################
# What datat summaries are useful regardless of your analysis of #
# interest? 
# Answer:
# 
# We calculate step lengths between locations
# tr.idv <- tr.idv %>%
#   mutate( sl_ = amt::step_lengths(.) )
#view
summary( tr.idv$sl_ )
# We need to ensure that sampling rate is approximately regular:
#view sampling rate:
amt::summarize_sampling_rate( tr.idv )
#from here we see that 2hrs mean may be a good resampling option
tr.slow <- tr.idv %>%  
  amt::track_resample( rate = hours(2),
      tolerance = minutes(20) )
#check
tr.slow
# resample subsets the data to the desired sampling frequency #
# in our case 2 hrs, and adds a 'burst-' column that groups sequence #
# of relocations within equal sampling rates. This accounts for #
# gaps in the data due to missing fixes. #
tr.slow$burst_
# How many groups or bursts do we have for our individual?
# Answer:
#

# We can plot step lengths and turning angles for each burst by:
steps_by_burst( tr.slow )%>% 
  right_join( tr.slow, by = "burst_" ) %>% 
  ggplot( . ) +
  geom_density( aes( y = sl_, fill = as.factor(jday)), alpha = 0.4 ) +
    #xlab("Step length" ) + xlab("") + 
  theme_bw( base_size = 19 )  +
  xlim( 0, 0.004 ) + ylim(0, 7000 )

# plot turning angles:
steps_by_burst( tr.slow ) %>% 
  right_join( tr.slow, by = "burst_" ) %>% 
  gglot( ., aes( x = ta_, fill = as.factor(jday) ) ) +
  geom_bar(stat="identity") +
  coord_polar() +
  ylab("Turning angle") + xlab("") + 
  theme_bw( base_size = 19)  

############# end of script  ###########################################