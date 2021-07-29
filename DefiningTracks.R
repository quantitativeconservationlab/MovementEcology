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
install.packages( "sf" )
# load packages relevant to this script:
library( sp )
library( tidyverse ) #easy data manipulation
library( amt ) #creating tracks from location data
library( sf ) #handling spatial data
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

#import GPS data# 
# records are stored as separate CSV files for each individual
## We therefore create a function that imports multiple files:
load_data <- function( path ){
  myfiles <- dir( path, pattern = '\\.csv', full.names = TRUE )
  for( i in 1:length(myfiles) ){
    mydata <- read.csv( file = myfiles[i], strip.white =TRUE, #removes white spaces  
              header = TRUE, 
              colClasses = c("serial" = "character") ) 
   ifelse( i == 1,
           df <- mydata, 
           df <- bind_rows(df, mydata) )
  } 
  return( df )
}

#apply function to import all files as list of databases:
dataraw <- load_data( paste0(datapath, 'allindvs/') )

#import trapping records with details of when radiotrackers were 
# fitted to the individuals
records <- read.csv( file = paste0( datapath,"survey_0.csv" ),
                     #replaces those values with NA, includes column heading
                     na.strings = c(""," ","NA"), header = TRUE )
head( records ); dim( records )

#import polygon of the NCA
NCA_Shape <- st_read("Z:/Common/QCLData/Habitat/NCA/GIS_NCA_IDARNGpgsSampling/BOPNCA_Boundary.shp")
##############

#######################################################################
######## cleaning data ###############################################
# Data cleaning is crucial to accurate analysis # 
# Trapping records provide info on when individuals were fitted #
# with transmitters.#
colnames( records )
# we keep transmitter id, date and sex
records <- records %>% dplyr::select( Telemetry.Unit.ID, Sex, StartDay )
#view
records
#convert date to correct format using lubridate
records$StartDate <- lubridate::mdy_hms( records$Date.and.Time, 
                                    tz = "MST")
#add a day so that we can ignore records from the trapping day #
# and start only with  those from the following day:
records$StartDate <- records$StartDate + lubridate::days(1)
#convert to day of year
records$StartDay <- lubridate::yday( records$StartDate )
#unit IDs were missing starting number of their serial number #
# we append those so we can match it to our GPS data
records$serial <- paste0( '894608001201',records$Telemetry.Unit.ID )

###################################################################
# Now we clean GPS data
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

#start by creating a new dataframe 
datadf <- dataraw 
#which columns do we have
colnames( datadf )
# Filter to remove inaccurate locations 
#remove some superfluous columns
datadf <- datadf %>% dplyr::filter( hdop < 10 ) %>%
  dplyr::filter( vdop < 10 ) %>%
  dplyr::filter( time_to_fix <= 20 ) %>% 
  dplyr::filter( nsats > 0 ) %>%
  dplyr::filter( lat > 0 ) %>% 
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

# we also add month and day of year information using lubridate
datadf <- datadf %>% 
  mutate( mth = lubridate::month(date),
          jday = lubridate::yday(date) )

# We need to remove records for fixes that were recorded before the #
# units were fitted to the animals so we append relevant information #
# from the records dataframe
datadf <- records %>%  dplyr::select( serial, Sex, StartDay ) %>% 
  right_join( datadf, by = "serial" )
#view
head( datadf ); dim( datadf )
#filter records to remove those early ones
datadf <- datadf %>% 
  group_by( serial ) %>% 
  filter( jday > StartDay ) %>% ungroup()
#view
head( datadf ); dim( datadf )
# Our individual IDs are cumbersome so we create a new id column
datadf$id <- group_indices( datadf, serial )

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
crstracks <- sp::CRS( "+proj=utm +zone=11" )
#We convert the NCA shapefile to the same projection as our tracks
NCA_Shape <- sf::st_transform( NCA_Shape, crstracks )
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

###########################################################################
###### Creating tracks, calculating step lengths and turning angles ######
#                     for a single individual                             #
######################
# We start by creating a track for a single individual:
tr.idv <- datadf %>% 
  #select data for one individual only
  dplyr::filter( id == 1 ) %>% 
  # remove duplicate times
  dplyr::filter( !duplicated( ts ) ) %>% 
  #make track. Note you can add additional columns to it
  amt::make_track(.y = lat, .x = lon, .t = ts, id = id, mth = mth, 
        jday = jday, sex = Sex, speed = speed, alt = alt, 
              #note that we give it the data CRS to start
              crs = crsdata )
#check it
class( tr.idv ); head( tr.idv ); dim( tr.idv )
# A common mistake when working with spatial data is forgetting to #
# set data to correct projection, which can introduce significant errors #
# in your analysis #
# Project object to UTM:
tr.idv <- amt::transform_coords( tr.idv, crstracks )
tr.idv
# Why do we change it to UTM?
# Answer: 
#
# We check our sampling rate:
tr.idv %>% amt::summarize_sampling_rate()
# What is it telling us?
# Answer:
#

# Resample track to high resolution frequency so that we can add a #
# burst_ id grouping fixes into separte paths #
# This accounts for gaps in the data due to missing fixes or uneven sampling
# rates #
tr.3 <- tr.idv %>%  
  amt::track_resample( rate = seconds(5),
                       tolerance = seconds(5) )
tr.3
# Convert resulting track to steps, while taking into account the grouping #
# set by bursts_:
tr.3 <- amt::steps_by_burst( tr.3 )
# How many groups or bursts do we have for our individual?
# Answer:
#
# Some analyses require independence of your fix locations. #
# Temporal autocorrelation of locations leads to underestimation of #
# home range size and bias in predictions of habitat selection, core area, #
# and intensity of resource use for those methods that rely on it. #
# We therefore need to create an additional thinned dataset that removes #
# autocorrelation in step lengths and turning angles for our individual track#

# Start by check autocorrelation of track based on direction and turning angles:
acf( tr.3[,"direction_p"] )
acf( tr.3[!is.na(tr.3[,'ta_']),'ta_'] )
# What do these plot tell us? 
# Answer:
# 
# Adjust sampling rate based on results from acf plots above#
tr.slow <- tr.idv %>%  
  amt::track_resample( rate = minutes(1),
                       tolerance = seconds(5) )
# Recalculate metrics and recheck autocorrelation
tr.slow <- steps_by_burst( tr.slow )
acf( tr.slow[,"direction_p"] )
acf( tr.slow[!is.na(tr.slow[,'ta_']),'ta_'] )
# What can you see in the new plots?
# Answer:
# 
#How much data did we loose with this resampling strategy?
tr.slow
dim(tr.idv)[1] - dim(tr.slow)[1] 
#Answer:
#
# Comment on what this means regarding sample size, etc
# Answer:
#
# We can plot step lengths and turning angles for each burst by:
tr.slow %>% 
  ggplot(.) +
  geom_density( aes( x = sl_, fill = as.factor(burst_)), alpha = 0.4 ) +
  xlab("Step length" ) + 
  #ylim( 0, 0.01 ) + xlim(0, 2000 ) +
  theme_bw( base_size = 19 )  +
  theme( legend.position = "none" )
#What does the plot tell us about the step lengths traveled by the individual?
# Answer:

#
# Turning angles:
tr.slow %>% 
  ggplot( ., aes( x = ta_, y = burst_ ) ) +
  geom_bar(stat="identity") +
  coord_polar() +
  ylab("Turning angle") + xlab("") + 
  theme_bw( base_size = 19)  
# Is there any evidence of biased movements for this individual?
# Answer:
# 

#We now reduce our individual track based to match our reduced steps above 
tr.red <- tr.slow %>% select( x_= x1_, y_=y1_, t_=t1_ ) %>% 
  left_join( tr.idv, by = c('x_', "y_", "t_" ) )
# view
tr.red
# Note that the output is now a tibble. We turn it back to a track
tr.red <- tr.red %>%  
  amt::make_track(.y = y_, .x = x_, .t = t_, id = id, 
            mth = mth, jday = jday, speed = speed, alt = alt, 
            crs = crstracks )

##################################################################
#### Working with all individuals at the same time:           #####
########################################################################
trks <- datadf %>% 
  #make track. Note you can add additional columns to it
  amt::make_track(.y = lat, .x = lon, .t = ts, 
        #define columns that you want to keep
        id = id, sex = Sex, mth = mth,jday = jday, speed = speed, alt = alt, 
        #assign crs
        crs = crsdata )

# Reproject to UTM:
trks <- amt::transform_coords( trks, crstracks )
#Group and nest by individual IDs:
trks <- trks %>%  amt::nest( data = -"id" )
#view
trks
for( i in 1:dim(trks)[1]){
a <-  steps( trks$data[[i]] ) %>% 
    mutate( jday = lubridate::yday( t1_ ) ) %>% 
    group_by( jday ) %>% 
    summarise( sl_ = log( sum(sl_) ) ) %>% 
    ggplot(.) + theme_bw(base_size = 17) +
    labs( title = paste0('individual =', trks$id[i]) ) +
    geom_line( aes( y = sl_, x = jday))
print(a)
}

# How do we know which individuals started migrating North?
for( i in 1:dim(trks)[1]){
  a <- as_sf_points( trks$data[[i]] ) %>% 
    ggplot(.) + theme_bw(base_size = 17) +
    labs( title = paste0('individual =', trks$id[i]) ) +
    geom_sf(data = NCA_Shape, inherit.aes = FALSE ) +
    geom_sf() 
  print(a)
} 
# Which ones have migration paths?
# Answer:
#
# We have multiple types of data including detailed data for flights and #
# foraging paths, constant data when they are stationary, then a shift #
# to hourly tracking after they migrate. #

# Sampling rate for each individual by looping through 
# the data using purrr function map
sumtrks <- trks %>%  summarize( 
  map( data, amt::summarize_sampling_rate ) )
#view
sumtrks[[1]]

# Add steps by bursts
trks.all <- trks %>% mutate(
  steps = map( data, function(x) 
    x %>%  track_resample( rate = seconds(5), 
    tolerance = seconds(5)) %>% 
    steps_by_burst() ) )
#view
trks.all


# plot autocorrelation for step lengths for all individuals
par( mfrow = c( 3,3 ) )
for( i in 1:dim(trks.all)[1] ){
  #you can modify the lag.max according to your data 
  idd <- trks.all$id[i]
  x <- pull( trks.all[["steps"]][[i]], direction_p )
  x <- x[!is.na(x)]
  acf( x, lag.max = 300,
       main = paste0("individual = ",idd ) )
}
# What would be a reasonable rate to resample at?
# Answer:
# 
# I choose 25min
trks.all <- trks.all %>% 
  mutate(red = map(data, function(x ) x %>%  
                     track_resample( rate = minutes(25),
                                     tolerance = minutes(5) ) ) )
#view
trks.all



# We can now unnest the dataframes of interest
#Starting with the raw data:
trks.comp <- trks.all %>% select( id, data ) %>% 
  unnest( cols = data ) 
head( trks.comp )

# Now the reduced tracks:
trks.red <- trks.all %>% select( id, red ) %>% 
  unnest( cols = red ) 
head( trks.red )

#############################################################################
# Saving relevant objects and data ---------------------------------
#save hourly detection dataframe with weather predictors
# write.csv(x = det_df, 
#           #ensure that you save it onto your datafolder
#           file = paste0( datapath, 'stoc_det_df.csv'), 
#           row.names = FALSE )
#save workspace in case we need to make changes
save.image( "TracksWorkspace.RData" )

############### END OF SCRIPT ########################################

############# end of script  ###########################################