#################################################################
# Script developed by Jen Cruz to clean and format location data #
# We also convert cleaned location data to tracks using atm package   # 
# We rely heavily on amt getting started vignette here:       #
# https://cran.r-project.org/web/packages/amt/vignettes/p1_getting_started.html#
#                                                               #
# Data are Prairie Falcon locations collected during Spring/Summer #
# of 2021 at Morley Nelson Birds of Prey NCA.                      #
# Data were collected for multiple individuals and at #
# different frequencies including 2 sec intervals when the individuals#
# were moving (every 2-3 days), and 30min? fixes otherwise to define #
# breeding season home range. # Frequency shifted to hourly once #
# individuals left their breeding grounds. #
#################################################################

################## Prep. workspace ###############################

# Install new packages from "CRAN" repository if you don't have them. # 
install.packages( "tidyverse" ) #actually a collection of packages 
install.packages( "sp" )
install.packages( "amt" )
install.packages( "sf" )

# load packages relevant to this script:
library( sp )
library( tidyverse ) #easy data manipulation and plotting
# set option to see all columns and more than 10 rows
options( dplyr.width = Inf, dplyr.print_min = 100 )
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
# Note that the path will be different for your.#
datapath <- "Z:/Common/PrairieFalcons/"

#import GPS data# 
# Fixes are stored as separate CSV files for each individual
## We therefore create a function that imports multiple files at once:
load_data <- function( path ){
  # extract all file names in the folder
  myfiles <- dir( path, pattern = '\\.csv', full.names = TRUE )
  for( i in 1:length(myfiles) ){
    mydata <- read.csv( file = myfiles[i], 
              #remove white spaces  
              strip.white =TRUE, 
              #include column headers
              header = TRUE, 
              # read the serial column as a character instead of number:
              colClasses = c("serial" = "character") ) 
  # create df for first file and append rows for other files
   ifelse( i == 1,
           df <- mydata, 
           df <- bind_rows(df, mydata) )
  } 
  #return df for all individuals
  return( df )
}

#apply function to import all files as list of databases:
dataraw <- load_data( paste0(datapath, 'allindvs/') )
#Note that the files are all in a subdirectory

# Import trapping records with details of when radiotrackers were 
# fitted to the individuals
records <- read.csv( file = paste0( datapath,"survey_0.csv" ),
                     #replaces those values with NA
                     na.strings = c(""," ","NA"), 
                     # include column headings
                     header = TRUE )
#check
head( records ); dim( records )

#import polygon of the NCA as sf spatial file:
NCA_Shape <- sf::st_read("Z:/Common/QCLData/Habitat/NCA/GIS_NCA_IDARNGpgsSampling/BOPNCA_Boundary.shp")
##############

#######################################################################
######## cleaning data ###############################################
# Data cleaning is crucial for accurate analysis # 
# Trapping records provide info on when individuals were fitted #
# with transmitters.#
colnames( records )
# we keep transmitter id, date and sex
records <- records %>% dplyr::select( Telemetry.Unit.ID, Sex, 
                                      Date.and.Time )
#view
records
#convert date to correct format using lubridate
records$StartDate <- lubridate::mdy_hms( records$Date.and.Time, 
                                    tz = "MST")
# Add a day so that we can ignore records from the trapping day #
# and start only with  those from the following day:
records$StartDate <- records$StartDate + lubridate::days(1)
#convert to day of year
records$StartDay <- lubridate::yday( records$StartDate )
#unit IDs were missing starting number of their serial number #
# we append those so we can match it to the GPS serial IDs:
records$serial <- paste0( '894608001201',records$Telemetry.Unit.ID )

#check 
head( records); dim( records)

# Using serial ID unique to each individual found in df, add territory column
# In Territory column each serial ID is linked to its corresponding territory
records <- records %>%
  mutate(territory = case_when(
    endsWith(serial, "47221") ~ "SG",
    endsWith(serial, "47775") ~ "CRW",
    endsWith(serial, "47874") ~ "SDTP",
    endsWith(serial, "48120") ~ "PR_II",
    endsWith(serial, "46751") ~ "HHGS_DS",
    endsWith(serial, "46983") ~ "HHGS_US",
    endsWith(serial, "47197") ~ "Mac",
    endsWith(serial, "48229") ~ "CRW_new",
    endsWith(serial, "48377") ~ "CFR",
  ))

unique(records$territory)

###################################################################
# Clean GPS data
# GPS units often provide information on the quality of the fixes they #
# obtained.#
# The units from Cellular track technologies provide HDOP, VDOP and #
# time to fix information # 
# Start by viewing what those look like in the dataset #

hist( dataraw$vdop, breaks = 50 )
hist( dataraw$hdop, breaks = 50 )
hist( dataraw$time_to_fix )

# Remove 2D fixes and fixes where HDOP or VDOP ≥10 following #
# D’eon and Delparte (2005).#
# Also those where time to fix > 20min or with 0 satellites:

#start by creating a new dataframe to store cleaned location records:
datadf <- dataraw 
#which columns do we have?
colnames( datadf )
# Filter to remove inaccurate locations
datadf <- datadf %>% dplyr::filter( hdop < 10 ) %>%
  dplyr::filter( vdop < 10 ) %>%
  dplyr::filter( time_to_fix <= 20 ) %>% 
  dplyr::filter( nsats > 0 ) %>%
  dplyr::filter( lat > 0 ) %>% 
  #remove superfluous columns
  dplyr::select( -inactivity, -geo, -data_voltage, -solar_current, 
                 -solar_charge )

#view
head( datadf ); dim( datadf )
#How many rows did we remove?
# Answer: 
#
dim( dataraw ) - dim( datadf )
# What % of data did we loose?
# Answer:
# 
# We also need to set a time column containing date and time information #
# in POSIX format (as required by amt)#
# We rely on lubridate for this. If you haven't used lubridate before #
# go here: https://cran.r-project.org/web/packages/lubridate/vignettes/lubridate.html
# to learn more about how to easily manipulate time and dates in R #
# Data are stored in year, month, day, hour, minute, second format in our data. 
# We define correct format with lubridate 
datadf$date <- lubridate::ymd_hms( datadf$GPS_YYYY.MM.DD_HH.MM.SS,
              tz = "MST" )
# and create new column where we convert it to posixct
datadf$ts <- as.POSIXct( datadf$date )
#view
head( datadf ); dim( datadf )

# # check if any data are missing
# all( complete.cases( datadf ) )
# # none so we can move on

# we also add month and day of year information using lubridate
datadf <- datadf %>% 
  mutate( mth = lubridate::month(date),
          jday = lubridate::yday(date) )

# We need to remove records for fixes that were recorded before the #
# units were fitted to the animals so we append relevant information #
# from the records dataframe. We do that by combining datadf to records df#
datadf <- records %>%  dplyr::select( serial, territory, Sex, StartDay ) %>% 
  right_join( datadf, by = "serial" )
#view
head( datadf ); dim( datadf )
# Then using StartDay to filter records, removing those that occurred#
#  earlier when unit was turned on, but not fitted to animal #
datadf <- datadf %>% 
  group_by( serial ) %>% 
  dplyr::filter( jday > StartDay ) %>% ungroup()
#view
head( datadf ); dim( datadf )
# serial IDs are cumbersome so we create a new individual ID column:
datadf$id <- group_indices( datadf, serial )

##################################################################
### Define coordinate system and projection for the data ######
# location data were recorded using WGS84 in lat long #
# We use the epsg code to define coordinate system for our sites #
# How? Google epsg WGS84 # First result should  take you here: #
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
# We choose zone 11:
crstracks <- sp::CRS( "+proj=utm +zone=11" )
#We convert the NCA shapefile to the same projection as our tracks
NCA_Shape <- sf::st_transform( NCA_Shape, crstracks )
# We are now ready to make tracks using atm package
#We first check sample size #
table( datadf$id )
# How many individuals have we dropped so far?
# 
# We can also get an idea of the data collection for each individual
# by plotting histograms
#sampling duration
ggplot( datadf, aes( x = jday, group = id ) ) +
  theme_classic( base_size = 15 ) +
  geom_histogram( ) +
  facet_wrap( ~ id )
#speeds travelled
ggplot( datadf, aes( x = speed, group = id ) ) +
  theme_classic( base_size = 15 ) +
  geom_histogram( ) +
  facet_wrap( ~ id )

# What do the histograms tell you about the nature of the data #
# Sample size, intensity for different individuals? #
# Answer:
#
#Why is the first bar on the speed histograms so tall?
#Answer:
#
# do we need to remove data based on these?
#Answer:
#
#######################################################################
###### Creating tracks, calculating step lengths and turning angles ###
####              for all individuals at the same time:           #####
########################################################################
#amt requires us to turn data into tracks for further analyses.
trks <- datadf %>% 
  #make track. Note you can add additional columns to it
  amt::make_track( .y = lat, .x = lon, .t = ts, 
    #define columns that you want to keep, relabel if you need:
    id = id, territory = territory,
    sex = Sex, mth = mth,jday = jday, speed = speed, alt = alt, 
    #assign correct crs
    crs = crsdata )

# Reproject to UTM to convert lat lon to easting northing:
trks <- amt::transform_coords( trks, crstracks )
#Turn into a tibble list by groupping and nest by individual IDs:
trks <- trks %>%  amt::nest( data = -"id" )
#view
trks

# Remember we have multiple types of data including detailed data for flights #
# 3 times a week, 30min fixes during the day, then hourly fixes during #
# migration. We start by focusing on data during breeding season. #
# That means we need to remove migration locations.
# How do we know when individuals started migrating North?
# We plot overall paths for each individual:
for( i in 1:dim(trks)[1]){
  a <- as_sf_points( trks$data[[i]] ) %>% 
    ggplot(.) + theme_bw(base_size = 17) +
    labs( title = paste0('individual =', trks$territory[i]) ) +
    geom_sf(data = NCA_Shape, inherit.aes = FALSE ) +
    geom_sf() 
  print(a)
} 
# Which ones have migration paths?
# Answer:
#
# Any ideas on how to remove migration data?
# Answer:
# 
# Here we rely on NCA polygon, removing records that exist East of the #
# NCA. We can extra the extent of a polygon:
sf::st_bbox( NCA_Shape )
#Then use the Eastern-most coordinate to filter out data 
xmax <- as.numeric(st_bbox(NCA_Shape)$xmax) #627081.5
#subset those tracks less than as breeding and those > as migrating:
trks <- trks %>% mutate(
  breeding = map( data, ~ filter(., x_ < xmax ) ),
  migrating = map( data, ~ filter(., x_ >= xmax ) ) )

#view
trks
# Note we created two other groups of tibbles for the breeding season
# and migrating season #
# Plot step lengths
for( i in 1:dim(trks)[1]){
  a <-  steps( trks$breeding[[i]] ) %>% 
  #a <-  steps( trks$migrating[[i]] ) %>% 
    mutate( jday = lubridate::yday( t1_ ) ) %>% 
    group_by( jday ) %>% 
    summarise( sl_ = log( sum(sl_) ) ) %>% 
    ggplot(.) + theme_bw(base_size = 17) +
    labs( title = paste0('individual =', trks$territory[i]) ) +
    geom_line( aes( y = sl_, x = jday))
  print(a)
}

# We focus on breeding season data:
# Estimate sampling rate for each individual by looping through 
# data using purr function map
sumtrks <- trks %>%  summarize( 
  map( breeding, amt::summarize_sampling_rate ) )
#view
sumtrks[[1]]

# Add tibbles with added step lengths calculated by bursts from #
# breeding season data:
trks.all <- trks %>% mutate(
  steps = map( breeding, function(x) 
    x %>%  track_resample( rate = seconds(5), 
                           tolerance = seconds(5)) %>% 
      steps_by_burst() ) )
#view
trks.all


# plot autocorrelation for step lengths for all individuals
par( mfrow = c( 2,3 ) )
for( i in 1:dim(trks.all)[1] ){
  #extract individual ids
  idd <- trks.all$id[i]
  #use tibbles we calculated in steps
  x <- pull( trks.all[["steps"]][[i]], direction_p )
  #remove missing data
  x <- x[!is.na(x)]
  #calculate autocorrelation function:
  acf( x, lag.max = 300,
       main = paste0( "individual = ",idd ) )
  #Note you can modify the lag.max according to your data 
}
# What would be a reasonable rate to resample at?
# Answer:
# 
# I choose 30min
trks.all <- trks.all %>% 
  mutate(red = map(breeding, function(x ) x %>%  
               track_resample( rate = minutes(30),
               tolerance = minutes(5) ) ) )
#view
trks.all

# We can now unnest the dataframes of interest
#Starting with all breeding season data
trks.breed <- trks.all %>% select( id, breeding ) %>% 
  unnest( cols = breeding ) 
head( trks.breed )

# Now breeding season data, without autocorrelation:
trks.red <- trks.all %>% select( id, red ) %>% 
  unnest( cols = red ) 
head( trks.red )

# Last all migration data:
trks.mig <- trks.all %>% select( id, migrating ) %>% 
  unnest( cols = migrating ) 
head( trks.mig )
###############

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
# burst_ id grouping fixes into separate paths #
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
# Temporal autocorrelation of locations leads to underestimation in #
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

#Reduce individual tracks based on selected steps above 
tr.red <- tr.slow %>% select( x_= x1_, y_=y1_, t_=t1_ ) %>% 
  left_join( tr.idv, by = c('x_', "y_", "t_" ) )
# view
tr.red
# Note that the output is a tibble. Turn it back to a track:
tr.red <- tr.red %>%  
  amt::make_track(.y = y_, .x = x_, .t = t_, id = id, 
            mth = mth, jday = jday, speed = speed, alt = alt, 
            crs = crstracks )

#############################################################################
# Saving relevant objects and data ---------------------------------
#save hourly detection dataframe with weather predictors
# write.csv(x = det_df, 
#           #ensure that you save it onto your datafolder
#           file = paste0( datapath, 'stoc_det_df.csv'), 
#           row.names = FALSE )
#save workspace in case we need to make changes
save.image( "TracksWorkspace.RData" )

########## end of save #########################
############### END OF SCRIPT ########################################