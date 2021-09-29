



# load packages relevant to this script:
library( sp )
library( tidyverse ) #easy data manipulation and plotting
# set option to see all columns and more than 10 rows
options( dplyr.width = Inf, dplyr.print_min = 100 )
library( amt ) #creating tracks from location data
library( sf ) #handling spatial data

# Load or create data ----------------------------------------------------------
# Clean your workspace to reset your R environment
rm( list = ls() )

# Set working directory. This is the path to your Rstudio folder for this 
# project. If you are in your correct Rstudio project then it should be:
getwd()
# if so then:
workdir <- getwd()

# set path to where you can access your data #
# Note that the path will be different for your.#
datapath <- "Z:/Common/PrairieFalcons/"

# Import GPS data 
# Fixes are stored as separate CSV files for each individual
# We therefore create a function that imports multiple files at once:
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
    # create df for first file and attach rows for other files
    ifelse( i == 1,
            df <- mydata, 
            df <- bind_rows(df, mydata) ) ## I'm not sure what's going on here
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
# Data cleaning is crucial to accurate analysis # 
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
# D’eon and Delparte (2005); Poessel et al. (2016).#
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
head( datadf ); dim( datadf ) ## Whats's the difference between date and ts? They look the same..

# check if any data are missing
all( complete.cases( datadf ) )
# none so we can move on

# we also add month and day of year information using lubridate
datadf <- datadf %>% 
  mutate( mth = lubridate::month(date),
          jday = lubridate::yday(date) )

# We need to remove records for fixes that were recorded before the #
# units were fitted to the animals so we append relevant information #
# from the records dataframe. We do that by combining datadf to records df#
datadf <- records %>%  dplyr::select( serial, Sex, StartDay ) %>% 
  right_join( datadf, by = "serial" )
#view
head( datadf ); dim( datadf )
# Then using StartDay to filter records, removing those that occurred#
#  earlier when unit was turned on, but not fitted to animal #
datadf <- datadf %>% 
  group_by( serial ) %>% 
  dplyr::filter( jday > StartDay ) %>% ungroup() # why group and ungroup?
#view
head( datadf ); dim( datadf )
# serial IDs are cumbersome so we create a new individual ID column:
datadf$id <- group_indices( datadf, serial )

datadf <- datadf %>%
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
# Maybe an easier way to do this?

unique(datadf$territory)
colnames(datadf)

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
table( datadf$territory )
# How many individuals have we dropped so far?
# 1
# We can also get an idea of the data collection for each individual
# by plotting histograms
ggplot( datadf, aes( x = jday, group = territory ) ) +
  theme_classic( base_size = 15 ) +
  geom_histogram( ) +
  facet_wrap( ~ territory )

# What do the histograms tell you about the nature of the data #
# Sample size, intensity for different individuals? #
# Answer:
#

#amt requires us to turn data into tracks for further analyses.
trks <- datadf %>% 
  #make track. Note you can add additional columns to it
  amt::make_track(.y = lat, .x = lon, .t = ts, 
                  #define columns that you want to keep, relabel if you need:
                  id = id, sex = Sex, mth = mth,jday = jday, speed = speed, alt = alt, territory = territory, 
                  #assign correct crs
                  crs = crsdata )

# Reproject to UTM to convert lat lon to easting northing:
trks <- amt::transform_coords( trks, crstracks )
#Turn into a tibble list by grouping and nest by individual IDs:
trks <- trks %>%  amt::nest( data = -"territory" )
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
# Answer: they all have points outside the NCA
#
# Any ideas on how to remove migration data?
# Answer:
# 
# Here we rely on NCA polygon, removing records that exist East of the #
# NCA. We can extra the extent of a polygon:
sf::st_bbox(NCA_Shape)
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
    mutate( jday = lubridate::yday( t1_ ) ) %>% #what is t1_?
    group_by( jday ) %>% #why group by jday?
    summarise( sl_ = log( sum(sl_) ) ) %>% #what's happening here? sl = step lengths
    ggplot(.) + theme_bw(base_size = 17) +
    labs( title = paste0('individual =', trks$territory[i]) ) +
    geom_line( aes( y = sl_, x = jday))
  print(a)
}

# We focus on breeding season data:
# Estimate sampling rate for each individual by looping through 
# data using purrr function map
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
par( mfrow = c( 2,3 ) ) # what's this?
for( i in 1:dim(trks.all)[1] ){
  #extract individual ids
  idd <- trks.all$territory[i]
  #use tibbles we calculated in steps
  x <- pull( trks.all[["steps"]][[i]], direction_p )
  #remove missing data
  x <- x[!is.na(x)]
  #calculate autocorrelation function:
  acf( x, lag.max = 300,
       main = paste0( "individual = ",idd ) )
  #Note you can modify the lag.max according to your data # lag is in minutes?
}
# What is ACF ??????????????????????????????????????????????????????????????????
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
# difference between rate and tolerance?????????????????????????????????????????

# We can now unnest the dataframes of interest
#Starting with all breeding season data
# What is nesting and unesting??????????????????????????????????????????????????
trks.breed <- trks.all %>% select( territory, breeding ) %>% 
  unnest( cols = breeding ) 
head( trks.breed )

# Now breeding season data, without autocorrelation:
trks.red <- trks.all %>% select( territory, red ) %>% 
  unnest( cols = red ) 
head( trks.red )

# Last all migration data:
trks.mig <- trks.all %>% select( territory, migrating ) %>% 
  unnest( cols = migrating ) 
head( trks.mig )
###############
#........................................ Lost after here.......................

###########################################################################
###### Creating tracks, calculating step lengths and turning angles ######
#                     for a single individual                             #
######################
# We start by creating a track for a single individual:
tr.idv <- datadf %>% 
  #select data for one individual only
  dplyr::filter( id == 1 ) %>% # tried to replace this line with filter( territory == SG ) but it didn't work?????
  # remove duplicate times
  dplyr::filter( !duplicated( ts ) ) %>% 
  #make track. Note you can add additional columns to it
  amt::make_track(.y = lat, .x = lon, .t = ts, id = id, mth = mth, 
                  jday = jday, sex = Sex, speed = speed, alt = alt, territory = territory, 
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
# What is it telling us? Can't remember
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
# What do these plot tell us? Not sure
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

#Reduce individual tracks based on selected steps above 
tr.red <- tr.slow %>% select( x_= x1_, y_=y1_, t_=t1_ ) %>% 
  left_join( tr.idv, by = c('x_', "y_", "t_" ) )
# view
tr.red
# Note that the output is a tibble. Turn it back to a track:
tr.red <- tr.red %>%  
  amt::make_track(.y = y_, .x = x_, .t = t_, id = id, 
                  mth = mth, jday = jday, speed = speed, alt = alt, territory = territory,
                  crs = crstracks )

#############################################################################
# Saving relevant objects and data ---------------------------------
#save hourly detection dataframe with weather predictors
# write.csv(x = det_df, 
#           #ensure that you save it onto your datafolder
#           file = paste0( datapath, 'stoc_det_df.csv'), 
#           row.names = FALSE )
#save workspace in case we need to make changes
save.image( "ER_TracksWorkspace.RData" )

########## end of save #########################
############### END OF Tracks ########################################

# load packages relevant to this script:
library( tidyverse ) #easy data manipulation
library( amt )
library( sp )
library( lubridate ) #easy date manipulation

#####################################################################
## end of package load ###############

###################################################################
#### Load or create data -----------------------------------------
# Clean your workspace to reset your R environment. #
rm( list = ls() )

# # Set working directory. This is the path to your Rstudio folder for this 
# # project. If you are in your correct Rstudio project then it should be:
# getwd()
# # if so then:
# workdir <- getwd()

# load workspace 
load( "ER_TracksWorkspace.RData" )
###############################################################
##### Comparing different estimators for occurrence      #######
###   distributions for all individuals at once.          ####
## We evaluate Minimum Convex Polygons (MCP),             ####
### Kernel density estimators (KDE) and autocorrelated    ####
### KDEs (AKDE). The latter does not need prior removal of ###
#  autocorrelated data.                                     ##
## Individuals are often sampled for different time periods ##
# so we also standardize time periods to evaluate the       ##
### effects of sampling period on home range estimates.     ##
##############################################################

# We start with thinned data, required for MCP and KDE methods:
hrred <- trks.red %>%  amt::nest( data = -"territory" ) %>% 
  mutate(
    hr_mcp = map(data, ~ hr_mcp(., levels = c(0.5, 0.9)) ),
    hr_kde = map(data, ~ hr_kde(., levels = c(0.5, 0.9)) ),
    n = map_int( data, nrow )
  )  
#view
hrred 

#plot home ranges
#select tibble 
hrred %>%
  #choose one home range method at a time
  hr_to_sf( hr_kde, territory, n ) %>% 
  #hr_to_sf( hr_mcp, territory, n ) %>% 
  #plot with ggplot
  ggplot( . ) +
  theme_bw( base_size = 17 ) + 
  geom_sf() +
  #plot separate for each indvidual
  facet_wrap( ~territory )

#We can see large variation of home range sizes between individuals#
# during the breeding season. BUT, we know our sampling wasn't #
# consistent. To account for our variable sampling, we can plot #
# estimated occurrence distributions weekly. #

#To do this, we first need to work out week of the year. #
trks.red <- trks.red %>%  
  mutate( wk = lubridate::week( t_) ) 

# recalculate n and homerange estimates
hr_wk <- trks.red %>%  
  # we nest by id and week
  nest( data = -c(territory, wk)) %>%
  mutate( n = map_int(data, nrow) ) %>% 
  #remove weeks without enough points
  filter( n > 15 ) %>% 
  mutate( #now recalculate weekly home range
    hr_mcp = map(data, ~ hr_mcp(., levels = c(0.5, 0.9)) ),
    hr_kde = map(data, ~ hr_kde(., levels = c(0.5, 0.9)) ))

hr_wk

# How many points are enough? 
# Answer: 
#
# How many weeks of data did you lose by removing weeks 
# without enough points?
# Answer:
# 
colnames(hr_wk)
head(ids)
unique(ids)

#plot weekly home ranges
#define a vector with individual ids
ids <- unique( hr_wk$territory )
# this way you can loop through each individual
for( i in 1:length(ids)){
  wp <- hr_wk %>% filter( territory == ids[i] ) %>% 
    #choose one home range method at a time
    hr_to_sf( hr_kde, wk, n ) %>% 
    #hr_to_sf( hr_mcp, wk, n ) %>% 
    #plot with ggplot
    ggplot( . ) +
    theme_bw( base_size = 15 ) + 
    geom_sf() +
    labs( main = ids[i] ) +
    #plot separate for each indvidual
    facet_wrap( ~wk )
  # prints each individual separately
  print( wp )
}
# Try plotting this with the other estimator. 
# what differences do you see between them?
# Answer: 
#
for( i in 1:length(ids)){
  wp <- hr_wk %>% filter( territory == ids[i] ) %>% 
    #choose one home range method at a time
    #hr_to_sf( hr_kde, wk, n ) %>% 
    hr_to_sf( hr_mcp, wk, n ) %>% 
    #plot with ggplot
    ggplot( . ) +
    theme_bw( base_size = 15 ) + 
    geom_sf() +
    labs( main = ids[i] ) +
    #plot separate for each indvidual
    facet_wrap( ~wk )
  # prints each individual separately
  print( wp )
}
# We can also calculate the weekly area. We take weekly home ranges#
# remove tracking data and convert to long dataframe
hr_area <- hr_wk %>%  select( -data ) %>% 
  pivot_longer( hr_mcp:hr_kde, names_to = "estimator", 
                values_to = "hr" )
#view
hr_area
# The we calculate area for each method 
hr_area <- hr_area %>%  
  mutate( hr_area = map( hr, ~hr_area(.)) ) %>% 
  unnest( cols = hr_area )
#convert area in m^2 to area in km^2 
hr_area$area_km <- hr_area$area / 1e6

#plot 
hr_area %>% 
  #choose desired level 
  filter( level  == 0.5 ) %>% 
  ggplot( aes(col = as.character(territory), y = area_km, x = wk )) + 
  geom_line(size = 1.5) + 
  geom_point(size = 4) +
  theme_light(base_size = 15) + 
  facet_wrap( ~estimator, nrow = 2, 
              scales = "free_y" )
# Comment on this graph # I'm pretty confused by this graph -- why are weeks not whole numbers?
# It seems nearly impossible that PR II had a week with a 600 km home range at 0.5, something might be wrong?
#Answer:
# 
# Try a different level and comment on new output
hr_area %>% 
  #choose desired level 
  filter( level  == 0.9 ) %>% 
  ggplot( aes(col = as.character(territory), y = area_km, x = wk )) + 
  geom_line(size = 1.5) + 
  geom_point(size = 4) +
  theme_light(base_size = 15) + 
  facet_wrap( ~estimator, nrow = 2, 
              scales = "free_y" )
# Answer: Something still feels very off about the area for PR II
#

# we can also calculate mean and CIs of area in Km for each id:
ci_wk <- hr_area %>% group_by(estimator, territory, level ) %>% 
  summarise( m = mean(area_km), 
             se = sd(area_km) / sqrt(n()), 
             me = qt(0.975, n() - 1) * se, 
             lci = m - me, uci = m + me)
#view
ci_wk %>% filter( level  == 0.9 )

# plot home ranges of each individual average across weeks:
ci_wk %>% 
  #choose desired level
  filter( level  == 0.9 ) %>% 
  ggplot(.) + 
  geom_pointrange(aes(x = as.character(territory), y = m, 
                      ymin = lci, ymax = uci, col = estimator), 
                  position = position_dodge2(width = 0.5)) +
  ylim( 0,2000 ) +
  theme_light(base_size = 15 )  
# Comment on output
# Answer:
#

# What other plots would be relevant for your specific question?
# Answer:
# Adapt to data here if possible. If not, provide details of what #
# is missing from this dataset to achieve the desired plot?
# Answer: 
# 

#######################################################################
############ Home range overlap ###############################
# amt currently implements methods reviewed by Fieberg & Kochany (2005) #
# hr: proportion of home range of instance i that overlaps with the home #
# range of instance j. This measure does not rely on a UD and is #
# directional (i.e., HRi,j≠HRj,i) and bound between 0 (no overlap) #
# and 1 (complete overlap) #
# phr: Is the probability of instance j being located in the home range #
# of instance i. phr is also directional and bounded between 0 (no #
# overlap) and 1 (complete overlap) #
# vi: The volumetric intersection between two UDs.#
# ba: The Bhattacharyya’s affinity between two UDs. #
# udoi: A UD overlap index. # 
# hd: Hellinger’s distance between two UDs. #


##############   Home range  for one individual ####################
#compare same methods for a single individual #
# using data with no autocorrelation again:
red_kde <- amt::hr_kde( tr.red, levels = c(0.5, 0.9) )
red_mcp <- amt::hr_mcp( tr.red, levels = c(0.5, 0.9) )

# Now estimate home range with a continuous-time movement model:
# options are "iid": for uncorrelated independent data, 
#  "bm": Brownian motion, "ou": Ornstein-Uhlenbeck process,
# "ouf": Ornstein-Uhlenbeck forage process, 
# "auto": uses model selection with AICc to find bets model
red_akde <- amt::hr_akde( x = tr.red, 
                          model = fit_ctmm( tr.idv, "ouf" ),
                          levels = c(0.5, 0.9))

# This takes foreverrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
# Also, didn't really work based on the errors in subsequent lines
# Comment on the choice of model you used for the akde?
# Answer:
#

# Plot output of different methods, including locations:
plot( red_kde, #add.relocations = FALSE, 
      lwd = 2, col = 'red' )
plot( red_akde,  add = TRUE, #add.relocations = TRUE,
      col = 'orange', lwd = 2 )
plot( red_mcp, add.relocations = TRUE, add = TRUE, 
      col = 'blue', lwd = 2 )

plot(red_akde) # error in plot

#Estimate areas for each method
amt::hr_area( red_kde )
#amt::hr_area( red_akde ); 
amt::hr_area( red_mcp )
# comment on the results
# Answer:
#

# Estimate home ranges but now use data that has not been thinned:
#Kernel 
idv_kde <- amt::hr_kde( tr.idv, levels = c(0.5, 0.9) )
#Autocorrelated kernel
red_akde <- amt::hr_akde( x = tr.idv, 
                          model = fit_ctmm( tr.idv, "ouf" ),
                          levels = c(0.5, 0.9))
#Minimum convex polygon
idv_mcp <- amt::hr_mcp( tr.idv, levels = c(0.5, 0.9) )

#We can also estimate corresponding areas
amt::hr_area( idv_kde );amt::hr_area( red_akde ); amt::hr_area( idv_mcp )

#view
par( mfrow = c(1,1))
plot( idv_kde )
plot( idv_akde, add.relocations = FALSE, add = TRUE, lty = 3  )
plot( idv_mcp, add.relocations = FALSE, add = TRUE, lty = 2 )
# What was the influence of removing temporal autocorrelation on 
# results for each method?
# Answer:
#

###########################################################
### Save desired results                                  #
save.image( "ER_homerangeresults.RData" )
############# end of script  ###########################################











