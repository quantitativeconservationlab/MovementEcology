#################################################################
# Script developed by Jen Cruz to format cleaned GPS points of 5sec   #
# resolution to even sample rates and to add available points for      #
# continuous tracks to then estimate #
# habitat selection using SSF or  iSSFs          # 
# We rely heavily on amt getting started vignette here:       #
# https://cran.r-project.org/web/packages/amt/vignettes/p1_getting_started.html#
#                                                               #
# Use locations of Prairie Falcon were collected during Spring/Summer #
# of 2021 at Morley Nelson Birds of Prey NCA.                      #
# Data were collected for 10 individuals and multiple #
# fix rates. We resample the 5 secs rate collected 3 times a week #

# Predictors are vegetation cover metrics from Rangeland Analysis Platform #
# https://rangelands.app/products/  #
# for western USA from satellite imagery processed at 30 x 30 m resolution #

# Habitat types including cropland and water  can also be obtained here:#
# https://www.mrlc.gov/data     #

# Elevation metrics using elevatr: #
# https://cran.r-project.org/web/packages/elevatr/vignettes/introduction_to_elevatr.html

# we use tmap for spatial data visualization. Learn more at:
# https://r-tmap.github.io/tmap-book/index.html #

# Note that this code takes a lot of time to run #
##################################################################

################## Prep. workspace ###############################
#install packages
install.packages( "tmap") 
install.packages( "corrplot" )
install.packages( "elevatr")
# load packages relevant to this script:
library( tidyverse ) #easy data manipulation and plotting
# set option to see all columns and more than 10 rows
options( dplyr.width = Inf, dplyr.print_min = 100 )
library( amt ) #handling tracks from location data
library( sf ) #handling vector data
library( raster ) #handle raster data
library( tmap ) #visualize raster and vector data together
library( corrplot )
#library( elevatr ) #to fetch elevation data
library( RColorBrewer)
## end of package load ###############

###################################################################
#### Load or create data -----------------------------------------
# Clean your workspace to reset your R environment. #
rm( list = ls() )

# If this is not the first time working on this script load workspace
# to pick up where you left off
#load( "DataPrep_5secSteps.RData" )

#if you are starting new then load your data:
#import polygon of the NCA as sf spatial file:
NCA_Shape <- sf::st_read(  "Data/BOPNCA_Boundary.shp" )

#import nest locations
nests <- sf::st_read(  "Data/nests.shp" )
#We cropped the vegetation cover raster to our study area so that #
# the image could be shared via github, which has size restrictions #
# Load the cropped raster here:
cover_NCA <- raster::stack( "Data/RAPcover2021_NCA.img" )

#import akde ranges you created which includes  thinned (30min) data
#akde_all <- read_rds( "Data/akde_all" )

#load high resolution data for comparison
trks.breed <- read_rds( "Data/trks.breed" )
#for comparison
trks.thin <- read_rds( "Data/trks.thin" )

#########
#####################################################################
###    select tracks that  are long enough to provide meaningful #
# paths...also check that no 30 min single points were left behind #
#######################################################################
head( trks.breed)
#create reference dataframe
iddf <- trks.breed %>% 
  group_by( id ) %>% 
  dplyr::select( id, territory, sex ) %>% 
  slice(1)

iddf

#create nest buffer
nest_buffer <- nests %>% 
  dplyr::select( territory= terrtry ) %>% 
  st_buffer(750) %>% 
  right_join( iddf, by = "territory" )

#view
head(nest_buffer)

#extract territory IDs
terids <- unique( nest_buffer$territory )
#add rowid 
trks.5sec <- trks.breed %>% 
  dplyr::mutate( rowid = row_number())
#create a points only sf
breed_sf <- trks.5sec %>% 
  dplyr::select( territory, rowid, x_, y_ ) %>% 
  amt::as_sf_points()
# I think the calculation needs to be individual specific
#try for 1
buf <-  nest_buffer %>% dplyr::filter( territory == "SG" )
bpnts <- breed_sf %>% dplyr::filter( territory == "SG" )
b <- st_difference( bpnts, buf )
forage.breed <- b
#now loop through the rest of the individuals
for( i in terids[2:length(terids)] ){
  buf <-  nest_buffer %>% dplyr::filter( territory == i )
  bpnts <- breed_sf %>% dplyr::filter( territory == i )
  #this function keeps nonoverlapping points only
  b <- st_difference( bpnts, buf )
  #append updated points to original
  forage.breed <- bind_rows(forage.breed, b )
}

trks.5sec <- trks.5sec %>% 
  dplyr::filter( rowid %in% forage.breed$rowid )
#the burst id restarts for each individual
# we need to create indiv specific ones
#add counts 
trks.5sec <- trks.5sec %>% 
  mutate( burstid = paste(burst_, id, sep = "_")) %>% 
  add_count( burstid )

#check 
head(trks.5sec)

# #we potentially limit to 10 steps
# trks.5sec <- trks.5sec %>% 
#   dplyr::filter( n > 19 )

#check loses
dim(trks.breed); dim(trks.5sec)
### plot to see what they look like
table(trks.5sec$id)

#we start by joining 5sec fixes for each individual with straight lines
lines5secsteps <- trks.5sec %>%
  group_by( id, burstid ) %>%
  as_sf_points() %>%
  summarise( do_union = FALSE ) %>%
  st_cast("LINESTRING" ) 

#then we plot 30min lines vs 5 sec points for each individual at a time
for( i in unique(trks.5sec$id) ){
  tp<-  
ggplot() +
    theme_bw( base_size = 15 ) + 
    theme( legend.position = "none" ) +
     geom_sf( data = as_sf_points( trks.thin ) %>% filter( id == i ), 
              color = "black", size = 0.5 ) +
    geom_sf( data = lines5secsteps %>% filter( id == i ), 
             aes( color = as.factor(burstid) ),
             linewidth = 1.5 ) +
    # geom_sf(data = NCA_Shape, linewidth = 1.5,
    #         inherit.aes = FALSE, fill=NA ) +
    labs( title = i )
   print(tp)
 }

###put them all in the same map
ggplot() +
  theme_bw( base_size = 15 ) + 
  theme( legend.position = "bottom" ) +
  geom_sf( data = lines5secsteps, 
           aes( color = as.factor(id) ),
           linewidth = 1.2 ) +
  geom_sf( data = as_sf_points( trks.thin ), 
           aes(color = as.factor(id)), size = 0.5 ) +
  geom_sf(data = NCA_Shape, linewidth = 1.5,
          inherit.aes = FALSE, fill=NA ) 
  

#### try to get speed and altitude
head(trks.5sec)
# for( i in unique(trks.5sec$id) ){
#   tp<-  
i = 8
    ggplot() +
    theme_bw( base_size = 15 ) + 
    theme( legend.position = "right" ) +
    geom_sf( data = as_sf_points( trks.5sec) %>% filter( id == i ), 
             aes( color = alt ) , 
            #aes( color = speed ) , 
                 size = 0.5 ) +
      scale_color_gradient(low="blue", high="red")
    # geom_sf( data = lines5secsteps %>% filter( id == i ), 
    #           color = "black",linewidth = 1 ) #+
    # geom_sf(data = NCA_Shape, linewidth = 1.5,
    #         inherit.aes = FALSE, fill=NA ) +
#     labs( title = i )
#   print(tp)
# }

###################################################################
####  resample to coarser resolutions to see changes         ######
######################################################################
#########
#### 20 sec resampled #####
trks.20sec <- trks.5sec %>%   
  amt::nest( data = -"id" ) %>% 
  mutate(
    xxx = map( data, function(x) x %>% 
                    track_resample( rate = seconds(20), 
                                    tolerance = seconds(5) )) )%>% 
  dplyr::select( id, xxx ) %>% 
  unnest( cols = xxx ) 

#view
head(trks.20sec)
#add burst ID
trks.20sec <- trks.20sec %>% 
  dplyr::select( -n ) %>% 
  add_count( burstid ) 
#
table(trks.20sec$burstid)

#remove steps < 3
trks.20sec <- trks.20sec %>% 
  group_by( burstid ) %>% 
  dplyr::filter( n > 9 ) %>% ungroup()

#we start by joining 2min fixes for each individual with straight lines
lines20sec <- trks.20sec %>%
  group_by( id, burstid ) %>%
  as_sf_points() %>%
  summarise( do_union = FALSE ) %>%
  st_cast("LINESTRING" ) 

head(lines20sec)
#then we plot 30min lines vs 5 sec points for each individual at a time
for( i in unique(trks.5sec$id) ){
    tp<-  
#i = 8
ggplot() +
  theme_bw( base_size = 15 ) + 
  theme( legend.position = "none" ) +
  geom_sf( data = as_sf_points( trks.thin ) %>% dplyr::filter( id == i ), 
           color = "black", size = 0.5 ) +
 # geom_sf( data = lines1min %>% 
    geom_sf( data = lines20sec %>% 
             dplyr::filter( id == i ), 
           aes( color = as.factor(burstid) ),
           linewidth = 1 ) +
  geom_sf( data = nest_buffer %>% dplyr::filter( id == i ), 
           color = "black", fill = NA ) +
  # geom_sf(data = NCA_Shape, linewidth = 1.5,
  #         inherit.aes = FALSE, fill=NA ) +
  labs( title = i )
   print(tp)
 }

######
################################################################
###calculate steps ####
steps20sec <- trks.20sec %>% 
  steps_by_burst( keep_cols = 'start' ) 


head(steps20sec)
table(steps20sec$n)
lubridate::minute( steps20sec$t1_)

a <- steps20sec %>% 
  group_by( burstid ) %>% 
  dplyr::mutate( minmins = min(lubridate::minute( t1_)),
                 mins =  lubridate::minute( t1_) - minmins,
                 hrs  =  lubridate::hour( t1_) ) %>% 
  ungroup() %>% 
  dplyr::select( rowid, mins, hrs )

steps20sec <- left_join( steps20sec, a, by = "rowid" )
# for( i in unique(trks.20sec$id) ){
#    tp <-  
i = 6
steps20sec %>% as_sf_points() %>% dplyr::filter( id == i ) %>% 
    ggplot(.) +
    theme_bw( base_size = 5 ) + 
    theme( legend.position = "right" ) +
  #   geom_sf( aes( color = direction_p ) ) +
  # scale_color_gradient2(low="blue", mid = "grey", high="red") +
geom_sf( aes( color = hrs ) ) +
  scale_color_gradient(low="blue", high="orange") +
 geom_sf( data = nest_buffer %>% dplyr::filter( id == i ),
           color = "black", fill = NA ) +
facet_wrap( ~burstid )
    # geom_sf(data = NCA_Shape, linewidth = 1.5,
    #         inherit.aes = FALSE, fill=NA ) +
     labs( title = i )
 #   print(tp)
 # }

     
### plot time of day that the falcons are moving
ggplot( steps20sec ) +
  theme_bw( base_size = 10 ) + 
  geom_density( aes( hrs, color = sex ) , alpha = 0.5 ) +
  facet_wrap(~territory)

###########################################################
### Save desired results #
#save lines 20secs
write_rds( lines20sec, "Data/lines20sec" )
write_rds( steps20sec, "Data/steps20sec" )
#save recalculated steps at 

#save workspace if in progress
save.image( 'DataPrep_5secSteps.RData' )
############# end of script  ##################################