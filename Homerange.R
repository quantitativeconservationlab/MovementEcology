##################################################################
# Script developed by Jen Cruz to calculate home ranges          # 
# We rely on amt vignette here:                                  #
# https://cran.r-project.org/web/packages/amt/vignettes/p2_hr.html #
# as well as: Signer & Fieberg (2021) PeerJ9:e11031              #
# http://doi.org/10.7717/peerj.11031                             #
###################################################################

################## prep workspace ###############################7

# load packages relevant to this script:
library( tidyverse ) #easy data manipulation
# set option to see all columns and more than 10 rows
options( dplyr.width = Inf, dplyr.print_min = 100 )
library( amt )
#library( sf )

#####################################################################
## end of package load ###############

###################################################################
#### Load or create data -----------------------------------------
# Clean your workspace to reset your R environment. #
rm( list = ls() )

# load workspace 
#load( "homerangeresults.RData" )

#load the thinned (30min) data
trks.thin <- read_rds( "Data/trks.thin" )

#load 5 sec data for comparison
trks.breed <- read_rds( "Data/trks.breed" )

#import polygon of the NCA as sf spatial file:
NCA_Shape <- sf::st_read( "Data/BOPNCA_Boundary.shp" )

###############################################################
##### Comparing different estimators for occurrence      #######
###   distributions for all individuals at once.          ####
## We evaluate Minimum Convex Polygons (MCP),             ####
### Kernel density estimators (KDE).                         ##
## Individuals are often sampled for different time periods ##
# so we also standardize time periods to evaluate the       ##
### effects of sampling period on home range estimates.     ##
##############################################################
#check if object has coordinate system in the right format
get_crs( trks.thin )
#view data
head( trks.thin, 10)
tail(trks.breed )
# MCP and KDE rely on data with no autocorrelation. But
# how do we check for autocorrelation to determine if # 
# we need to thin data further?
# One approach is by calculating autocorrelation functions
# We do so for each individual using our 30min sampled data:
par( mfrow = c( 2,3 ) )
#based on direction
for( i in 1:dim(trks.thin)[1] ){
  #select x location of each individual
  x <- trks.thin %>% dplyr::filter( id == i ) %>% 
    #select( x_ )
    select( y_ )
    #select(speed )
  #calculate autocorrelation function:
  acf( x, lag.max = 40,
       main = paste0( "individual = ", i ) )
  #Note you can modify the lag.max according to your data 
}

# Are our data autocorrelated?
# Answer:
# 

#if you are struggling to determine that maybe we can compare
# against our 5 second data, which we know will be correlated:
par( mfrow = c( 2,3 ) )
for( i in 1:dim(trks.breed)[1] ){
  #select x location of each individual
  x <- trks.breed %>% dplyr::filter( id == i ) %>% 
    #select( x_ )
    select( y_ )
  #select(speed )
  #select( alt )
  #calculate autocorrelation function:
  acf( x, lag.max = 500,
       main = paste0( "individual = ", i ) )
  #Note we modified the lag.max 
}
#What do you see regarding autocorrelation?
#Answer:
# 

# How much would you have to thin your data to remove autocorrelation?
# Answer:
#

# we keep the 30min interval (as we will make use of autocorrelation)
# next week and calculate MCP and KDE for each individual:
ranges <- trks.thin %>% 
  #we group tibbles for each individual:
  nest( data = -"id" ) %>% 
  #then add estimates from two home range measures:
  mutate(
    #Minimum Convex Polygon
    hr_mcp = map(data, ~ hr_mcp(., levels = c(0.5, 0.95)) ),
    #Kernel density estimator
    hr_kde = map(data, ~ hr_kde(., levels = c(0.5, 0.95)) ),
    #also calculate the sample size for each individual
    n = map_int( data, nrow )
  )  
#view
ranges

#plot KDEs:
ranges %>%
  hr_to_sf( hr_kde, id, n ) %>% 
  ggplot( . ) +
  theme_bw( base_size = 17 ) + 
  geom_sf( aes(color= as.factor(id) ) )  +
  geom_sf(data = NCA_Shape, inherit.aes = FALSE, 
          fill = NA ) +
  theme( legend.position = "none" ) +
  facet_wrap( ~id )

#save the MCP only so that we can plot it against the KDE:
mcps <- ranges %>%
  hr_to_sf( hr_mcp, id, n )

#plot both methods:
#select tibble 
ranges %>%
  #choose one home range method at a time
  hr_to_sf( hr_kde, id, n ) %>% 
  #plot with ggplot
  ggplot( . ) +
  theme_bw( base_size = 17 ) + 
  geom_sf( aes( fill = as.factor(id)), 
           linewidth = 0.8, alpha = 0.6 ) +
  geom_sf( data = mcps, colour = "black",
         linewidth = 1, fill = NA ) +
  geom_sf(data = NCA_Shape, inherit.aes = FALSE, fill=NA ) +
  theme( legend.position = "none" ) +
  #plot separate for each individual
  facet_wrap( ~id )

#which individuals show the biggest variation between the two methods?
#Answer:
#

#Is there evidence that sample size is affecting differences among 
#individuals? How are you deciding your answer?
#Answer:
#

# Using KDE we can see large variation among individuals#
# during the breeding season. Is it real? We know our sampling wasn't #
# consistent. To account for our variable sampling, we can plot #
# estimated occurrence distributions weekly. #

# recalculate n and weekly range estimates for kdes:
hr_wk <- trks.thin %>%  
      # we nest by id and week
      nest( data = -c(id, wk, sex ) ) %>%
      mutate( n = map_int(data, nrow ) ) %>% 
  #remove weeks without enough points
  filter( n > 15 ) %>% 
  mutate( #now recalculate weekly home range
  hr_kde = map(data, ~ hr_kde(., levels = c(0.5, 0.95)) ))

#plot weekly ranges for each individual at a time so that
# we can focus on within-individual differences
#define a vector with individual ids to loop through
ids <- hr_wk %>% 
  group_by( id ) %>% 
  slice( 1 ) %>% 
  dplyr::select( id, sex )

# this way you can loop through each individual
for( i in ids$id ){
# #filter data for one individual at a time:
wp <- hr_wk %>% filter( id == ids$id[i] ) %>% 
#turn range estimates to sf objects, keeping relevant id details
hr_to_sf( hr_kde, id, sex, wk, n ) %>% 
  #plot with ggplot
  ggplot( . ) +
  theme_bw( base_size = 13 ) + 
  # color code weekly ranges by week
  geom_sf( aes( fill = as.factor(wk),color = as.factor(wk) ), 
           alpha = 0.3, linewidth = 1 ) +
  #add used locations from 5 sec data as a check:
  geom_sf( data = as_sf_points( trks.breed %>% 
                  filter( id == ids$id[i] ) ),
           size = 0.5 ) +
  #add NCA polygon
#  geom_sf(data = NCA_Shape, inherit.aes = FALSE, fill=NA ) +
  #add labels including individual ID and sex as title
  labs( title = paste( ids$id[i], ids$sex[i] ),
    fill = "week", color = "week", x = "lat" ) + 
  #choose legend location
  theme( legend.position = "bottom" ) +
  #plot separate for each week
  facet_wrap( ~wk )
 # prints each individual separately
   print( wp )
}

# to help us work out where we are in the season we extract initial 
# dates for each week
trks.thin %>% group_by(wk ) %>% 
  slice( 1 )
#How are the range distributions changing on a weekly basis?
# Answer:
#

#Let's compare how our 30m and 5 sec tracks compare tracks
ggplot( trks.breed, aes( x = x_, y = y_ ) ) +
  theme_bw( base_size = 15 ) +
  geom_point( size = 0.5 ) +
  geom_point( data = trks.thin, color = "red",
           size = 1 ) +
  facet_wrap( ~id, scales = "free" )

###### Estimating home range area ##########################
# We can also calculate the range size area.
hr_area <- ranges %>%  select( -data ) %>% 
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

#add sex attribute 
hr_area <- left_join( hr_area, ids, by = "id" )
#check
head(hr_area)
#plot 
hr_area %>% 
  #choose desired level 
  filter( level  == 0.95 ) %>% 
  ggplot( aes(x = as.character(id), y = area_km, 
              color = sex ) ) + 
  geom_point(size = 4) +
  theme_light(base_size = 15) + 
  facet_wrap( ~estimator, nrow = 2, 
              scales = "free_y" )
# Comment on this graph
#Answer:
# 
# Replot this with 50% area size.
# Add Code:
#
# How do individuals and sexes differ between 95 and 50% 
# breeding ranges? # what could you tentatively say about 
# their ecology based on these results # 
# Answer:
#

# optional question:
# Adapt area size code to plot weekly changes in area size 
# for the weekly ranges you subsetted 
# Code Answer: 
# 

#### end of area ##############

###########################################################
### Save desired results                                  #
#if you are still getting through this script then save 
# the workspace here so you don't have to rerun your code
save.image( 'homerangeresults.RData' )

# Once you are finished:
#save breeding ranges
write_rds( ranges, "Data/ranges" )
#save weekly home ranges
write_rds( hr_wk, "Data/hr_wk" )

# save the homework plots and upload them to github 
# Here:
# 

############# end of script  ###########################################