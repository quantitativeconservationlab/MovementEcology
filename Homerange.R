##################################################################
# Script developed by Jen Cruz to calculate home ranges          # 
# We rely on amt vignette here:                                  #
# https://cran.r-project.org/web/packages/amt/vignettes/p2_hr.html #
# as well as: Signer & Fieberg (2021) PeerJ9:e11031              #
# http://doi.org/10.7717/peerj.11031                             #
###################################################################

################## prep workspace ###############################


# load packages relevant to this script:
library( tidyverse ) #easy data manipulation
# set option to see all columns and more than 10 rows
options( dplyr.width = Inf, dplyr.print_min = 100 )
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
#load( "homerangeresults.RData" )

#load the thinned (30min) data
trks.thin <- read_rds( "Data/trks.thin" )

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
head( trks.thin)

# MCP and KDE rely on data with no autocorrelation. But
# how do we check for autocorrelation to determine if # 
# we need to thin data further?
# One approach is by calculating autocorrelation functions
# We do so for each individual using our 30min sampled data:
par( mfrow = c( 2,3 ) )
#based on direction
for( i in 1:dim(trks.thin)[1] ){
  #select x location of each individual
  x <- trks.thin %>% filter( id == i ) %>% 
    #select( x_ )
    #select( y_ )
    #select(speed)
    select( alt )
  #calculate autocorrelation function:
  acf( x, lag.max = 1000,
       main = paste0( "individual = ", i ) )
  #Note you can modify the lag.max according to your data 
}

# Are our data autocorrelated?
# Answer:
# 

# We start by calculating MCP and KDE for each individual:
ranges <- trks.thin %>% 
  #we group tibbles for each individual:
  nest( data = -"id" ) %>% 
  #then add estimates from two home range measures:
  mutate(
    #Minimum Convex Polygon
    hr_mcp = map(data, ~ hr_mcp(., levels = c(0.5, 0.9)) ),
    #Kernel density estimator
    hr_kde = map(data, ~ hr_kde(., levels = c(0.5, 0.9)) ),
    #also calculate the sample size for each individual
    n = map_int( data, nrow )
  )  
#view
ranges

#plot KDEs:
#select tibble 
ranges %>%
  #choose one home range method at a time
  hr_to_sf( hr_kde, id, n ) %>% 
  #hr_to_sf( hr_mcp, id, n ) %>% 
  #plot with ggplot
  ggplot( . ) +
  theme_bw( base_size = 17 ) + 
  geom_sf() +
  #plot separate for each individual
  facet_wrap( ~id )

#plot MCPs:
ranges %>%
  hr_to_sf( hr_mcp, id, n ) %>% 
  ggplot( . ) +
  theme_bw( base_size = 17 ) + 
  geom_sf() +
  facet_wrap( ~id )

#We can see large variation of home range sizes between individuals#
# during the breeding season. Is it real? We know our sampling wasn't #
# consistent. To account for our variable sampling, we can plot #
# estimated occurrence distributions weekly. #

# recalculate n and homerange estimates
hr_wk <- trks.thin %>%  
      # we nest by id and week
      nest( data = -c(id, wk) ) %>%
      mutate( n = map_int(data, nrow) ) %>% 
  #remove weeks without enough points
  filter( n > 10 ) %>% 
  mutate( #now recalculate weekly home range
  hr_mcp = map(data, ~ hr_mcp(., levels = c(0.5, 0.9)) ),
  hr_kde = map(data, ~ hr_kde(., levels = c(0.5, 0.9)) ))

#plot weekly home ranges
#define a vector with individual ids
ids <- unique( hr_wk$id )
# this way you can loop through each individual
for( i in 1:length(ids)){
  #convert point locations to sf object for each animal
  points <- as_sf_points( trks.thin %>% filter( id == ids[i] ))
  #extract data for one home range method at a time:
  #mcp
  wp1 <- hr_wk %>% filter( id == ids[i] ) %>% 
    hr_to_sf( hr_kde, id, wk, n ) %>% 
    filter( level == 0.9 )
  #now for mcp
  wp <- hr_wk %>% filter( id == ids[i] ) %>% 
  hr_to_sf( hr_mcp, id, wk, n ) %>% 
    filter( level == 0.9 ) %>% 
  #plot with ggplot
  ggplot( . ) +
  theme_bw( base_size = 15 ) + 
  # start with mcp which is part of the piping
  geom_sf(aes( fill = as.factor(wk) ) ) +
  # to add kde then define separate data object
  geom_sf( data = wp1, aes( colour = as.factor(wk)), 
           alpha = 0, size = 2 ) +
    #scale_x_continuous( breaks = c( -117.0,-116.0, -115.0 ) ) +  
  #add used locations:
  geom_sf( data = points ) +
  labs( title = ids[i], fill = "week", x = "lat") +
  #plot separate for each indvidual
  facet_wrap( ~wk )
 # prints each individual separately
  print( wp )
}

#plot tracks
ggplot( trks.thin, aes( x = x_, y = y_, color = as.factor(wk) ) ) +
  theme_bw( base_size = 15 ) +
  geom_point( size = 2 ) +
  facet_wrap( ~as.factor(id), scales = "free" )

###### Estimating home range area ##########################
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
  ggplot( aes(col = as.character(id), y = area_km, x = wk )) + 
  geom_line(size = 1.5) + 
  geom_point(size = 4) +
  theme_light(base_size = 15) + 
  facet_wrap( ~estimator, nrow = 2, 
              scales = "free_y" )
# Comment on this graph
#Answer:
# 
# Try a different level and comment on new output
# Answer:
#

# we can also calculate mean and CIs of area in Km for each id:
ci_wk <- hr_area %>% group_by(estimator, id, level ) %>% 
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
  geom_pointrange(aes(x = as.character(id), y = m, 
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

#### end of area ##############
##############   Home range  for one individual ####################

###########################################################
### Save desired results                                  #
#save breeding season data (not thinned)
write_rds( hr_wk, "hr_wk")
#save range area estimates
write_rds( ci_wk, "ci_wk" )
#save.image( 'homerangeresults.RData' )
############# end of script  ###########################################