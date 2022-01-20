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
load( "TracksWorkspace.RData" )
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
hrred <- trks.red %>%  amt::nest( data = -"id" ) %>% 
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
  hr_to_sf( hr_kde, id, n ) %>% 
  #hr_to_sf( hr_mcp, id, n ) %>% 
  #plot with ggplot
  ggplot( . ) +
  theme_bw( base_size = 17 ) + 
  geom_sf() +
  #plot separate for each indvidual
  facet_wrap( ~id )

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
      nest( data = -c(id, wk)) %>%
      mutate( n = map_int(data, nrow) ) %>% 
  #remove weeks without enough points
  filter( n > 15 ) %>% 
  mutate( #now recalculate weekly home range
  hr_mcp = map(data, ~ hr_mcp(., levels = c(0.5, 0.9)) ),
  hr_kde = map(data, ~ hr_kde(., levels = c(0.5, 0.9)) ))

# How many points are enough? 
# Answer: 
#
# How many weeks of data did you loose by removing weeks 
# without enough points?
# Answer:
# 

#plot weekly home ranges
#define a vector with individual ids
ids <- unique( hr_wk$id )
# this way you can loop through each individual
#for( i in 1:2){#length(ids)){
  wp1 <- hr_wk %>% filter( id == ids[i] ) %>% 
    #choose one home range method at a time
    hr_to_sf( hr_kde, id, wk, n ) %>% 
    filter( level == 0.9 )
  wp <- hr_wk %>% filter( id == ids[i] ) %>% 
  #choose one home range method at a time
  hr_to_sf( hr_mcp, id, wk, n ) %>% 
    filter( level == 0.9 ) %>% 
  #hr_to_sf( hr_mcp, wk, n ) %>% 
  #plot with ggplot
  ggplot( . ) +
  theme_bw( base_size = 15 ) + 
  geom_sf(aes( fill = as.factor(wk) ) ) +
  geom_sf( data = wp1, aes( colour = as.factor(wk)), alpha = 0 ) +
    scale_x_continuous( breaks = c( -116.0,-115.8, -115.6 ) ) +  
  labs( title = ids[i], fill = "week", x = "lat") +
  #plot separate for each indvidual
  facet_wrap( ~wk )
 # prints each individual separately
  print( wp )
}

# Try plotting this with the other estimator. 
# what differences do you see between them?
# Answer: 
#

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

# Comment on the choice of model you used for the akde?
# Answer:
#

# Plot output of different methods, including locations:
plot( red_kde, add = TRUE, #add.relocations = FALSE, 
      lwd = 2, col = 'red' )
plot( red_akde,  add = TRUE, #add.relocations = TRUE,
      col = 'orange', lwd = 2 )
plot( red_mcp, add.relocations = TRUE, add = TRUE, 
      col = 'blue', lwd = 2 )

#Estimate areas for each method
amt::hr_area( red_kde ); amt::hr_area( red_akde ); amt::hr_area( red_mcp )
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
amt::hr_area( idv_kde );amt::hr_area( idv_akde ); amt::hr_area( idv_mcp )

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
save.image( 'homerangeresults.RData' )
############# end of script  ###########################################