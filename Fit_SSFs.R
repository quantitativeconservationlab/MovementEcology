##################################################################
# Script developed by Jen Cruz to estimate SSFs and iSSFs          #
# approach derived from Fieberg et al. 2021 DOI: 10.1111/1365-2656.13441 #
# using code from Appendices B and C                             #
# also Muff et al. 2019 DOI: 10.1111/1365-2656.13087             #
# code here:                                                     #  
# https://conservancy.umn.edu/handle/11299/204737                #
#                                                                #
# We use landcover data from the National Geospatial Data Asset  #
# https://www.mrlc.gov                                           #
# Habitat predictors include 2018 estimates of sagebrush cover   #
###################################################################

################## prep workspace ###############################

# Clean your workspace to reset your R environment. #
rm( list = ls() )

# load packages relevant to this script:
library( tidyverse ) #easy data manipulation
# set option to see all columns and more than 10 rows
options( dplyr.width = Inf, dplyr.print_min = 100 )
library( amt )
library( sf )
library( terra ) # for raster manipulation
library( raster )
library( rasterVis ) #for raster visualization (of raster:: objects)
library( glmmTMB ) # for analysis

#####################################################################
## end of package load ###############

###################################################################
#### Load or create data -----------------------------------------

#load 30m steps estimated for all individuals
trks.steps <- read_rds( "trks.steps30" )
# load points also so that we can combine data
trks <- read_rds( "trks.thin" )
#view
head( trks.steps )

head( trks )
# #load range for all individuals estimated using ctmm 
# load( "../ctmm_akde_w.rda" ) 
# #check that it loaded the object
# class( akde.w )

#import polygon of the NCA as sf spatial file:
NCA_Shape <- sf::st_read("Z:/Common/QCLData/Habitat/NCA/GIS_NCA_IDARNGpgsSampling/BOPNCA_Boundary.shp")

# #import sagebrush raster with raster:
sagebrush <- raster::raster( "Z:/Common/QCLData/Habitat/NLCD/Sage_2007_2018/nlcd_sage_2018_mos_rec_v1.img" )
#visualize with either rasterVis or terra::plot:
rasterVis::levelplot( sagebrush )
#view
sagebrush
#######################################################################
######## preparing data ###############################################
# we prepare the predictor data similarly to how we did it for RSFs:
#get coordinates from shapefile
crstracks <- sf::st_crs( NCA_Shape )
# checking outline of NCA
sf::st_bbox( NCA_Shape )
# We define available habitat as area of NCA with a small buffer #
# around it and draw points from it #
#create a buffer around the NCA using outline of NCA and sf package:
# we are more generous than with our RSF analyses
NCA_buf <- NCA_Shape %>% sf::st_buffer( dist =1e4 )
#create a version that matches coordinates of the predictor raster:
NCA_trans <- sf::st_transform( NCA_buf, st_crs( sagebrush ) ) 
#crop raster to buffered NCA:
sage_cropped <- raster::crop( sagebrush, NCA_trans )
# Now that we have cropped it to the appropriate area it should be faster #
# to process 
#view
sage_cropped
#values greater than 100 are empty so replace with missing
sage_cropped[ sage_cropped > 100 ] <- NA
#Plot sagebrush
rasterVis::levelplot( sage_cropped )

#extract individual ids
ids <- unique(trks$territory)
# to create random steps, we start by nesting our data using purr:
steps_all <- trks.steps %>% nest( data = -"id" )
#view
steps_all
#we then estimate random steps
steps_all <- steps_all %>% 
          dplyr::mutate( rnd = lapply( data, function(x){
            amt::random_steps( x ) } ) )
# Note that by default the random_steps() function fits a tentative #
# gamma distribution to the observed step lengths and a tentative #
# von Mises distribution to the observed turn angles. #
# It then generates random available points by sampling step-lengths #
# and turn angles from these fitted distributions and combining these #
# random steps with the starting locations associated with each observed #
# movement step. #

#now unnest the new dataframes to make sure they worked
stepsdf <- steps_all %>% dplyr::select( id, rnd ) %>% 
  unnest( cols = rnd ) 
#view
head( stepsdf );dim( stepsdf )

# Before extracting data we need to turn it to sf object and #
# change projection to match the raster.

# We start by turning it to sf object, assigning the correct projection
steps_sf <- sf::st_as_sf( stepsdf, coords = c("x2_", "y2_"), 
                          crs = crstracks )
# Note that we use the end of the step coordinates, since we want to 
# focus on habitat selection, not movement.
steps_sf

# We then transform the crs:
steps_trans <- sf::st_transform( steps_sf, st_crs(sage_cropped) )

#extracting with raster we can used the sf object directly, you also 
# have the choice to use a buffer around each point if you want to increase 
# your resolution:
sage_30m <- raster::extract( x = sage_cropped, steps_trans,
                             method = "simple" )

#check
sage_30m
# What proportion of our data are missing values
sum( is.na( sage_30m ))/ length( sage_30m )

# We append our predictor estimates to the original steps tibble:
df_all <- cbind( stepsdf, sage_30m )
# Scale predictors 
df_all$sage_30m <- scale( df_all$sage_30m )
#replace missing values with mean, which is 0 after they have been scaled
df_all$sage_30m[ is.na(df_all$sage_30m) ] <- 0
# we also assign weights to available points to be much greater than #
# used points
df_all$weight <- 1000 ^( 1 - as.integer(df_all$case_ ) )
#check
head( df_all )
#### end data prep #############
###########################################################################
##### analyse data  ##########

#amt function doesn't take random effects or weights so we move to a more
# flexible package. We move straight into a model #
# that includes random intercepts and slopes as well as a fixed #
# large variance for the random intercepts, as recommended by Muff #
# et al. 2019 and weights of 1000 for available points #

# we start by defining the model without running it, which let's us
# fit the large variance to the random ID intercepts
m1.struc <- glmmTMB( case_ ~ sage_30m +  
                   #define random effects
                   ( 1| step_id_ ) + 
                   ( 1| id ) + ( 0 + sage_30m | id ), 
                   family = poisson, data = df_all, 
                   weights = weight, doFit=FALSE ) 

# fix variance
m1.struc$parameters$theta[ 1 ] <- log( 1e3 ) 
# tell it not to change variance
m1.struc$mapArg <- list( theta = factor( c(NA, 1:2) ) )

#then fit the model
m1 <- glmmTMB::fitTMB( m1.struc )
summary( m1 )

###########################################################
##### end ######

######################################################################
######## how would we partition data into behavioral states prior to #
# analyses ###############################################################
# for the used data we start by exploring it a bit more to refresh #
# our memory:
# we combine track dataframe to use information we have calculated 
# previously 
trks_all <- trks %>% 
  #select columns of interest
  dplyr::select( id, territory, sex, mth, jday,
                 alt, speed, x2_ = x_, y2_ = y_, t2_ = t_  )

trks_all <- trks_all %>% 
  # append to steps
  right_join( trks.steps, by = c("id", "x2_", "y2_", "t2_" ) )

#check
head( trks_all )
# we add week for easier visualization and subsetting 
trks_all <- trks_all %>%  
  #remove those bursts that have too few points
  amt::filter_min_n_burst( min_n = 5 ) %>% 
  #add week and hour columns 
    mutate( wk = lubridate::week( t2_ ), 
          hr =  lubridate::hour( t2_ ) )  

# we remind ourselves about over which weeks our data overlap
ggplot( trks_all, aes( x = jday, fill = as.factor(wk) ) ) +
  theme_classic( base_size = 15 ) +
  geom_histogram( alpha = 0.8 ) +
  facet_wrap( ~ id )
# from this we can see that we have an uneven sample size that we 
# need to deal with 

# we visualise step lengths and turning angles for each individual 
trks_all %>%   
  ggplot(.) +
#  geom_density( aes( x = sl_, fill = as.factor(wk) ), alpha = 0.6 ) +
  geom_density( aes( x = ta_, fill = as.factor(wk) ), alpha = 0.6 ) +
#  geom_density( aes( x = hr, fill = as.factor(wk)), alpha = 0.6 ) +
#  geom_density( aes( x = speed, fill = as.factor(wk) ), alpha = 0.6 ) +
  #geom_histogram( aes( x = sl_, fill = as.factor(wk) ) ) +
  #xlab("Step length" ) + 
  #ylim( 0, 0.01 ) + 
  #xlim(0, 300 ) +
  theme_bw( base_size = 19 )  +
  theme( legend.position = "none" ) +
  facet_wrap( ~id, scales = 'free' )

# We plot alternatives selections for step lengths and turning angles
# to choose parameters that may help to coarsely remove nesting 
# locations and travelling locations

ids <- unique( trks_all$id )
ters <- unique( trks_all$territory )
for( i in ters ){  
  #pull observation that you think belong to foraging
  a <- trks_all %>% dplyr::filter( territory == i ) %>% 
    #filter( ta_ < -0.2 | ta_ > 0.2 ) %>% 
    dplyr::filter( sl_ < 5000 ) %>% 
    dplyr::filter( speed > 1 | speed < 30 ) %>% 
    #dplyr::filter( hr > 6 | hr < 20 ) %>% 
    sf::st_as_sf(., coords = c("x2_", "y2_"), crs = crstracks )
  #pull all observations for that individual
  b <- trks_all %>% filter( territory == i ) %>% 
    sf::st_as_sf(., coords = c("x2_", "y2_"), crs = crstracks )
#compare
  c  <- ggplot( b ) +
    theme_bw( base_size = 15 ) +
    labs( title = i ) +
    theme( legend.position = "none" ) +
    geom_sf() +
    geom_sf( data = NCA_Shape, inherit.aes = FALSE, alpha = 0 ) +
    geom_sf( data = a, aes(col = as.factor( burst_ ) ), 
                           inherit.aes = FALSE )# +
    #facet_wrap( ~wk )

  print( i )
  print( c )
}

# We can use these parameters to subset the data but we would have #
# to recalculate steps #

##########################################################################
### Save desired results                                  #
# Save the steps dataframe with extracted raster values so that I 
# don't have to recreate it when estimating issfs 
write_rds( df_all, "df_all" )

write_rds( trks_all, "trks_all" )

#Also save the unscaled raster values as a csv
write.csv( sage_30m, "sage_30m_steps.csv", row.names = FALSE )

#save workspace if in progress
save.image( 'SSF_results.RData' )
############# end of script  ##################################