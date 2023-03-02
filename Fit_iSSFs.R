##################################################################
# Script developed by Jen Cruz to estimate iSSFs          #
# Approach largely followed Fieberg et al. 2021 DOI: 10.1111/1365-2656.13441 #
# specifically Appendices B and C                             #
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

#install.packages( "circular" )
# load packages relevant to this script:
library( tidyverse ) #easy data manipulation
# set option to see all columns and more than 10 rows
options( dplyr.width = Inf, dplyr.print_min = 100 )
library( amt )
library( sf )
library( glmmTMB )
library( circular ) #for plotting von mises distribution
library( raster )
#####################################################################
## end of package load ###############

###################################################################
#### Load or create data -----------------------------------------

# #load steps at 5sec resolution
trks <- read_rds( "Data/trks.steps" )
# #check
head(trks)
#OR if wanting to work with 30min resolution, load the 
# dataframe that we created in Fit_SSFs.R script:
dfraw <- read_rds( "df_all" )
#check
class( dfraw )
#import polygon of the NCA as sf spatial file:
NCA_Shape <- sf::st_read("Z:/Common/QCLData/Habitat/NCA/GIS_NCA_IDARNGpgsSampling/BOPNCA_Boundary.shp")
# #import sagebrush raster with raster:
sagebrush <- raster::raster( "Z:/Common/QCLData/Habitat/NLCD_new/TimeSeriesCover/Sagebrush/Sagebrush_2009_2020/rcmap_sagebrush_2020.img")


#######################################################################
######## preparing data ###############################################
# # now tks
 head( trks ); dim(trks)
# We cannot have missing values for the predictors so we 
# need to remove steps with missing ta_ 
#step ids are not unique to individuals so we create a unique id:
trks$id_step <- paste0( trks$id, trks$burst_ )
head( df);dim(df)

#we use either trks or df_all depending on our preference
df <- trks
#
#which steps have missing ta values:
df$id_step[ which( is.na(df$ta_) ) ]
#how many don't:
length( df$id_step[ which( !is.na(df$ta_) ) ] )
dim(df)
df <- df[which( !is.na(df$ta_) ), ]
head(df)
#which steps are 0 length
unique( df$id_step[ which( df$sl_ == 0) ] )

#recheck sample size:
table( df$id)
head(df)

### reextract sagebrush if using 5 sec resolution #######
## otherwise you already have it #
#get coordinates from shapefile
crstracks <- sf::st_crs( NCA_Shape )
#create a buffer around the NCA using outline of NCA and sf package:
# we are more generous than with our RSF analyses
NCA_buf <- NCA_Shape %>% sf::st_buffer( dist =1e4 )
#create a version that matches coordinates of the predictor raster:
NCA_trans <- sf::st_transform( NCA_buf, st_crs( sagebrush ) ) 
#crop raster to buffered NCA:
sage_cropped <- raster::crop( sagebrush, NCA_trans )
#values greater than 100 are empty so replace with missing
sage_cropped[ sage_cropped > 100 ] <- NA
#Plot sagebrush
rasterVis::levelplot( sage_cropped )
terra::plot(sage_cropped)

# to create random steps, we start by nesting our data using purr:
steps_all <- df %>% nest( data = -"id" )
#view
steps_all
#we then estimate random steps
steps_all <- steps_all %>% 
  dplyr::mutate( rnd = lapply( data, function(x){
    amt::random_steps( x ) } ) )
#now unnest the new dataframes to make sure they worked
stepsdf <- steps_all %>% dplyr::select( id, rnd ) %>% 
  unnest( cols = rnd ) 
# We start by turning it to sf object, assigning the correct projection
steps_sf <- sf::st_as_sf( stepsdf, coords = c("x2_", "y2_"), 
                          crs = crstracks )
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
df <- cbind( stepsdf, sage_30m )
df$sage_30m <- scale( df$sage_30m )
df$sage_30m[ is.na(df$sage_30m) ] <- 0

# we also assign weights to available points to be much greater than #
# used points
df$weight <- 1000 ^( 1 - as.integer(df$case_ ) )

#check
head( df)
### end adding sagebrush and weight ####
# we check sample sizes for each individual
table( df$id)
# #select ids those with poor sample size
# remids <- c( 1,3,7)
# # use those to remove their data from dataframe:
# df <- df %>%  
#   dplyr::filter( !(id %in% remids) )
# #check 
# dim(df)

### end data prep #######
######## analysis for single individual using amt #################
# Estimate empirical distributions for individual of interest:
emp_d_sl <- df %>% dplyr::filter( id == 4 ) %>%
  dplyr::select( sl_ ) %>% 
  fit_distr( ., dist_name = "gamma" )
#view
emp_d_sl
#now for turning angle using a von mises, circular distribution:
emp_d_ta <- df %>% dplyr::filter( id == 4 ) %>%
  dplyr::select( ta_ ) %>% 
  fit_distr( ., dist_name = "vonmises" )
#view
emp_d_ta

#Fit iSSF model to same individual:
m_sg <- df %>% dplyr::filter( id == 4 ) %>% 
  fit_issf( case_ ~ sage_30m + sl_ + log(sl_) + cos(ta_) +
              strata( step_id_ ), model = TRUE )
# Note the cosine of the turn angle and the log of the 
#step-length allow us to adjust/refine the parameters of our #
# tentative step-length and turn-angle distributions after #
# fitting our integrated step-selection model. #

#view results
summary( m_sg )

#### interpret results ###
# One option that can help interpret results is by calculating #
# relative selection strength for two locations that are equal #
# except for 1 of the predictors. We demonstrate that approach here: #

# start by creating new dataframes that contain details for each #
# of the two hypothetical locations that we want to compare: 
# We focus on differences in sagebrush #
s1 <- data.frame( 
  sage_30m = 0, 
  sl_ = 5000,
  ta_ = 0 )
# now a second dataframe with higher sagebrush
s2 <- data.frame( 
  sage_30m = 1, 
  sl_ = 5000,
  ta_ = 0 )

# now we use log_rss() to calculate log-RSS 
lr1 <- log_rss( m_sg, x1 = s1, x2 = s2 )
#  view
lr1$df

# we can shift to using a range for one of the distributions instead
# so that we  can make a plot of change:
s3 <- data.frame( 
  sage_30m = seq(from = -1, to = 3, length.out = 100), 
  sl_ = 5000,
  ta_ = 0 )

# Calculate log-RSS
lr2 <- log_rss(m_sg, s3, s1)

# Plot log-RSS using ggplot2
ggplot( lr2$df, aes(x = sage_30m_x1, y = log_rss) ) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  xlab("Sagebrush (SD)") +
  ylab("log-RSS vs Mean Sagebrush") +
  theme_bw( base_size = 16 )

# Plot RSS using ggplot2
ggplot( lr2$df, aes(x = sage_30m_x1, y = exp(log_rss) ) ) +
  geom_line(size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray30") +
  xlab("Sagebrush (SD)") +
  ylab("RSS vs Mean Sagebrush") +
  theme_bw( base_size = 16 )


###

#Assign the empirical distributions to our model:
m_sg$sl_ <- emp_d_sl
m_sg$ta_ <- emp_d_ta

# Now we can use coefficients associated with movement parameters
# to update our movement related distributions for the same individual. #
# the amt has an update function, which will work if you do not have #
# any interactions in your models. Otherwise, refer to appendix c in #
# Fieberg et al. 2021 to see the equations that you need to use to manually #
# update those parameters #
# Update step length distribution:
updated_sl <- update_sl_distr( m_sg, 
            beta_sl = 'sl_', beta_log_sl =  'log(sl_)' )
#update turning angle distribution:
updated_ta <- update_ta_distr( m_sg, 
              beta_cos_ta = "cos(ta_)" )
#view
updated_sl
updated_ta

#plot movement parameter distributions
# Start by reminding ourselves of step lengths of the individual
df %>% dplyr::filter( id == 4 ) %>%
  ggplot( . ) +
  theme_bw( base_size = 16 ) +
  geom_histogram( aes( x = sl_ ) )

#plot empirical gamma for step lengths
par( mfrow = c(2,1))
plot( density( dgamma( 1:300, shape = emp_d_sl$params$shape,
                       scale = emp_d_sl$params$scale ) ) , 
#      xlim = c( 0, 0.001)
      )
#plot updated gamma distribution:
plot( density( dgamma( 1:300, shape = updated_sl$params$shape,
                       scale = updated_sl$params$scale ) ),
              lwd = 3, 
#     xlim = c( 0, 0.001) 
     )

#histogram of turning angles for the individual:
df %>% dplyr::filter( id == 4 ) %>%
  ggplot( . ) +
  theme_bw( base_size = 16 ) +
  geom_histogram( aes( x = ta_ ) )

#Von mises empirical distribution of turning angles:
plot( density( dvonmises( seq(from = -1 * pi, to = pi, length.out = 100), 
                          mu = emp_d_ta$params$mu,
                          kappa = emp_d_ta$params$kappa ) ) )
#updated distribution:
plot( density( dvonmises( seq(from = -1 * pi, to = pi, length.out = 100), 
                          mu = updated_ta$params$mu,
                          kappa = updated_ta$params$kappa ) ) )
# Note that the kappa parameter in the updated distribution is negative,#
# which is not allowed. Maybe this is due to not removing the non-movement #
# locations ahead of creating steps? # 

# We can also use hypothetical locations to interpret how the individual #
# is moving. we use the updated distribution to estimate the likelihood #
# under a selection-free step-length distribution of taking a 5,000-m-step #
# which for our individual at 30 minute resolution is short, and a #
# long step at double that. # 

# estimate likelihood for short step:
short <- dgamma(50, 
                 shape = updated_sl$params$shape,
                 scale = updated_sl$params$scale)
#now for long step
long <- dgamma(200, 
                  shape = updated_sl$params$shape,
                  scale = updated_sl$params$scale)
# calculate selection:
short/long
# individual is 47 times more likely to take the shorter than the 
# longer step when all habitat conditions are the same
####end analysis of single individual ###
#######  fit movement model for all individuals in glmmTMB ############
#Start by replicating the approach for one individual only #
# the only difference is that we add weights:
#create subsetted dataframe
newdf <- df %>%  filter( id == 4 )
#define model structure 
m1.struc <- glmmTMB( case_ ~ sage_30m + 
                    #add movement parameters
                    sl_ + log(sl_) + 
                      cos(ta_) +
                    #define random effects
                    ( 1| step_id_ ),# + 
                    family = poisson, data = newdf, 
                    weights = weight,
                    doFit=FALSE ) 

# fix variance
m1.struc$parameters$theta[ 1 ] <- log( 1e3 ) 
# tell it not to change variance for strata
m1.struc$mapArg <- list( theta = factor( NA ) )

#then fit the model
m1 <- glmmTMB::fitTMB( m1.struc )
#view results
summary( m1 )
# results are similar

# After testing our model on a single individual we move to estimate #
# parameters for all individuals:
# for all individuals we can add the individual random intercepts:
m2.struc <- glmmTMB( case_ ~ sage_30m + 
                       #add movement parameters
                       sl_ + log(sl_) + cos(ta_) +
                       #define random effects
                       ( 1| step_id_ ) + 
                       ( 1| id ), 
                     family = poisson, data = df, 
                     weights = weight, doFit=FALSE ) 

# fix variance
m2.struc$parameters$theta[ 1 ] <- log( 1e3 ) 
# tell it not to change variance
m2.struc$mapArg <- list( theta = factor( c(NA, 1) ) )

#then fit the model
m2 <- glmmTMB::fitTMB( m2.struc )
summary( m2 )

# We can also add random slopes for movement and habitat parameters
m3.struc <- glmmTMB( case_ ~ sage_30m + 
                       #add movement parameters
                       sl_ + log(sl_) + cos(ta_) +
                       #define random effects
                       ( 1| step_id_ ) + 
                       ( 1| id ) +
                       ( 0 + sage_30m | id ) +
                       ( 0 + sl_ | id ) +
                       ( 0 + log(sl_) | id ) +
                       ( 0 + cos(ta_) | id ),
                     family = poisson, data = df, 
                     weights = weight, doFit=FALSE ) 

# fix variance
m3.struc$parameters$theta[ 1 ] <- log( 1e3 ) 
# tell it not to change variance
m3.struc$mapArg <- list( theta = factor( c(NA, 1:5) ) )

#then fit the model
m3 <- glmmTMB::fitTMB( m3.struc )
summary( m3 )

# it seems like our model was overly ambitious and couldn't estimate
# slopes for the step length parameters. We simplify it assuming #
# that they have similar slopes among individuals
m4.struc <- glmmTMB( case_ ~ sage_30m + 
                       #add movement parameters
                       sl_ + log(sl_) + cos(ta_) +
                       #define random effects
                       ( 1| step_id_ ) + 
                       ( 1| id ) +
                       ( 0 + sage_30m | id ) +
                       ( 0 + cos(ta_) | id ),
                     family = poisson, data = df, 
                     weights = weight, doFit=FALSE ) 

# fix variance
m4.struc$parameters$theta[ 1 ] <- log( 1e3 ) 
# tell it not to change variance
m4.struc$mapArg <- list( theta = factor( c(NA, 1:3) ) )

#then fit the model
m4 <- glmmTMB::fitTMB( m4.struc )
summary( m4 )
# this model has a lower AIC than the one that doesn't include #
# any random slopes and so it seems like we can account for individual #
# specific differences between selection for sagebrush and movement #
# in relation to turning angles #


##########################################################################
### Save desired results                                  #

#save workspace if in progress
save.image( 'iSSF_results.RData' )
############# end of script  ##################################