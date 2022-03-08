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

install.packages( "circular" )
# load packages relevant to this script:
library( tidyverse ) #easy data manipulation
# set option to see all columns and more than 10 rows
options( dplyr.width = Inf, dplyr.print_min = 100 )
library( amt )
library( sf )
library( glmmTMB )
library( circular ) #for plotting von mises distribution

#####################################################################
## end of package load ###############

###################################################################
#### Load or create data -----------------------------------------

#load 30m steps estimated for all individuals
trks <- read_rds( "trks_all" )
#check
class(trks)
#now import dataframe with sagebrush values
dfraw <- read_rds( "df_all" )
#check
class( dfraw )
#import polygon of the NCA as sf spatial file:
NCA_Shape <- sf::st_read("Z:/Common/QCLData/Habitat/NCA/GIS_NCA_IDARNGpgsSampling/BOPNCA_Boundary.shp")


#######################################################################
######## preparing data ###############################################
dim(dfraw)
#check dataframe
head( df );dim(df)
# now tks
head( trks ); dim(trks)

# Remember in our prior script we created the random steps, which were #
# fit using step-length and turning angle distributions. #
df <- dfraw
# we check sample sizes for each individual
table( df$id)
#select ids those with poor sample size
remids <- c( 1,3,7)
# use those to remove their data from dataframe:
df <- df %>%  
  dplyr::filter( !(id %in% remids) )
#check 
dim(df)

# We cannot have missing values for the predictors so we 
# need to remove steps with missing ta_ 
#step ids are not unique to individuals so we create a unique id:
df$id_step <- paste0( df$id, df$step_id_ )
head( df);dim(df)

#which steps have missing ta values:
df$id_step[ which( is.na(df$ta_) ) ]
#how many don't:
length( df$id_step[ which( !is.na(df$ta_) ) ] )

# which steps have missing sl values:
df$id_step[ which(is.na(df$sl_)) ] 
#none so we focus on ta only

# we record ids for those with missing ta values:
rem <- df$id_step[ which(is.na(df$ta_)) ] 

#remove them from dataframe
df <- df %>%  
  dplyr::filter( !(id_step %in%  rem ) )

dim(df)
#recheck sample size:
table( df$id)

### end data prep ####
######### analysis for single individual using amt #################
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
plot( density( dgamma( 1:30000, shape = emp_d_sl$params$shape,
                       scale = emp_d_sl$params$scale ) ) , 
#      xlim = c( 0, 0.001)
      )
#plot updated gamma distribution:
plot( density( dgamma( 1:30000, shape = updated_sl$params$shape,
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
short <- dgamma(5000, 
                 shape = updated_sl$params$shape,
                 scale = updated_sl$params$scale)
#now for long step
long <- dgamma(10000, 
                  shape = updated_sl$params$shape,
                  scale = updated_sl$params$scale)
# calculate selection:
short/long
# individual is 3.5 times more likely to take the shorter than the 
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
# results are similar, except that sl_ parameter is not 'significant'

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