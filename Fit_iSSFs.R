##################################################################
# Script developed by Jen Cruz to estimate iSSFs          #
# approach derived from Fieberg et al. 2021 and Signer et al. 2019 #
# using code from Appendices B and C for the single individual    #
# example and for all individuals using code from  Muff et al. 2019: #
# https://conservancy.umn.edu/handle/11299/204737                #
#                                                                #
# Vegetation cover  was downloaded from Rangeland Analysis Platform #
# https://rangelands.app/products/ for 2021 and includes        #
# % cover for shrub, perennial herbaceous, annual herbaceous    #
# tree, litter and bare ground                                   #
# coordinate system is WGS84 EPSG:4326, spatial resolution is 30m #
# and was extracted at two scales                                #
# Prairie Falcon data was thinned to 30minutes for 9 individuals #
# tracked in 2021.                                               #
###################################################################

################## prep workspace ###############################
#install relevant packages
install.packages( "circular" )

# load packages relevant to this script:
library( tidyverse ) #easy data manipulation
# set option to see all columns and more than 10 rows
options( dplyr.width = Inf, dplyr.print_min = 100 )
library( amt )
library( glmmTMB )
library( circular ) #for plotting von mises distribution
#####################################################################
## end of package load ###############

###################################################################
#### Load or create data -----------------------------------------
# Clean your workspace to reset your R environment. #
rm( list = ls() )
#load 30m steps estimated for all individuals and habitat 
# variables extracted for each step
df_steps <- read_rds( "Data/df_steps" )

#load the scaled data so we don't have to do it again
df_scl <- read_rds( "Data/df_scl" )

#######################################################################
######## visualizing data ###############################################

#check
head( df_scl )
# We cannot have missing values for the predictors so we 
# need to remove steps with missing ta_ 
#which steps have missing ta values:
df_steps$step_id_[ is.na( df_steps$ta_ )  ]
# do we have any zero step lengths:
df_steps$step_id_[ which( df_steps$sl_ == 0) ] 

#recheck sample size:
table( df_steps$territory )

#To check whether there is evidence of individual differences 
#in habitat selection we could categorize the vegetation metrics for
#plotting purposes. We create categories in a new plotting dataframe:
df_plot <- df_steps %>% 
  dplyr::mutate( peren_cat = ifelse( perennial_30m < 20, "low",
  ifelse( perennial_30m >  40, "high", "medium" ) ),
  annual_cat = ifelse( annual_30m < 20, "low",
    ifelse( annual_30m > 40, "high", "medium" ) ),
  shrub_cat = ifelse( shrub_30m < 10, "low",
    ifelse( shrub_30m > 20, "high", "medium" ) ))

# Start by plotting differences in step lengths between low,medium
# and high percentages of perennial herbaceous vegetation.
#before plotting we remove missing values 
df_plot %>% filter( !is.na(peren_cat) ) %>% 
ggplot( . ) +
  theme_bw( base_size = 15 ) +
  geom_density( aes( x = sl_, 
                     fill = case_, group = case_ ),
                alpha = 0.5  )  +
  xlim( 0,20000 ) +
  facet_wrap( ~ peren_cat, ncol = 1 )

#repeat for annual
df_plot %>% filter( !is.na(annual_cat) ) %>% 
  ggplot( . ) +
  theme_bw( base_size = 15 ) +
  geom_density( aes( x = sl_, 
                     fill = case_, group = case_ ),
                alpha = 0.5  )  +
  xlim( 0,20000 ) +
  facet_wrap( ~ annual_cat, ncol = 1 )

#repeat for shrub
df_plot %>% filter( !is.na(shrub_cat) ) %>% 
  ggplot( . ) +
  theme_bw( base_size = 15 ) +
  geom_density( aes( x = sl_, 
                     fill = case_, group = case_ ),
                alpha = 0.5  )  +
  xlim( 0,20000 ) +
  facet_wrap( ~ shrub_cat, ncol = 1 )

# How do you interpret these results?
# Answer: 
# 

# Now we look at turning angles 
df_plot %>% filter( !is.na(peren_cat) ) %>% 
  ggplot( . ) +
  theme_bw( base_size = 15 ) +
  geom_density( aes( x = ta_, 
                     fill = case_, group = case_ ),
                alpha = 0.5  )  +
  facet_wrap( ~ peren_cat, ncol = 1 )
#repeat for annual
df_plot %>% filter( !is.na(annual_cat) ) %>% 
  ggplot( . ) +
  theme_bw( base_size = 15 ) +
  geom_density( aes( x = ta_, 
                     fill = case_, group = case_ ),
                alpha = 0.5  )  +
  facet_wrap( ~ annual_cat, ncol = 1 )

#repeat for shrub
df_plot %>% filter( !is.na(shrub_cat) ) %>% 
  ggplot( . ) +
  theme_bw( base_size = 15 ) +
  geom_density( aes( x = ta_, 
                     fill = case_, group = case_ ),
                alpha = 0.5  )  +
  facet_wrap( ~ shrub_cat, ncol = 1 )

# How do you interpret these results?
# Answer: 
# 
# For homework try visualizing plots with 100m scale and 
# use results to decide which scale to analyse
# 
### end visualizing prep #######
######## analysis for single individual using amt #################
# Remember that amt does not allow random effects...but we still 
# want to explore features that it does offer. Another alternative #
# to using random effects would be to run analysis separate for each #
# individual so we demonstrate by running analysis for one individual #
# subset dataframe to one
df_one <- df_scl %>% dplyr::filter( id == 3 )

#Fit iSSF model to same individual:
m1_sg <- df_one %>% 
  #remember that we remove the intercept
  fit_issf( case_ ~ -1 + annual_30m + perennial_30m + shrub_30m +
              sl_ + log(sl_) + cos(ta_) +
              #here we specify that it is a conditional model:
              strata( step_id_ ), model = TRUE )
# Note the cosine of the turn angle and the log of the 
#step-length allow us to adjust/refine the parameters of our #
# tentative step-length and turn-angle distributions after #
# fitting our integrated step-selection model. #

#view results
summary( m1_sg )

# There is no evidence of the need to modify the gamma distribution
# of the step length, or the von misses for the turning angle.
# How do we know this?
# Answer:
# which habitat predictors seem to matter to this falcon?
# Answer:
# 

# From our visualizations earlier we did expect individuals 
# to use different movements in different habitats. Here we test 
# out this theory by adding interaction terms:
m2_sg <- df_one %>% 
  #remember that we remove the intercept
  fit_issf( case_ ~ -1 + annual_30m + perennial_30m + shrub_30m +
              sl_ + log(sl_) + cos(ta_) +
              #add interactions between movements and habitats
              log(sl_):annual_30m + log(sl_):perennial_30m + 
              log(sl_):shrub_30m +
              cos(ta_):annual_30m + cos(ta_):perennial_30m +
              cos(ta_):shrub_30m + 
              #here we specify that it is a conditional model:
              strata( step_id_ ), model = TRUE )

#view results
summary( m2_sg )
# What do these results tell you? Did the model converge?
# Answer:
# 
# we try reduced models
m3_sg <- df_one %>% 
  #remember that we remove the intercept
  fit_issf( case_ ~ -1 + annual_30m + perennial_30m + shrub_30m +
              sl_ + log(sl_) + cos(ta_) +
              #add interactions between movements and habitats
              log(sl_):annual_30m + #log(sl_):perennial_30m + 
              #log(sl_):shrub_30m +
              cos(ta_):annual_30m + #cos(ta_):perennial_30m +
              #cos(ta_):shrub_30m + 
              #here we specify that it is a conditional model:
              strata( step_id_ ), model = TRUE )

#view results
summary( m3_sg )

m4_sg <- df_one %>% 
  #remember that we remove the intercept
  fit_issf( case_ ~ -1 + annual_30m + perennial_30m + shrub_30m +
              sl_ + log(sl_) + cos(ta_) +
              #add interactions between movements and habitats
              log(sl_):perennial_30m +
              cos(ta_):perennial_30m +
              #here we specify that it is a conditional model:
              strata( step_id_ ), model = TRUE )

#view results
summary( m4_sg )

m5_sg <- df_one %>% 
  #remember that we remove the intercept
  fit_issf( case_ ~ -1 + annual_30m + perennial_30m + shrub_30m +
              sl_ + log(sl_) + cos(ta_) +
              #add interactions between movements and habitats
              log(sl_):shrub_30m + cos(ta_):shrub_30m + 
              #here we specify that it is a conditional model:
              strata( step_id_ ), model = TRUE )

#view results
summary( m5_sg )

#we can compare model results
AIC( m1_sg); AIC( m2_sg); AIC( m3_sg); AIC( m4_sg); AIC( m5_sg)
# Interpret these AIC results. Which model has the most support?
# Can we be confident in all models?
# Answers:
#
#

# We end by adding a quadratic term to the top model
m6_sg <- df_one %>% 
  fit_issf( case_ ~ -1 + annual_30m + perennial_30m + shrub_30m +
              I( annual_30m^2 ) +
              sl_ + log(sl_) + cos(ta_) +
              strata( step_id_ ), model = TRUE )
summary( m6_sg)
AIC( m1_sg); AIC( m2_sg)
# Was it supported?
# Answer:
# 
#### Visualize results for top model ###
# One option that can help interpret results is by calculating #
# relative selection strength for two locations that are equal #
# except for 1 of the predictors. We demonstrate that approach here: #

# start by creating new dataframes that contain details for each #
# of the two hypothetical locations that we want to compare: 
# We focus on differences in sagebrush #
s1 <- data.frame( 
  annual_30m = 0, 
  perennial_30m = 0,
  shrub_30m = 0,
  sl_ = 1000,
  ta_ = 0 )
# now a second dataframe with higher sagebrush
s2 <- data.frame( 
  annual_30m = 1, 
  perennial_30m = 0,
  shrub_30m = 0,
  sl_ = 1000,
  ta_ = 0 )

# now we use log_rss() and our top model to calculate log-RSS 
lr1 <- log_rss( m1_sg, x1 = s1, x2 = s2 )
#  view
lr1$df

# we can shift to using a range for one of the distributions instead
# so that we  can make a plot of change:
#first we need to know what the individual experience
min( df_one$annual_30m ); max(df_one$annual_30m )
# use those values for prediction
s3 <- data.frame( 
  annual_30m = seq(from = round( min( df_one$annual_30m ),0), 
            to = round(max( df_one$annual_30m ),0), length.out = 100), 
  perennial_30m = 0,
  shrub_30m = 0,
  sl_ = 1000,
  ta_ = 0 )

# Calculate log-RSS for our top model
lr2 <- log_rss(m2_sg, s3, s1)

#view
head( lr2$df )

# Plot log-RSS using ggplot2
ggplot( lr2$df, aes(x = annual_30m_x1, y = log_rss) ) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  xlab("Annual (SD)") +
  ylab("log-RSS vs Mean Annual") +
  theme_bw( base_size = 16 )
#the line depicts the log-RSS between each value of s1 relative to s2.
# Remember that these locations only differ in their values of annual

# To plot the changes in relative selection strength instead #
# we exponentiate the log_rss
ggplot( lr2$df, aes(x = annual_30m_x1, y = exp(log_rss) ) ) +
  geom_line(size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray30") +
  xlab("Annual (SD)") +
  ylab("RSS vs Mean Annual") +
  theme_bw( base_size = 16 )


########## compare our empirical movement parameters to those estimated #
# from the top model:
# Estimate empirical distributions for individual of interest:
emp_d_sl <- df_one %>%
  dplyr::select( sl_ ) %>% 
  fit_distr( ., dist_name = "gamma" )
#view
emp_d_sl
#now for turning angle using a von mises, circular distribution:
emp_d_ta <- df_one %>% 
  dplyr::select( ta_ ) %>% 
  fit_distr( ., dist_name = "vonmises" )
#view
emp_d_ta

#Assign the empirical distributions to our model:
m1_sg$sl_ <- emp_d_sl
m1_sg$ta_ <- emp_d_ta

# Now we can use coefficients associated with movement parameters
# to update our movement related distributions for the same individual. #
# the amt has an update function, which will work if you do not have #
# any interactions in your models. Otherwise, refer to appendix c in #
# Fieberg et al. 2021 to see the equations that you need to use to manually #
# update those parameters #
# Update step length distribution:
updated_sl <- update_sl_distr( m1_sg, 
            beta_sl = 'sl_', beta_log_sl =  'log(sl_)' )
#update turning angle distribution:
updated_ta <- update_ta_distr( m1_sg, 
              beta_cos_ta = "cos(ta_)" )
#view
updated_sl
updated_ta

#plot movement parameter distributions
# Start by reminding ourselves of step lengths of the individual
  ggplot( df_one ) +
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
ggplot( df_one ) +
  theme_bw( base_size = 16 ) +
  geom_histogram( aes( x = ta_ ) )

#Von mises empirical distribution of turning angles:
plot( density( circular::dvonmises( seq(from = -1 * pi, to = pi, length.out = 100), 
                          mu = emp_d_ta$params$mu,
                          kappa = emp_d_ta$params$kappa ) ) )
#updated distribution:
plot( density( circular::dvonmises( seq(from = -1 * pi, to = pi, length.out = 100), 
                          mu = updated_ta$params$mu,
                          kappa = updated_ta$params$kappa ) ) )

# We can also use hypothetical locations to interpret how the individual #
# is moving. we use the updated distribution to estimate the likelihood #
# under a selection-free step-length distribution of taking a step #
# which for our individual is short, and a long step. # 

# estimate likelihood for short step:
short <- dgamma( 100, 
                 shape = updated_sl$params$shape,
                 scale = updated_sl$params$scale )
#now for long step
long <- dgamma( 5000, 
                  shape = updated_sl$params$shape,
                  scale = updated_sl$params$scale )
# calculate selection:
short/long
# individual is 11.9 times more likely to take the shorter than the 
# longer step when all habitat conditions are the same
####end analysis of single individual ###
#########################################################################
#######  fit movement model for all individuals in glmmTMB ############
# for all individuals firstly ignoring individual differences:
m1 <- glmmTMB( case_ ~ 0 + annual_30m + perennial_30m + shrub_30m + 
                       #add movement parameters
                       sl_ + log(sl_) + cos(ta_) +
                       #define random effects
                       ( 1| step_id_ ), 
                     family = poisson, data = df_scl, 
                     #define weights
                     weights = weight, 
                     #tell it not to change variance for step level
                     map = list( theta = factor( c(NA ) ) ),
                     #fix variance for step level random intercept
                     start = list( theta = c( log( 1e3 ) ) )
                     ) 
#view
summary( m1 )
# what happened?
# Answer:
# 
# What are possible solutions? 
# Answer
df_scl$sl_ = df_scl$sl_ /1000 

#Next by adding random slopes for movement and habitat parameters
m2 <- glmmTMB( case_ ~ 0 + annual_30m + perennial_30m + shrub_30m +  
                       #add movement parameters
                       sl_ + log(sl_) + cos(ta_) +
                       #define random effects
                       ( 1| step_id_ ) + 
                       ( 0 + annual_30m | id ) +
                       ( 0 + perennial_30m | id ) +
                       ( 0 + shrub_30m | id ) +
                      ( 0 + sl_ | id ) +
                       ( 0 + log(sl_) | id ) +
                       ( 0 + cos(ta_) | id ),
              family = poisson, data = df_scl, 
               #define weights
               weights = weight, 
               #tell it not to change variance for step level
               map = list( theta = factor( c(NA, 1:6 ) ) ),
               #fix variance for step level random intercept
               start = list( theta = c( log( 1e3 ),0,0,0,0,0,0 ) )
               ) 
#view
summary( m2 )

# We have another go at interactions
m3 <- glmmTMB( case_ ~ 0 + annual_30m + perennial_30m + shrub_30m +  
                       #add movement parameters
                  sl_ +log(sl_) + cos(ta_) +
                 #add interactions between movements and habitats
                 log(sl_):annual_30m + log(sl_):perennial_30m + 
                 log(sl_):shrub_30m +
                 cos(ta_):annual_30m + cos(ta_):perennial_30m +
                 cos(ta_):shrub_30m +
                 #define random effects
                 ( 1| step_id_ ) + 
                 ( 0 + annual_30m | id ) +
                 ( 0 + perennial_30m | id ) +
                 ( 0 + shrub_30m | id ) +
                 ( 0 + sl_ | id ) +
                 ( 0 + log(sl_) | id ) +
                 ( 0 + cos(ta_) | id ),
               family = poisson, data = df_scl, 
               #define weights
               weights = weight, 
               #tell it not to change variance for step level
               map = list( theta = factor( c(NA, 1:6 ) ) ),
               #fix variance for step level random intercept
               start = list( theta = c( log( 1e3 ),0,0,0,0,0,0 ) )
)


summary( m3 )

# Which model is our top model? 
# How do you interpret results from it?
# Answer:
# 

##############################################################################
################### visualizing top model results ############
#pull out random effects at the id level #
ran.efs <- ranef( m2 )$cond$id
#note that we don't want the ones at the step level

#pull out fixed effects
fix.efs <- fixef( m2 )$cond
#view
fix.efs

#we need to add the fixed effect to the random for each vegetation 
# and exponentiate our results
rss <- ran.efs
rss[,1 ] <- exp( rss[,1] + fix.efs[1] )
rss[,2 ] <- exp( rss[,2] + fix.efs[2] )
rss[,3 ] <- exp( rss[,3] + fix.efs[3] )
rss[,4 ] <- exp( rss[,4] + fix.efs[4] )
rss[,5 ] <- exp( rss[,5] + fix.efs[5] )
rss[,6 ] <- exp( rss[,5] + fix.efs[6] )
#create id column
rss$id <- as.numeric(  rownames( rss ) )
#view
round(rss,2)
# now extract additional details from our steps dataframe to combine 
# with our results
iddf <- df_steps %>% 
  group_by( id, territory, sex ) %>% 
  summarise( annual_30m_mean = mean( annual_30m, na.rm = TRUE),
             perennial_30m_mean = mean( perennial_30m, na.rm = TRUE),
             shrub_30m_mean = mean( shrub_30m, na.rm = TRUE),
             mean_sl =  mean(sl_, na.rm = TRUE),
             log_sl = log( mean(sl_, na.rm = TRUE) ),
             cos_ta = cos( mean(ta_, na.rm = TRUE) ) 
  )
iddf
#combine with our resource selection strength estimates
iddf <- left_join( iddf, rss, by = "id" )

#plot results
ggplot( iddf ) +
  theme_classic( base_size = 15 ) +
  labs( x = "Mean annual cover (%)", 
        y = "Resource selection strength" ) +
  geom_point( aes( x = annual_30m_mean , y = annual_30m, color = sex ) ) +
  geom_hline( yintercept = 1, linewidth = 1 )

ggplot( iddf ) +
  theme_classic( base_size = 15 ) +
  labs( x = "Mean perennial cover (%)", 
        y = "Resource selection strength" ) +
  geom_point( aes( x = perennial_30m_mean , y = perennial_30m, color = sex ) ) +
  geom_hline( yintercept = 1, linewidth = 1 )

ggplot( iddf ) +
  theme_classic( base_size = 15 ) +
  labs( x = "Mean shrub cover (%)", 
        y = "Resource selection strength" ) +
  geom_point( aes( x = shrub_30m_mean , y = shrub_30m, color = sex ) ) +
  geom_hline( yintercept = 1, linewidth = 1 )

ggplot( iddf ) +
  theme_classic( base_size = 15 ) +
  labs( x = "Mean perennial cover (%)", 
        y = "Resource selection strength" ) +
  geom_point( aes( x = perennial_30m_mean , y = perennial_30m, color = sex ) ) +
  geom_hline( yintercept = 1, linewidth = 1 )


# What else could we add to our analysis?
# Answer: 
# 
##########################################################################
### Save desired results                                  #

#save workspace if in progress
save.image( 'iSSF_results.RData' )
############# end of script  ##################################