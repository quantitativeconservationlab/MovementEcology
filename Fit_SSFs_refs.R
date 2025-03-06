##################################################################
# Script developed by Jen Cruz to estimate SSFs and iSSFs          #
# using code from Muff et al. 2019 DOI: 10.1111/1365-2656.13087  #
# code here:                                                     #  
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
##################################################################
# Script developed by Jen Cruz to estimate SSFs and iSSFs          #
# approach derived from Fieberg et al. 2021 and Signer et al. 2019 #
# using code from Appendices B and C                             #
# also vignette here:
# https://conservancy.umn.edu/server/api/core/bitstreams/63727072-87b1-4b35-b81c-8fd31b8f1e57/content #
# Vegetation cover  was downloaded from Rangeland Analysis Platform #
# https://rangelands.app/products/ for 2021 and includes        #
# % cover for shrub, perennial herbaceous, annual herbaceous    #
# tree, litter and bare ground                                   #
# coordinate system is WGS84 EPSG:4326, spatial resolution is 30m #
#                                                                #
# Prairie Falcon data was thinned to 30minutes for 9 individuals #
# tracked in 2021 and uses NAD83 UTM zone 11N +                   #
# which is the same as the NCA polygon                           #
###################################################################

################## prep workspace ###############################
# load packages relevant to this script:
library( tidyverse ) #easy data manipulation
# set option to see all columns and more than 10 rows
options( dplyr.width = Inf, dplyr.print_min = 100 )
library( amt )
library( glmmTMB ) # for analysis
library( circular ) #for plotting von mises distribution

#####################################################################
## end of package load ###############

###################################################################
#### Load or create data -----------------------------------------

# Clean your workspace to reset your R environment. #
rm( list = ls() )
#load 30m steps estimated for all individuals 
df_steps <- read.csv( "Data/df_steps30.csv" )
# load 30m steps with scaled predictors
df_scl <- read.csv( "Data/df_scl.csv" )

#######################################################################
######## preparing data ###############################################
#view data
head( df_steps)
#recheck sample size:
table( df_steps$id )

#To check whether there is evidence of individual differences 
#in habitat selection we could categorize the vegetation metrics for
#plotting purposes. We create categories in a new plotting dataframe:
df_plot <- df_steps %>% 
  dplyr::mutate( peren_cat = ifelse( perennial < 20, "low", "high" ),
   annual_cat = ifelse( annual < 30, "low","high" ),
    shrub_cat = ifelse( shrub < 10, "low","high" ) )

# Start by plotting differences in step lengths between low,medium
# and high percentages of perennial herbaceous vegetation.
#before plotting we remove missing values 
df_plot %>% filter( !is.na(peren_cat) ) %>% 
  ggplot( . ) +
  theme_bw( base_size = 10 ) +
  geom_density( aes( x = sl_, 
                     fill = case_, group = case_ ),
                alpha = 0.5  )  +
  xlim( 0,20000 ) +
  facet_wrap( id ~ peren_cat, scales = "free" )

#repeat for annual
df_plot %>% filter( !is.na(annual_cat) ) %>% 
  ggplot( . ) +
  theme_bw( base_size = 10 ) +
  geom_density( aes( x = sl_, 
          fill = case_, group = case_ ),
                alpha = 0.5  )  +
  xlim( 0,20000 ) +
  facet_wrap( id ~ annual_cat, scales = "free" )

#repeat for shrub
df_plot %>% filter( !is.na(shrub_cat) ) %>% 
  ggplot( . ) +
  theme_bw( base_size = 10 ) +
  geom_density( aes( x = sl_, 
                     fill = case_, group = case_ ),
                alpha = 0.5  )  +
  xlim( 0,20000 ) +
  facet_wrap( id ~ shrub_cat, scales = "free" )

# How do you interpret these results?
# Answer: 
# 

# Now we look at turning angles 
df_plot %>% filter( !is.na(peren_cat) ) %>% 
  ggplot( . ) +
  theme_bw( base_size = 10 ) +
  geom_density( aes( x = ta_, 
                     fill = case_, group = case_ ),
                alpha = 0.5  )  +
  facet_wrap( id ~ peren_cat, scales = "free" )
#repeat for annual
df_plot %>% filter( !is.na(annual_cat) ) %>% 
  ggplot( . ) +
  theme_bw( base_size = 10 ) +
  geom_density( aes( x = ta_, 
                     fill = case_, group = case_ ),
                alpha = 0.5  )  +
  facet_wrap( id ~ annual_cat, scales = "free" )

#repeat for shrub
df_plot %>% filter( !is.na(shrub_cat) ) %>% 
  ggplot( . ) +
  theme_bw( base_size = 15 ) +
  geom_density( aes( x = ta_, 
                     fill = case_, group = case_ ),
                alpha = 0.5  )  +
  facet_wrap( id ~ shrub_cat, scales = "free" )

# How do you interpret these results?
# Answer: 
# 
### end visualizing prep #######
#########################################################################
#######  fit movement model for all individuals in glmmTMB ############
# for all individuals firstly ignoring individual differences:
m1 <- glmmTMB( case_ ~ 0 + 
                 #add habitat predictors
                 annual + perennial + shrub + 
                 #add movement parameters
                 sl_ + log_sl_ + cos_ta_ +
                 # add movement interactions with habitat
                 log_sl_*annual + cos_ta_*annual +
                 log_sl_*perennial + cos_ta_*perennial +
                 log_sl_*shrub + cos_ta_*shrub +
                 #add random intercept for step id (stratum)
                 #to ensure pairing of random steps to their used step
                 ( 1| step_id_ ) +
                 #define random slopes for habitat
                 ( 0 + annual | id ) +
                 ( 0 + perennial | id ) +
                 ( 0 + shrub | id ),
               family = poisson, data = df_scl, 
               #define weights
               weights = weight, 
               #tell it not to change variance for step level
               map = list( theta = factor( c(NA, 1:3 ) ) ),
               #fix variance for step level random intercept
               start = list( theta = c( log( 1e3 ),0,0,0 ) )
) 
#view
summary( m1 )

### We compare against a simpler model that assumes no interactions
m2 <- glmmTMB( case_ ~ 0 + annual + perennial + shrub + 
                 #add movement parameters
                 sl_ + log_sl_ + cos_ta_ +
                 #define random effects
                 ( 1| step_id_ ) +
                 # add random slopes for habitat variables
                 ( 0 + annual | id ) +
                 ( 0 + perennial | id ) +
                 ( 0 + shrub | id ),
               family = poisson, data = df_scl, 
               #define weights
               weights = weight, 
               #tell it not to change variance for step level
               map = list( theta = factor( c(NA, 1:3 ) ) ),
               #fix variance for step level random intercept
               start = list( theta = c( log( 1e3 ),0,0,0 ) )
) 
#view
summary( m2 )

# compare models:
anova(m1,m2)
#Which model had the most support?
# Answer:
# 
# How do you interpret results of model comparison?
# Answer:
#

################### visualizing top model results ############
# select top model
mr <- m2
#pull out random effects at the id level #
ran.efs <- ranef( mr )$cond$id
#note that we don't want the ones at the step level

#pull out fixed effects
fix.efs <- fixef( mr )$cond
#view
fix.efs

#we need to add the fixed effect to the random for each vegetation 
# and exponentiate our results
#make sure that the random and fixed effect order match 
rss <- ran.efs
rss[,1 ] <- exp( rss[,1] + fix.efs[1] )
rss[,2 ] <- exp( rss[,2] + fix.efs[2] )
rss[,3 ] <- exp( rss[,3] + fix.efs[3] )
#create id column
rss$id <- as.numeric(  rownames( rss ) )
#view
round(rss,2)
# now extract additional details from our steps dataframe to combine 
# with our results
iddf <- df_steps %>% 
  group_by( id, territory, sex ) %>% 
  summarise( annual_mean = mean( annual, na.rm = TRUE),
             perennial_mean = mean( perennial, na.rm = TRUE),
             shrub_mean = mean( shrub, na.rm = TRUE),
    #remember that we converted sl to km in df_scl last week
             mean_sl =  mean(sl_/1000, na.rm = TRUE) ,
             log_sl = log( mean(sl_/1000, na.rm = TRUE) ),
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
  geom_point( aes( x = annual_mean , y = annual, color = sex ) ) +
  #we add population-level average
  geom_hline( yintercept = exp(fix.efs[1]), linewidth = 1 ) +
  geom_hline( yintercept = 1, lty = 2 )

ggplot( iddf ) +
  theme_classic( base_size = 15 ) +
  labs( x = "Mean perennial cover (%)", 
        y = "Resource selection strength" ) +
  geom_point( aes( x = perennial_mean , y = perennial, color = sex ) ) +
  geom_hline( yintercept = exp(fix.efs[2]), linewidth = 1 ) +
  geom_hline( yintercept = 1, lty = 2 )

ggplot( iddf ) +
  theme_classic( base_size = 15 ) +
  labs( x = "Mean shrub cover (%)", 
        y = "Resource selection strength" ) +
  geom_point( aes( x = shrub_mean , y = shrub, color = sex ) ) +
  geom_hline( yintercept = exp(fix.efs[3] ), linewidth = 1 ) +
  geom_hline( yintercept = 1, lty = 2 )

# Could we have missinterpreted habitat selection if we had ignored
# Random slopes?
# How?
# Answer:
#
########### visualize the movement distributions #######
# We calculate the tentative distributions from empirical data 
# Start with step length fitted as a gamma with shape and scale parameters
emp_d_sl <- df_scl %>% 
  dplyr::select( sl_ ) %>% 
  #fit a gamma distribution using empirical data
  amt::fit_distr(., dist_name = "gamma" )

#Fit a von misses to the turning angles for that individual
emp_d_ta <- df_scl %>% 
  dplyr::select( ta_ ) %>% 
  #use the amt fit_dist function
  amt::fit_distr(., dist_name = "vonmises" )

#view
emp_d_sl
emp_d_ta

# update sl distribution parameters
b_log_sl <- fix.efs[5]
b_sl <-  fix.efs[4]
#update sl distribution
updated_sl <- update_gamma( emp_d_sl,
                              beta_sl = b_sl,          
                              beta_log_sl = b_log_sl )
#update ta distribution parameters
b_costa <- fix.efs[6]
#update turning angle distribution:
updated_ta <- update_vonmises( emp_d_ta,
                                 beta_cos_ta = b_costa )

# view results
updated_sl
# Are any of the parameters negative? If so then the model is ill fitted. 
# Tav Avgar recommends to try a different step-length distribution
# include different interractions 
# remove non-movement steps (based on a step-length threshold )
# resample data to coarser resolution

updated_ta
#Is the von Mises concentration parameter (kappa) negative? #
# If so Tal Avgar indicates that the adjusted turn angle distribution
# is centred at pi (180) rather than 0, meaning that the animal is 
# more likely to turn back. 

###### plot changes in step length distribution ###
# data.frame for plotting
plot_sl <- data.frame(x = rep(NA, 100))

#check step lenghts of chosen individual to choose your x
hist(df_scl[ which(df_scl$id == 2),'sl_'])

# x-axis is sequence of possible step lengths
plot_sl$x <- seq(from = 0, to = 10, length.out = 100)

# y-axis is the probability density under the given gamma distribution
# For the updated distribution when habitat is low
plot_sl$tentative <- dgamma(
  x = plot_sl$x,
  shape = emp_d_sl$params$shape,
  scale = emp_d_sl$params$scale)
#when habitat is high
plot_sl$updated <- dgamma(
  x = plot_sl$x,
  shape = updated_sl$params$shape,
  scale = updated_sl$params$scale)

# Pivot from wide data to long data
plot_sl <- plot_sl %>% 
  pivot_longer(cols = -x)

tail(plot_sl)
# Plot
ggplot(plot_sl, aes(x = x, y = value, color = factor(name))) +
  geom_line(size = 1) +
  xlab("Step Length (m)") +
  ylab("Probability Density") +
  scale_color_manual(name = "Distribution", 
                     breaks = c("tentative", "updated"),
                     values = c("blue", "orange")) +
  theme_bw()

#How did the distribution change with model results?
# Answer:
#

############### save section ######################################
#save workspace if in progress
save.image( 'iSSF_refs_results.RData'  )

################## end of script #####################################