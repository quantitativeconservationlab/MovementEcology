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
#install relevant packages
install.packages( "circular" )

# load packages relevant to this script:
library( tidyverse ) #easy data manipulation
# set option to see all columns and more than 10 rows
options( dplyr.width = Inf, dplyr.print_min = 100 )
library( amt )
library( circular ) #for plotting von mises distribution

#####################################################################
## end of package load ###############

###################################################################
#### Load or create data -----------------------------------------

# Clean your workspace to reset your R environment. #
rm( list = ls() )
#load 30m steps estimated for all individuals and habitat 
# variables extracted for each step
df_steps <- read.csv( "Data/df_steps30.csv" )

#######################################################################
######## preparing data ###############################################
#view data
head( df_steps)
#create vector of potential predictors
prednames <- c( "annual", "perennial", "shrub" )

#check for missing values
colSums( is.na( df_steps[,prednames] ) )

# Scale predictors 
#create new dataframe to hold scaled predictors, while keeping 
# unscaled ones for plotting later
df_scl <- df_steps

#scale only those columns:
df_scl[, prednames] <- apply( df_scl[,prednames], 2, scale )
#view
head( df_scl)
#now replace small amount of missing values with 0,which represents 
# the mean for scaled predictors
df_scl$annual[is.na(df_scl$annual)] <- 0
df_scl$perennial[is.na(df_scl$perennial)] <- 0
df_scl$shrub[is.na(df_scl$shrub)] <- 0

#for the movement parameters we also calculate the log of sl and the cos of ta"
df_scl$log_sl_ <- log( df_scl$sl_ )
df_scl$cos_ta_ <- cos( df_scl$ta_ )
# we also turn our step lengths to km instead of meters
df_scl$sl_ <- df_scl$sl_ / 1000
#check
hist(df_scl$sl_ )
# we also assign weights to available points to be much greater than used points
df_scl$weight <- 1000 ^( 1 - as.integer(df_scl$case_ ) )
#check
head( df_scl )

#### end data prep #############
###########################################################################
#####   running STEP SELECTION FUNCTIONS          ##########
# We start with a traditional step selection function
#use step_id as random intercept to account for conditional likelihood
# we fit separate models for each individual
mi1 <- df_scl %>% dplyr::filter( id == 1 ) %>% 
  fit_ssf( #response variable
          case_ ~ 
            #habitat variable
            annual + perennial + shrub +
            #stratum to ensure random steps match to each point
             strata( step_id_ ), 
          model = TRUE )   

summary( mi1 )

#now we fit the same model for all individuals:
mall <- df_scl %>% 
    nest( data = -id ) %>% 
  dplyr::mutate( ssf = lapply( data, function(x) {
    x %>%  amt::fit_ssf(  case_ ~ annual + perennial + shrub +
      strata( step_id_ ) )
  } ) )

mall
#we clean up and combine results to get the average selection
# across all individuals
d2 <- mall %>% 
  dplyr::mutate( coef = map( ssf, 
  ~broom::tidy(.x$model) ) ) %>% 
  dplyr::select( id, coef ) %>% 
  unnest( cols = c(coef) ) %>% 
  dplyr::mutate( id = factor(id) ) %>%
  dplyr::group_by( term )%>% 
  dplyr::summarize( 
    mean = mean( estimate ), 
    #calculate 95% CIs
    ymin = mean - 1.96 *sd(estimate), 
    ymax = mean + 1.96 *sd(estimate) )

d2$x <- 1:nrow( d2 )

# visualizing model results #
# extract coefficients for each individual
coefsall <- mall %>% 
  dplyr::mutate( coef = map( ssf, 
                             ~broom::tidy(.x$model) ) ) %>% 
  dplyr::select( id, coef ) %>% 
  unnest( cols = c(coef) ) %>% 
  dplyr::mutate( id = factor(id),
                 conf.low = estimate - 1.96 * std.error,
                 conf.high = estimate + 1.96 * std.error )
head(coefsall)

#we plot individual differences 
pall <- coefsall %>%
  ggplot(., aes(x = term, y = estimate, 
                group = id, col = id ) ) +
  #add individual results
  geom_pointrange( aes( ymin = conf.low, 
                        ymax = conf.high ),
      position = position_dodge( width = 0.7 ), size = 0.8 ) +
  #draw line at 0
  geom_hline( yintercept = 0, lty = 2 ) +
  #start with population level averages we calculated earlier
  geom_rect( mapping = aes(xmin = x - 0.4, xmax = x + 0.4, 
                          ymin = ymin, ymax = ymax ), 
             data = d2, 
            inherit.aes = FALSE, fill = "grey90", alpha = 0.5) +
  geom_segment(mapping = aes(x = x - 0.4, xend = x + 0.4, 
                             y = mean, yend = mean ), 
               data = d2, inherit.aes = FALSE, size = 1 ) +
  #Add the labels to each axis
  labs(x = "Habitat", y = "Relative Selection Strength") + 
  theme_light()

pall
# How do you interpret the results from this figure?
# Answer:
# 
#

# could selection be due to the amount of habitat available for each 
# individual?
# To answer this question we extract additional details from our steps
# dataframe including sex and the average amount of veg cover available
# for each individual

#we create a new id df with those summary values
iddf <- df_steps %>% 
  group_by( id, territory, sex ) %>% 
  summarise( annual_mean = mean( annual, na.rm = TRUE),
             perennial_mean = mean( perennial, na.rm = TRUE),
             shrub_mean = mean( shrub, na.rm = TRUE)
             )
iddf
#turn into a long format to combine with coefs dataframe
id_long <- iddf %>% 
  pivot_longer( cols = ends_with( "mean" ),
              names_to = "term",
              values_to = "cover" )
#view
head( id_long )
#modify term to match coefs by removing _mean from labels
id_long$term <- str_split_i( id_long$term, "_",1 )
#view
head( id_long )
#turn id to factor to match
id_long$id <- as.factor( id_long$id)
#combine with our resource selection strength estimates
coefs_df <- left_join( id_long, coefsall, by = c("id", "term" ) )
#view
head(coefs_df)

#plot resource selection strength by vegetation cover 
ggplot( coefs_df,aes( x = cover , y = estimate, color = id ) ) +
  theme_classic( base_size = 15 ) +
  labs( x = "Mean cover (%)", 
        y = "Resource selection strength" ) +
  geom_point() +
  geom_errorbar( aes( ymin = conf.low, 
                      ymax = conf.high ) ) +
  geom_hline( yintercept = 0, linewidth = 1, lty = 2 ) + 
  facet_wrap( ~term, scales = "free", ncol = 1 )


#### what do you interpret from this plot?
# is amount of vegetation influencing results?
# Answer:
#
# What about sex?
# Answer:
#

##### end of ssf analysis #####
############################################################
#### iSSF analysis                                  #####
#####################################################################
## We saw that there are differences in how individuals are selecting#
# habitat based on our previous analysis BUT we do not know yet #
# the relationship between habitat and how individuals move. To #
# explore those we shift to iSSFs using the same data            #
################## single individual iSSF ########################
# 
# For homework choose a different one by modifying code below:
mi <- df_scl %>% dplyr::filter( id == 2 ) %>% 
  fit_issf( #response variable
    case_ ~ 
      #add habitat variables
      annual + perennial + shrub +
      #add movement variables
      log_sl_ + cos_ta_ + sl_ +
      # add movement interactions with shrub
      log_sl_:shrub + cos_ta_:shrub +
       log_sl_:perennial + cos_ta_:perennial +
       log_sl_:annual + cos_ta_:annual +
      #add stratum to ensure random steps are matched to corresponding used step
      strata( step_id_ ), model = TRUE )

summary( mi )

# We calculate the tentative distributions from empirical data 
# for that same individual
# Start with step length fitted as a gamma with shape and scale parameters
emp_d_sl <- df_scl %>% 
  #select step lengths for that individual
  dplyr::filter( id == 2 ) %>% 
  dplyr::select( sl_ ) %>% 
  #fit a gamma distribution using empirical data
  amt::fit_distr(., dist_name = "gamma" )

#Fit a von misses to the turning angles for that individual
emp_d_ta <- df_scl %>% 
  #select turning angles for that individual
  dplyr::filter( id == 2) %>% 
  dplyr::select( ta_ ) %>% 
  #use the amt fit_dist function
  amt::fit_distr(., dist_name = "vonmises" )

#Assign the empirical distributions to model object:
mi$sl_ <- emp_d_sl
mi$ta_ <- emp_d_ta
# view
mi$sl_
mi$ta_

# Now we can use coefficients associated with movement parameters
# to update our movement related distributions for the same individual. #
# Refer to appendix c in Fieberg et al. 2021 to see the equations that
# are used to update distribution parameters #

# we need to first relabel coefficients since we have interactions
# Start by extracting coefficients of the model 
b <- coef( mi )
#choose significant interaction with habitat
# Here I choose annual.
# Change it depending on your individual results
# Modify code accordingly:
summary(mi)
b_log_l <- b["log_sl_"] 
b_log_h <- b["log_sl_"] + b["annual:log_sl_"] 
b_sl <-  b["sl_"] 
# Update step length distribution to the baseline when shrubs don't interact 
# with step length:
updated_sl_l <- update_gamma( mi$sl_, 
                  beta_sl = b_sl,          
                  beta_log_sl = b_log_l )
# Update step length distribution of how habitat alters step distribution
updated_sl_h <- update_gamma( mi$sl_, 
                              beta_sl = b_sl,          
                              beta_log_sl = b_log_h )

#view estimated parameters
updated_sl_l;updated_sl_h
# Are any of the parameters negative? If so then the model is ill fitted. 
# Tal Avgar recommends to try a different step-length distribution
# include different interactions 
# remove non-movement steps (based on a step-length threshold )
# resample data to coarser resolution

#For turning angle, choose significant interaction if present for that 
# individual 
# Modify code accordingly:
b_costa_l <- b["cos_ta_"]
b_costa_h <- b[ "cos_ta_" ] + b["perennial:cos_ta_"]
#update turning angle distribution:
updated_ta_l <- update_vonmises( mi$ta_,
                               beta_cos_ta = b_costa_l )
updated_ta_h <- update_vonmises( mi$ta_,
                                beta_cos_ta = b_costa_h )

#View results
updated_ta_l
updated_ta_h
#Is the von Mises concentration parameter (kappa) negative? #
# If so Tal Avgar indicates that the adjusted turn angle distribution
# is centred at pi (180) rather than 0, meaning that the animal is 
# more likely to turn back. 

# What results did you get for your individual? 
# Answer:
#

###### plot changes in step length distribution ###
# data.frame for plotting
plot_sl <- data.frame(x = rep(NA, 100))

#check step lenghts of chosen individual to choose your x
hist(df_scl[ which(df_scl$id == 2),'sl_'])

# x-axis is sequence of possible step lengths
plot_sl$x <- seq(from = 0, to = 10, length.out = 100)

# y-axis is the probability density under the given gamma distribution
# For the updated distribution when habitat is low
plot_sl$updated_l <- dgamma(
  x = plot_sl$x,
  shape = updated_sl_l$params$shape,
  scale = updated_sl_l$params$scale)
#when habitat is high
plot_sl$updated_h <- dgamma(
  x = plot_sl$x,
  shape = updated_sl_h$params$shape,
  scale = updated_sl_h$params$scale)

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
                     breaks = c("updated_l", "updated_h"),
                     values = c("blue", "orange")) +
  theme_bw()

#How did the distribution change with model results?
# Answer:
#

####  Plot turning angle distribution changes ###
# data.frame for plotting
plot_ta <- data.frame(x = rep(NA, 100))

# x-axis is sequence of possible step lengths
plot_ta$x <- seq(from = -1 * pi, to = pi, length.out = 100)

# y-axis is the probability density under the given von Mises distribution
# For low habitat
plot_ta$updated_l <- circular::dvonmises(
  x = plot_ta$x, 
  mu = updated_ta_l$params$mu,
  kappa = updated_ta_l$params$kappa)

# For high habitat
plot_ta$updated_h <- circular::dvonmises(
  x = plot_ta$x, 
  mu = updated_ta_h$params$mu,
  kappa = updated_ta_h$params$kappa)

# Pivot from wide data to long data
plot_ta <- plot_ta %>% 
  pivot_longer(cols = -x)

tail(plot_ta)
# Plot
ggplot(plot_ta, aes(x = x, y = value, color = factor(name))) +
  geom_line(size = 1) +
  coord_cartesian(ylim = c(0, 0.25)) +
  xlab("Relative Turn Angle (radians)") +
  ylab("Probability Density") +
  scale_x_continuous(breaks = c(-pi, -pi/2, 0, pi/2, pi),
                     labels = c(expression(-pi, -pi/2, 0, pi/2, pi))) +
  scale_color_manual(name = "Distribution", 
                     breaks = c("tentative", "updated"),
                     values = c("blue", "orange")) +
  theme_bw()

#How do you interpret output for these last two plots?
# Did habitat and movement parameters interact to alter how #
# Prairie falcons use their landscape ?
# Answer:
#
# Add answer for the new individual that you tried for homework
#
#
##### end of single individual issf ####
########################################################################
######### All individual issfs ####################################

# Repeat the analysis for all individuals at the same time #
miall2 <- df_scl %>% 
  nest( data = -id ) %>% 
  dplyr::mutate( issf = lapply( data, function(x) {
    x %>%  amt::fit_issf(  case_ ~ 
             #add habitat variables
             annual + perennial + shrub +
             #add movement variables
             log_sl_ + cos_ta_ + sl_ +
             # add movement interactions with shrub
             log_sl_:shrub + cos_ta_:shrub +
             log_sl_:perennial + cos_ta_:perennial +
             log_sl_:annual + cos_ta_:annual +
  #add stratum to ensure random steps are matched to corresponding used step
            strata( step_id_ ), model = TRUE )
  } ) )

miall2

# Extract coefficients for all indviduals
coefs_issf2  <- miall2 %>% 
  dplyr::mutate( coef = map( issf, 
                             ~broom::tidy(.x$model) ) ) %>% 
  dplyr::select( id, coef ) %>% 
  unnest( cols = c(coef) ) %>% 
  dplyr::mutate( id = factor(id),
                 conf.low = estimate - 1.96 * std.error,
                 conf.high = estimate + 1.96 * std.error )
#view results
coefs_issf2

#average across individuals
d4 <- miall2 %>% 
  dplyr::mutate( coef = map( issf, 
                             ~broom::tidy(.x$model) ) ) %>% 
  dplyr::select( id, coef ) %>% 
  unnest( cols = c(coef) ) %>% 
  dplyr::mutate( id = factor(id) ) %>%
  dplyr::group_by( term )%>% 
  dplyr::summarize( 
    mean = mean( estimate ), 
    #calculate 95% CIs
    ymin = mean - 1.96 *sd(estimate), 
    ymax = mean + 1.96 *sd(estimate) )

d4$x <- 1:nrow( d4 )

# Plot individual differences and population averages 
pissfs2 <- coefs_issf2 %>%
  ggplot(., aes(x = term, y = estimate, 
                group = id, col = id ) ) +
  #add individual results
  geom_pointrange( aes( ymin = conf.low, 
                        ymax = conf.high ),
                   position = position_dodge( width = 0.7 ), size = 0.8 ) +
  #draw line at 0
  geom_hline( yintercept = 0, lty = 2 ) +
  #start with population level averages we calculated earlier
  geom_rect( mapping = aes(xmin = x - 0.4, xmax = x + 0.4, 
                           ymin = ymin, ymax = ymax ), 
             data = d4, 
             inherit.aes = FALSE, fill = "grey90", alpha = 0.5) +
  geom_segment(mapping = aes(x = x - 0.4, xend = x + 0.4, 
                             y = mean, yend = mean ), 
               data = d4, inherit.aes = FALSE, size = 1 ) +
  #Add the labels to each axis
  labs(x = "Predictors", y = "Relative Selection Strength") + 
  theme_light()

pissfs2

### what are the weaknesses of averaging separate individual models #
# in this way? #
# Answer:
# 

# What is a posible solution to the weakness you stated above?
# Answer:
# 
##########################################################################
### Save desired results   #
#we save the scaled dataframe so that we can use it for our random effects
write.csv(df_scl, "Data/df_scl.csv"  )
#save workspace if in progress
save.image( 'SSF_results.RData'  )

############# end of script  ##################################