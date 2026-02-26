##################################################################
# Script developed by Jen Cruz to estimate resource selection     #
# functions. Code adapted from atm vignette:                      #
#https://cran.r-project.org/web/packages/amt/vignettes/p3_rsf.html #
# Vegetation cover  was downloaded from Rangeland Analysis Platform #
# https://rangelands.app/products/ for 2021 and includes        #
# % cover for shrub, perennial herbaceous, annual herbaceous    #
# tree, litter and bare ground                                   #
# coordinate system is WGS84 EPSG:4326,                          #
# We perform analyses in amt (Signer et al. 2019) and glmmTMB    #
# Prairie Falcon data was thinned to 30minutes for 9 individuals #
# tracked in 2021.                                               #
###################################################################

################## prep workspace ###############################
# we will be using new packages:
install.packages( "glmmTMB" )
# load packages relevant to this script:
library( tidyverse ) #easy data manipulation
# set option to see all columns and more than 10 rows
options( dplyr.width = Inf, dplyr.print_min = 100 )
library( amt )
library( glmmTMB ) # for analysis

#####################################################################
## end of package load ###############

###################################################################
#### Load or create data -----------------------------------------
# Clean your workspace to reset your R environment. #
rm( list = ls() )

#if you already started the analysis then load workspace:
#load('RSFresults.RData')
#load our clean data frame to assess 1st order selection
df_sa <- read.csv( "Data/df_sa.csv" )
#load data for 2nd order selection
df_hr <- read.csv( "Data/df_hr.csv" )
#import polygon of the NCA as sf spatial file:
NCA_Shape <- sf::st_read("Data/BOPNCA_Boundary.shp")
#this one has the same CRS as our used/available points
#######################################################################
######## preparing data ###############################################
#create vector of predictors taking advantage of naming commonality 
# to automatically extract them:
prednames <- grep('0m', colnames(df_sa), value = TRUE)
#view
prednames
# Scale predictors create new dataframes to hold scaled predictors, while keeping 
# unscaled ones for plotting later
sa_scl <- df_sa
#scale only those columns:
sa_scl[, prednames] <- apply( sa_scl[,prednames], 2, scale )
#view
head( sa_scl)
# why do we scale predictors?
# Answer:
#

#now check for missing values
colSums( is.na( sa_scl[,prednames] ) )
#no missing values in this instance. 
# we also assign weights to available points to be much greater than used points
sa_scl$weight <- 1000 ^( 1 - as.integer( sa_scl$case_ ) )
#check
head( sa_scl )

# We repeat the process for second order selection where
# available points were extracted within each individual's range
#extract individual id numbers:
idnos <- sort( unique( df_hr$territory )) 
#duplicate dataframe
hr_scl <- df_hr
#scale only those columns:
hr_scl[, prednames] <- apply( hr_scl[,prednames], 2, scale )
#view
head( hr_scl)
#now check for missing values
colSums( is.na( hr_scl[,prednames] ) )
#no missing values in this instance. 
# we also assign weights to available points to be much greater than used points
hr_scl$weight <- 1000 ^( 1 - as.integer( hr_scl$case_ ) )
#check
head( hr_scl )

#########################################################################
#################### Population-level RSF #########################
# We want to determine use within the NCA assuming 9 prairie falcons #
# are a representative sample of the population. #

# When would this not be the case? #
# Answer:
#

# We start with our finest resolution of predictors:
msa_100m <- glmmTMB( case_ ~ 1 + annual_100m + perennial_100m +
                      shrub_100m,
                    family = binomial(), data = sa_scl, 
                    weights = weight ) 
#view results
summary( msa_100m )
# How do you interpret the summary table of results?
# Answer:
#
# Now we look at 500m resolution
msa_500m <- glmmTMB( case_ ~ 1 + annual_500m + perennial_500m +
                      shrub_500m,
                   family = binomial(), data = sa_scl, 
                   weights = weight ) 

summary( msa_500m )
#Did inference change at the bigger scale?
# Answer:
#

#which scale has the most support? We compare model fit using AIC
anova( msa_100m, msa_500m)

# Which scale is most supported by model selection?
# Answer:
#

# Interpreting results of the top model ###
# we start by exponentiating the coefficients:
exp( glmmTMB::fixef( msa_100m )$cond )
# this reflects the relative selection strength for choosing each
# vegetation cover when the remaining vegetation covers are kept 
# at their mean values
# Thus prairie falcons are 1.45 times more likely to choose
# shrub with cover that is 1 SD higher when annual and perennial are 
# kept at their mean

# To remind ourselves what the SD for our predictor is
apply( df_sa[,prednames], 2, sd )
# And now the mean values for each habitat:
apply( df_sa[,prednames], 2, mean )

# So a Prarie Falcon will be 1.4 times more likely to use an area with 10 %
# shrub than 2.8 % shrub when annual is 19% and perennial is 13.5 %

#What would be the equivalent statement for annual?
# Answer: 
#

# plot differences in distribution between used and available #
# locations for our predictor of choice. To plot on the real scale we #
# combine unscaled data first:
ggplot( df_sa ) +
  theme_bw( base_size = 15 ) +
  geom_density( aes( x = shrub_100m, 
                     fill = case_, group = case_ ),
                alpha = 0.5  ) 
#now for perennial
ggplot( df_sa ) +
  theme_bw( base_size = 15 ) +
  geom_density( aes( x = perennial_100m, 
                     fill = case_, group = case_ ),
                alpha = 0.5  )  
#for annual:
ggplot( df_sa ) +
  theme_bw( base_size = 15 ) +
  geom_density( aes( x = annual_100m, 
                     fill = case_, group = case_ ),
                alpha = 0.5  )  

#What do these plots tell us about how prairie falcons select habitat?
# Is it reasonable to assume that all individuals are selecting 
# habitat similarly?
# Answer:
# 
#########
################ 2nd order RSFs ####################
######
# We replicate our approach at the home-range scale starting with 100m#
mhr_100m <- glmmTMB( case_ ~ 1 + annual_100m + perennial_100m +
                       shrub_100m,
                     family = binomial(), data = hr_scl, 
                     weights = weight ) 
#view results
summary( mhr_100m )
exp( glmmTMB::fixef( mhr_100m )$cond )

# Now we look at 500m
mhr_500m <- glmmTMB( case_ ~ 1 + annual_500m + perennial_500m +
                       shrub_500m,
                     family = binomial(), data = hr_scl, 
                     weights = weight ) 

summary( mhr_500m )
exp( glmmTMB::fixef( mhr_500m )$cond )

#We compare models using AIC
anova( mhr_100m, mhr_500m)

# Which scale is most supported by model selection?
# Answer:
#
# Is the scale the same as that chosen for the study area analysis?
# If not, why would that be?
# Answer:
# 


# Interpret results of the top model ###
exp( glmmTMB::fixef( mhr_500m )$cond )
# To remind ourselves what the SD for our predictor is
apply( df_hr[,prednames], 2, sd )
# And now the mean values for each habitat:
apply( df_hr[,prednames], 2, mean )

#Using the same approach as we did for study area rephrase selection
# for the most significant predictor?
# Answer:
# 
#

# we also plot differences in distribution between used and available #
# locations for our predictor of choice. This time we look at #
# potential differences among individual ranges vs used habitat: 
ggplot( df_hr ) +
  theme_bw( base_size = 15 ) +
  geom_density( aes( x = shrub_500m, 
                     fill = case_, group = case_ ),
                alpha = 0.5  ) +
    facet_wrap( ~ territory )

# Are all individuals using shrub in higher proportions than what 
# is available inside their range? 
# Describe and contrast selection for each individual:
# Answer: 
# 
#

# Are their ranges filled with similar amounts
# Describe differences here:
# 
# 
#
# Tally individuals using more shrub cover than what 
# is available in their range:
# Answer:
# 
#

# For homework also add similar figures for the other
# two vegetation types at the same scale of the top model
# and interpret differences and similarities
# Add code and responses here:
#
#

# Did the interpretation of which habitats are selected by #
# Prairie falcons differ between the 1st order and 2 order selection?
# Answer:
#



###########################################################
### Save desired results                                  #
# we can save the movement model results
#save workspace if in progress
save.image( 'RSFresults.RData' )
############# end of script  ##################################