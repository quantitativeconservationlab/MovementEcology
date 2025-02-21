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

#load our clean data frame to assess 1st order selection
df_sa <- read.csv( "Data/df_sa.csv" )
#load data for 2nd order selection
df_hr <- read.csv( "Data/df_hr.csv" )
#import polygon of the NCA as sf spatial file:
NCA_Shape <- sf::st_read("Z:/Common/QCLData/Habitat/NCA/GIS_NCA_IDARNGpgsSampling/BOPNCA_Boundary.shp")
#this one has the same CRS as our used/available points
#######################################################################
######## preparing data ###############################################
#create vector of predictors taking advantage of naming commonality 
# to automatically extract them:
prednames <- grep('0m', colnames(df_sa), value = TRUE)
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
idnos <- sort( unique( df_scl$territory )) 
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
# We want to determine use within the NCA assuming 10 individuals #
# is a representative sample. When would this be the case? #

# When would it not be the case? #
# Answer:
#

# We start with our finest resolution of predictors at 30 x 30 m cells:
msa_30m <- fit_rsf( case_ ~  1 + annual_30m + perennial_30m +
                     shrub_30m, data = sa_scl ) 

#view results
summary( msa_30m )

#as is best practice we add weights. Note that the weights cannot be 
# added to fit_rsf() function so we shift to the more flexible
# glmmTMB which allows weights, random effects and different distributions
msa_30m <- glmmTMB( case_ ~  1 + annual_30m + perennial_30m +
                      shrub_30m, 
                    data = sa_scl, 
        #we define binomial distribution for the response and add weights
            family = binomial(), weights = weight ) 

#view results
summary( msa_30m )

#now that we are happy with the approach we replicate for coarser scales
msa_100m <- glmmTMB( case_ ~ 1 + annual_100m + perennial_100m +
                      shrub_100m,
                    family = binomial(), data = sa_scl, 
                    weights = weight ) 
#view results
summary( msa_100m )
#Did the estimates for the coefficients change at this coarser scale?
# Answer:
#
# Now we look at 500m
msa_500m <- glmmTMB( case_ ~ 1 + annual_500m + perennial_500m +
                      shrub_500m,
                   family = binomial(), data = sa_scl, 
                   weights = weight ) 

summary( msa_500m )
#What about at our biggest scale?
# Answer:
#

#which scale has the most support? We compare model fit using AIC
anova( msa_100m, msa_30m, msa_500m)

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

# So a Prarie will be 1.45 times more likely to use an area with 10 %
# shrub than 2.5 % shrub when annual is 13. 6% and perennial is 19 %

# we also plot differences in distribution between used and available #
# locations for our predictor of choice. To plot on the real scale we #
# combine unscaled data first:
ggplot( df_sa ) +
  theme_bw( base_size = 15 ) +
  geom_density( aes( x = shrub_100m, 
                     fill = case_, group = case_ ),
                alpha = 0.5  ) #+
  #facet_wrap( ~ territory )
#compare number of points 
table( df_all$territory)

ggplot( df_sa ) +
  theme_bw( base_size = 15 ) +
  geom_density( aes( x = perennial_500m, 
                     fill = case_, group = case_ ),
                alpha = 0.5  )  

ggplot( df_sa ) +
  theme_bw( base_size = 15 ) +
  geom_density( aes( x = annual_500m, 
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
# We replicate our approach at the home-range scale#
#start with finest scale 
mhr_30m <- glmmTMB( case_ ~  1 + annual_30m + perennial_30m +
                      shrub_30m, 
                    data = hr_scl, 
                    #we define binomial distribution for the response and add weights
                    family = binomial(), weights = weight ) 

#view results
summary( mhr_30m )
# we exponentiate coefficients for easier interpretation:
exp( glmmTMB::fixef( mhr_30m )$cond )

# We follow with 100 m scale:
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
anova( mhr_100m, mhr_30m, mhr_500m)

# Which scale is most supported by model selection?
# Answer:
#

# Interpret results of the top model ###
exp( glmmTMB::fixef( mhr_500m )$cond )
# To remind ourselves what the SD for our predictor is
apply( df_hr[,prednames], 2, sd )
# And now the mean values for each habitat:
apply( df_hr[,prednames], 2, mean )

#How do you interpret selection for the most used predictor?
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