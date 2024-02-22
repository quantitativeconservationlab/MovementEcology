##################################################################
# Script developed by Jen Cruz to estimate resource selection     #
# functions. Code adapted from atm vignette:                      #
#https://cran.r-project.org/web/packages/amt/vignettes/p3_rsf.html #
# Vegetation cover  was downloaded from Rangeland Analysis Platform #
# https://rangelands.app/products/ for 2021 and includes        #
# % cover for shrub, perennial herbaceous, annual herbaceous    #
# tree, litter and bare ground                                   #
# coordinate system is WGS84 EPSG:4326, spatial resolution is 30m #
# We perform analyses in amt (Signer et al. 2019) and glmmTMB    #
# following Muff et al. (2019)                                  #
# go here: https://conservancy.umn.edu/handle/11299/204737      #
# for detailed code used by Muff.                               #
# Prairie Falcon data was thinned to 30minutes for 9 individuals #
# tracked in 2021.                                               #
###################################################################

################## prep workspace ###############################
# we will be using new packages:
install.packages( "glmmTMB" )

# load packages relevant to this script:
library( tidyverse ) #easy data manipulation
library( amt )
library( glmmTMB ) # for analysis

#####################################################################
## end of package load ###############

###################################################################
#### Load or create data -----------------------------------------
# Clean your workspace to reset your R environment. #
rm( list = ls() )

#load the thinned (30min) data for all individuals
df_all <- read_rds( "Data/df_all" )
#load ranges for all individuals:
#akde_all <- read_rds( "Data/akde_all" )

#import polygon of the NCA as sf spatial file:
#NCA_Shape <- sf::st_read("Z:/Common/QCLData/Habitat/NCA/GIS_NCA_IDARNGpgsSampling/BOPNCA_Boundary.shp")

#######################################################################
######## preparing data ###############################################

###############  Single individual Example ##################
###########
#extract data for individual of interest:
df_one <- df_all %>% filter( id == 4 )
head( df_one)
#create vector of predictors taking advantage of naming commonality 
# to automatically extract them:
prednames <- grep('0m', colnames(df_one), value = TRUE)
#check for correlation but also whether the values sum to one or close 
# to. See: https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.4256
# for reasons why that is an issue.
corrplot::corrplot( round(cor(df_one[prednames]),1 ), method = "number" )
#sum covariates at appropriate scale and then plot results:
#for 30m
hist(apply( df_one[ ,prednames[1:3] ], 1,sum ))
#for 500m
hist(apply( df_one[ ,prednames[4:6] ], 1,sum ))

# Scale predictors 
#create new dataframe to hold scaled predictors, while keeping 
# unscaled ones for plotting later
df_scl <- df_one

#scale only those columns:
df_scl[, prednames] <- apply( df_scl[,prednames], 2, scale )
#view
head( df_scl)
#now check for missing values
colSums( is.na( df_scl[,prednames] ) )
#no missing values in this instance. 
# we also assign weights to available points to be much greater than used points
df_scl$weight <- 1000 ^( 1 - as.integer(df_scl$case_ ) )
#check
head( df_scl )


##### analyse data  ##########
#We can use fit_rsf, which is just a wrapper around 
#stats::glm with family = binomial(link = "logit").

# starting with single individual:
m1 <- df_scl %>% 
  fit_rsf( case_ ~ 0 + annual_30m + perennial_30m +
                          shrub_30m ) %>% 
  summary()
# Based on what you have learnt, what is missing from this analysis #
# that may be biasing results?
# Answer:
# 
# we rerun the same model but using glmmTMB to make sure our results 
# are comparable
m1.1 <- glmmTMB( case_ ~  0 + annual_30m + perennial_30m +
                   shrub_30m,
                 family = binomial(), data = df_scl ) 
#view
summary( m1.1 )

# Are the results comparable?
# Answer:
# 
# We now include weights in the model
m1.w <- glmmTMB( case_ ~  0 + annual_30m + perennial_30m +
                   shrub_30m,
                 family = binomial(), data = df_scl, 
                 weights = weight ) 
#view
summary( m1.w )

# we repeat the weighted model with the larger area
m2.w <- glmmTMB( case_ ~  0 + annual_500m + perennial_500m +
                   shrub_500m,
                 family = binomial(), data = df_scl, 
                 weights = weight ) 
#view
summary( m2.w )

# computation was fairly fast for this individual so we move on to a 
#population-level analyses

#for homework try a different individual. 
# Interpret results here:
# Answer:
#

##### end single indiv analyses ####
#########################################################################
#################### Population-level RSF #########################
# We want to determine use within the NCA assuming 10 individuals #
# is a representative sample. When would this be the case? #
# When would it not be the case? #

#check for correlation but also whether the values sum to one or close 
# to. See: https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.4256
# for reasons why that is an issue.
corrplot::corrplot( round(cor(df_all[prednames]),1 ), method = "number" )
#sum covariates at appropriate scale and then plot results:
#for 30m
hist(apply( df_all[ ,prednames[1:3] ], 1,sum ))
#for 500m
hist(apply( df_all[ ,prednames[4:6] ], 1,sum ))

# Scale predictors 
#create new dataframe to hold scaled predictors, while keeping 
# unscaled ones for plotting later
df_scl <- df_all
#scale only those columns:
df_scl[, prednames] <- apply( df_scl[,prednames], 2, scale )
#view
head( df_scl)
#now check for missing values
for( i in 1:length(prednames)){
  a <- sum( is.na( df_scl[, prednames[i]] ) )
  print( prednames[i])
  print(a)
}
#no missing values in this instance. 
# we also assign weights to available points to be much greater than used points
df_scl$weight <- 1000 ^( 1 - as.integer(df_scl$case_ ) )
#check
head( df_scl )

#remember available points were extracted within each individual's 
#range

#extract individual id numbers:
idnos <- sort( unique( df_scl$territory )) 
# Start with amt package  using the NCA as available habitat
mp1 <- df_scl %>%  amt::fit_rsf( case_ ~ sage_30m ) %>% 
            summary()

#amt function doesn't take weights so we move to a more flexible #
# package for more accurate analyses
# We focus on comparing models between our two scales:
mp_30m <- glmmTMB( case_ ~  0 + annual_30m + perennial_30m +
                     shrub_30m,
                   family = binomial(), data = df_scl, 
                weights = weight ) 

summary( mp_30m )

# run analysis for the coarser scale:
mp_500m <- glmmTMB( case_ ~ 0 + annual_500m + perennial_500m +
                      shrub_500m,
                   family = binomial(), data = df_scl, 
                   weights = weight ) 

summary( mp_500m )

#which scale has the most support? How do we choose?
# Answer:
#

# Interpreting results ###
# we start by exponentiating the coefficients:
exp( glmmTMB::fixef( mp_500m )$cond )
# this reflects the relative selection strength for choosing each
# vegetation cover when the remaining vegetation covers are kept 
# at their mean values
# Thus prairie falcons are 1.05 times more likely to choose
# shrub with cover that is 1 SD higher when annual and perennial are 
# kept at their mean

# To remind ourselves what the SD for our predictor is
apply( df_all[,prednames], 2, sd )
# And now the mean values for each habitat:
apply( df_all[,prednames], 2, mean )

# So a Prarie will be 1.05 times more likely to use an area with 17 %
# shrub than 9 % shrub when annual is 13% and perennial is 24 %

# we also plot differences in distribution between used and available #
# locations for our predictor of choice. To plot on the real scale we #
# combine unscaled data first:
ggplot( df_all ) +
  theme_bw( base_size = 15 ) +
  geom_density( aes( x = shrub_500m, 
                     fill = case_, group = case_ ),
                alpha = 0.5  ) +
  facet_wrap( ~ territory )

ggplot( df_all ) +
  theme_bw( base_size = 15 ) +
  geom_density( aes( x = annual_500m, 
                     fill = case_, group = case_ ),
                alpha = 0.5  )  +
  facet_wrap( ~ territory )

ggplot( df_all ) +
  theme_bw( base_size = 15 ) +
  geom_density( aes( x = perennial_500m, 
                     fill = case_, group = case_ ),
                alpha = 0.5  )  +
  facet_wrap( ~ territory )


#What do these plots tell us about indiviudual differences?
# Is it reasonable to assume that all individuals are selecting 
# habitat similarly?
# Answer:
# 


###########################################################
### Save desired results                                  #
# we can save the movement model results
#save workspace if in progress
save.image( 'RSFresults.RData' )
############# end of script  ##################################