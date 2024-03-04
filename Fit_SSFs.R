##################################################################
# Script developed by Jen Cruz to estimate SSFs and iSSFs          #
# approach derived from Fieberg et al. 2021 and Signer et al. 2019 #
# using code from Appendices B and C                             #
# also Muff et al. 2019 DOI: 10.1111/1365-2656.13087             #
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

################## prep workspace ###############################
#install relevant packages
# #install.packages( "peakRAM" )
# install.packages("INLA",repos=c(getOption("repos"),
#   INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

# load packages relevant to this script:
library( tidyverse ) #easy data manipulation
# set option to see all columns and more than 10 rows
options( dplyr.width = Inf, dplyr.print_min = 100 )
library( amt )
library( glmmTMB ) # for analysis
#library( INLA )
#to enable maximum RAM
#library( peakRAM )

#####################################################################
## end of package load ###############

###################################################################
#### Load or create data -----------------------------------------

# Clean your workspace to reset your R environment. #
rm( list = ls() )
#load 30m steps estimated for all individuals and habitat 
# variables extracted for each step
df_steps <- read_rds( "Data/df_steps" )
#view
head( df_steps )

#import polygon of the NCA as sf spatial file:
#NCA_Shape <- sf::st_read("Z:/Common/QCLData/Habitat/NCA/GIS_NCA_IDARNGpgsSampling/BOPNCA_Boundary.shp")
#######################################################################
######## preparing data ###############################################
#view data
head( df_steps)
#create vector of predictors taking advantage of naming commonality 
# to automatically extract them:
prednames <- grep('0m', colnames(df_steps), value = TRUE)

#check for missing values
colSums( is.na( df_steps[,prednames] ) )
#check for correlation but also whether the values sum to one or close 
# to. See: https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.4256
# for reasons why that is an issue.
# we have missing values so we start by removing them 
preddf <- df_steps[!is.na(df_steps[prednames[1]]), prednames]
preddf <- preddf[!is.na(preddf[prednames[2]]), ]
preddf <- preddf[!is.na(preddf[prednames[3]]), ]
preddf <- preddf[!is.na(preddf[prednames[4]]), ]
preddf <- preddf[!is.na(preddf[prednames[5]]), ]
preddf <- preddf[!is.na(preddf[prednames[6]]), ]
dim(preddf); dim(df_steps);head(preddf)

corrplot::corrplot( round(cor( preddf ) , 1 ), method = "number" )
#sum covariates at appropriate scale and then plot results:
#for 30m
hist(apply( preddf[ ,prednames[1:3] ], 1,sum ))
#for 100m
hist(apply( preddf[ ,prednames[4:6] ], 1,sum ))

#now check whether there is evidence of individual differences 
#in habitat selection 
ggplot( df_steps ) +
  theme_bw( base_size = 15 ) +
  geom_density( aes( x = perennial_30m, 
                     fill = case_, group = case_ ),
                alpha = 0.5  )  +
  facet_wrap( ~ territory )
ggplot( df_steps ) +
  theme_bw( base_size = 15 ) +
  geom_density( aes( x = annual_30m, 
                     fill = case_, group = case_ ),
                alpha = 0.5  )  +
  facet_wrap( ~ territory )
ggplot( df_steps ) +
  theme_bw( base_size = 15 ) +
  geom_density( aes( x = shrub_30m, 
                     fill = case_, group = case_ ),
                alpha = 0.5  )  +
  facet_wrap( ~ territory )

#For homework check differences at the other spatial scale
# What did you find?
# Answer: 
#
# Scale predictors 
#create new dataframe to hold scaled predictors, while keeping 
# unscaled ones for plotting later
df_scl <- df_steps

#scale only those columns:
df_scl[, prednames] <- apply( df_scl[,prednames], 2, scale )
#view
head( df_scl)
#now check for missing values
colSums( is.na( df_scl[,prednames] ) )
#replace missing values with 0
df_scl$annual_30m[ is.na(df_scl$annual_30m) ] <- 0
df_scl$annual_100m[ is.na(df_scl$annual_100m) ] <- 0
df_scl$perennial_30m[ is.na(df_scl$perennial_30m) ] <- 0
df_scl$perennial_100m[ is.na(df_scl$perennial_100m) ] <- 0
df_scl$shrub_30m[ is.na(df_scl$shrub_30m) ] <- 0
df_scl$shrub_100m[ is.na(df_scl$shrub_100m) ] <- 0

# we also assign weights to available points to be much greater than used points
df_scl$weight <- 1000 ^( 1 - as.integer(df_scl$case_ ) )
#check
head( df_scl )

#### end data prep #############
###########################################################################
##### analyse data  ##########

# amt function doesn't take random effects or weights so we move to a more
# flexible package. We move straight into a model #
# that includes random intercepts and slopes as well as a fixed #
# large variance for the random intercepts, as recommended by Muff #
# et al. 2019 and weights of 1000 for available points #

#df_one <- df_scl %>% filter( id == 2 )
# we start by defining the model without running it, which let's us
# fit the large variance to the random ID intercepts
#p1 <- peakRAM( 
m1 <- glmmTMB( case_ ~ 0 + annual_30m + perennial_30m + 
                                shrub_30m +
             #use step_id as random effect to simulate conditional likelihood
             ( 1| step_id_ ), 
             family = poisson, data = df_scl,
             weights = weight,
              start = list(theta = log( 1e3 ) ), 
             map = list( theta = factor( c(NA) ) ) )  

summary( m1 )

# rerun model with random effects that account for individual
# differences in habitat selection
m2 <- glmmTMB( case_ ~  -1 + annual_30m + perennial_30m + shrub_30m +
               #define random effects
               ( 1 | step_id_ ) + 
               ( 0 + annual_30m  | id ) +
               ( 0 + perennial_30m  | id ) +
               ( 0 + shrub_30m | id ) ,
               #set family to Poisson
               family = poisson, data = df_scl,
               #define weights
               weights = weight, 
               #tell it not to change variance for step level
               map = list( theta = factor( c(NA, 1:3 ) ) ),
               #fix variance for step level random intercept
               start = list( theta = c( log( 1e3 ),0,0,0) )
                )

summary( m2 )

#now test the larger scale
m3 <- glmmTMB( case_ ~  -1 + annual_100m + perennial_100m + shrub_100m +
                 #define random effects
                 ( 1 | step_id_ ) + 
                 ( 0 + annual_100m  | id ) +
                 ( 0 + perennial_100m  | id ) +
                 ( 0 + shrub_100m | id ) ,
               #set family to Poisson
               family = poisson, data = df_scl,
               #define weights
               weights = weight, 
               #tell it not to change variance for step level
               map = list( theta = factor( c(NA, 1:3 ) ) ),
               #fix variance for step level random intercept
               start = list( theta = c( log( 1e3 ),0,0,0) )
)

summary( m3 )

#which scale was better??? 
# Answer:
#
##############################################################################
################### visualizing top model results ############
###################
# We start by visualizing population-level effects:
#pull out confidence intervals for your top model
cis <- confint( m2 )
#change to df
cis <- as.data.frame( cis )
#view
cis
#extract row names as column
cis$preds <- rownames( cis )
#edit column names
colnames(cis )[1:3] <- c("L", "H", "Mean" )
#keep only fixed effects 
cis <- cis[1:3, ]
cis
#plot fixed effects (population-level effects)
ggplot( data = cis, aes( x = preds, y  = exp(Mean) ) ) + 
  theme_classic( base_size = 16) +
  geom_point( size = 3 ) +
  geom_errorbar( aes( ymin = exp(L), ymax = exp(H) ), linewidth = 1 ) +
  labs(x = element_blank(), y = "Relative selection strengh") +
  geom_hline( yintercept = 1, linewidth = 1 )

# Now we extract random effects so we can assess differences between 
# individuals 
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
#create id column
rss$id <- as.numeric(  rownames( rss ) )
#view
rss
# now extract additional details from our steps dataframe to combine 
# with our results
iddf <- df_steps %>% 
  group_by( id, territory, sex ) %>% 
  summarise( annual_30m_mean = mean( annual_30m, na.rm = TRUE),
             perennial_30m_mean = mean( perennial_30m, na.rm = TRUE),
             shrub_30m_mean = mean( shrub_30m, na.rm = TRUE)
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

###########################################################
##### end ######
##########################################################################
### Save desired results   #
#we save the scaled dataframe so that we can use it for the iSSFs
write_rds( df_scl, "Data/df_scl"  )
#save workspace if in progress
save.image( 'SSF_results.RData'  )

############# end of script  ##################################