##################################################################
# Script developed by Jen Cruz to estimate resource selection     #
# functions. Code adapted from atm vignette:                      #
#https://cran.r-project.org/web/packages/amt/vignettes/p3_rsf.html #
# file:///G:/My%20Drive/Teaching/MovementEcology/vignettes/jane13441-sup-0001-AppendixA.html #
# but modified to remove the need to convert the predictor rasters #
# to a different CRS. Rasters were imported and manipulated using #
# the terra package, while sf was used for spatial objects.       #
# We use landcover data from the National Geospatial Data Asset  #
# https://www.mrlc.gov                                           #
# Habitat predictors include 2018 estimates of sagebrush cover  #
# We also perform analyses in glmmTMB following Muff et al.2019 #
# go here: https://conservancy.umn.edu/handle/11299/204737      #
# for detailed code discussed in the manuscript                 #
###################################################################

################## prep workspace ###############################

# Clean your workspace to reset your R environment. #
rm( list = ls() )
# we will be using new packages:
install.packages( "terra" )
install.packages( "glmmTMB" )
install.packages( "raster" )
install.packages( "rasterVis" )

# load packages relevant to this script:
library( tidyverse ) #easy data manipulation
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

#load the thinned (30min) data for all individuals
trks.thin <- read_rds( "Data/trks.thin" )
#load range for a single individual:
akde_amt <- read_rds( "Data/akde_auto" )
akde_all <- read_rds( "Data/akde_all" )

#import polygon of the NCA as sf spatial file:
NCA_Shape <- sf::st_read("Z:/Common/QCLData/Habitat/NCA/GIS_NCA_IDARNGpgsSampling/BOPNCA_Boundary.shp")

# define location of your raster file

# #import sagebrush raster:
sagebrush <- raster::raster( "Z:/Common/QCLData/Habitat/NLCD_new/TimeSeriesCover/Sagebrush/Sagebrush_2009_2020/rcmap_sagebrush_2020.img"  )
#view
sagebrush
# You can also plot it with the terra package:
terra::plot( sagebrush )
#######################################################################
######## preparing data ###############################################

#get coordinates from shapefile
crstracks <- sf::st_crs( NCA_Shape )
# checking outline of NCA
sf::st_bbox( NCA_Shape )
# We define available habitat as area of NCA with a small buffer #
# around it and draw points from it #
#create a buffer around the NCA using outline of NCA and sf package:
NCA_buf <- NCA_Shape %>% sf::st_buffer( dist =1e4 )
#create a version that matches coordinates of the predictor raster:
NCA_trans <- sf::st_transform( NCA_buf, st_crs( sagebrush ) ) 
#compare outline of trasformed polygon:
sf::st_bbox( NCA_trans )
#check extent
terra::ext( NCA_trans )
#crop raster to buffered NCA if you have a raster object:
sage_cropped <- raster::crop( sagebrush, NCA_trans )
# Now that we have cropped it to the appropriate area it should be faster #
# to process eventhough we are still using raster #
# compare extent with original
terra::ext(sage_cropped )
terra::ext( sagebrush )
#view
sage_cropped
#values greater than 100 are empty so replace with missing
sage_cropped[ sage_cropped > 100 ] <- NA
#Plot sagebrush
terra::plot( sage_cropped, main = "Sagebrush (2020)" )
#or 
rasterVis::levelplot( sage_cropped )
#Note all the missing data for that year!
#Plot density plot
terra::density( sage_cropped )
#Note that the max observed percentage of sagebrush < 40 %

# Note that transforming (projecting) raster data is fundamentally #
#different from transforming vector data. Vector data can be transformed #
#and back-transformed without loss in precision and without changes in #
#the values. This is not the case with raster data. In each transformation#
#the values for the new cells are estimated in some fashion. Therefore, #
#if you need to match raster and vector data for analysis, #
#you should generally transform the vector data. # 

##########################################################################
###############  Single individual Example ##################
###########
#extract data for individual of interest:
idv.thin <- trks.thin %>% filter( id == 2 )
# Extract available points within range of individual and #
# add used points from track:
r_one <- amt::random_points( akde_amt, n = (dim(idv.thin)[1]*50), 
                             presence = idv.thin )
#view
r_one
#plot
plot( r_one )
# How many available points were created?
# Answer:
# Is this enough?
# Answer: 
#
# Instead of using the amt::extract_covariates() function we convert our #
# r_one object to sf so that we can change the coordinate system and #
# extract values from the raster without loosing accuracy
class( r_one )
#convert to sf object defining coordinate column
r_one_sf <- sf::st_as_sf( r_one, coords = c("x_", "y_"), 
                          crs = crstracks )
#now transform to predictor crs:
r_one_trans <- sf::st_transform( r_one_sf, st_crs(sage_cropped) )

#now we can extract values from our predictor rasters while keeping them in #
# their original crs #
# As spatial scale we used original resolution of raster of 30 m grid cells

#if you were going to extract using terra, if your object was terra:
# #convert to dataframe for use in terra
# df_one <- r_one_trans %>%
#   dplyr::mutate( x = sf::st_coordinates(.)[,1],
#                  y = sf::st_coordinates(.)[,2]) %>% 
#   st_drop_geometry() %>% 
#   as.data.frame()
# #view
# head( df_one ); dim( df_one )
# sage_30m <- terra::extract( x = sage_cropped, y = df_one[ ,c("x","y") ],
#                            method = "simple" )

#extracting with raster we can used the sf object directly, you also 
# have the choice to use a buffer around each point if you want to increase 
# your resolution:
sage_30m <- raster::extract( x = sage_cropped, r_one_trans,
                            method = "simple" )

#check
sage_30m
# for 30 min intervals in the data, coarse resolutions may be appropriate
# We try using a buffer of 160m and extract the mean value of cells that #
# fall within that radius surrounding each point
# Trapping webs for ground squirrels have a 100 m radius, and the extra 50 m
# account for movements outside the web
sage_300m <- raster::extract( x = sage_cropped, r_one_trans,
                             method = "simple", buffer = 150,
                             fun  = mean, na.rm = TRUE )

#check
sage_300m
# What proportion of our data are missing values
sum( is.na( sage_30m ))/ length( sage_30m )
sum( is.na( sage_300m ))/ length( sage_300m )
# how do we deal with missing values in our analysis?
# Answer:
# 
# if your number of missing values is small, one approach is to use the mean #
# value, by averaging across your known values to extract the mean
# this becomes problematic if you have a lot of missing data. 
# What else can you do?
# Answer:
# 
# We append our predictor estimates to the original spatial tibble
df_one <- cbind( r_one, sage_30m, sage_300m )

# Scale predictors 
df_one$sage_30m <- scale( df_one$sage_30m )
df_one$sage_300m <- scale( df_one$sage_300m )
#replace missing values with mean, which is 0 after they have been scaled
df_one$sage_30m[ is.na(df_one$sage_30m) ] <- 0
df_one$sage_300m[ is.na(df_one$sage_300m) ] <- 0
# we also assign weights to available points to be much greater than used points
df_one$weight <- 1000 ^( 1 - as.integer(df_one$case_ ) )
#check
head( df_one )


##### analyse data  ##########
#We can use fit_rsf, which is just a wrapper around 
#stats::glm with family = binomial(link = "logit").

# starting with single individual:
m1 <- df_one %>% fit_rsf( case_ ~ sage_30m ) %>% 
  summary()
# Based on what you have learnt, what is missing from this analysis #
# that may be biasing results?
# Answer:
# 
# we rerun the same model but using glmmTMB to make sure our results 
# are comparable
m1.1 <- glmmTMB( case_ ~ sage_30m, family = binomial(), data = df_one ) 
#view
summary( m1.1 )

# Are the results comparable?
# Answer:
# 
# We now include weights in the model
m1.w <- glmmTMB( case_ ~ sage_30m, family = binomial(), data = df_one, 
                 weights = weight ) 
#view
summary( m1.w )

# we repeat the weighted model with sagebrush estimated over a larger area
m2.w <- glmmTMB( case_ ~ sage_300m, family = binomial(), data = df_one, 
                 weights = weight ) 
#view
summary( m2.w )

# computation was fairly fast for this individual so we move on to a 
#population-level analyses
##### end single indiv analyses ####
#########################################################################
#################### Population-level RSF #########################
# We want to determine use within the NCA assuming 10 individuals #
# is a representative sample. When would this be the case? #
# When would it not be the case? #

#use buffer to define available area and the tracks for used points:
# we specify how many available points we want
# but which is our available area? We start by choosing the entire
# NCA. What assumptions are we making with this choice?
r_all <- random_points( NCA_buf, n = (dim(trks.thin)[1] * 10 ), 
                        presence = trks.thin )
r_all
plot( r_all )

#let's check that the tracks fall inside the cropped area
spr <- sf::st_as_sf( trks.thin, coords = c("x_", "y_"), 
              crs = crstracks )
spr <- sf::st_transform( spr, st_crs(sage_cropped) )
spr <- as( st_geometry( spr), Class="Spatial")
rasterVis::levelplot( sage_cropped, margin = FALSE ) +
  latticeExtra::layer(sp.lines( spr, col = "yellow",
                                lwd = 3 ) )
  
#convert to sf object defining coordinate column
r_all_sf <- sf::st_as_sf( r_all, coords = c("x_", "y_"), 
                          crs = crstracks )
#now transform to predictor crs:
r_all_trans <- sf::st_transform( r_all_sf, st_crs(sage_cropped) )

#extract predictor at 30m resolution:
sage_all_30m <- raster::extract( x = sage_cropped, r_all_trans,
                             method = "simple" )

#check
sage_all_30m
#Now at 300m following process for one individual above
sage_all_300m <- raster::extract( x = sage_cropped, r_all_trans,
                              method = "simple", buffer = 150,
                              fun  = mean, na.rm = TRUE )

# We append our predictor estimates to the original spatial tibble
df_all <- cbind( r_all, sage_all_30m, sage_all_300m )
#check
head(df_all)
names(df_all)
# What proportion of our data are missing values
sum( is.na( sage_all_30m ))/ length( sage_all_30m )
sum( is.na( sage_all_300m ))/ length( sage_all_300m )
#simplify names of columns
names(df_all)[4:5] <- c( "sage_30m", "sage_300m" )
# Scale predictors 
df_all$sage_30m <- scale( df_all$sage_30m )
df_all$sage_300m <- scale( df_all$sage_300m )

#replace missing values with mean, which is 0 after they have been scaled
df_all$sage_30m[ is.na(df_all$sage_30m) ] <- 0
df_all$sage_300m[ is.na(df_all$sage_300m) ] <- 0
# we also assign weights to available points to be much greater than used points
df_all$weight <- 1000 ^( 1 - as.integer(df_all$case_ ) )
#check
head( df_all )

#alternatively we can choose to extract available from the estimated #
# ranges that we obtained for each individual
#extract individual id numbers:
idnos <- sort( unique( trks.thin$id )) 
#check your homerange data
akde_all
#now get random points
df_inds <- akde_all %>% 
mutate( 
  rsf_pnts =  map( hr_akde_all,
            ~ random_points(., n = nrow(.$data)*10, presence=.$data) ) )
#view
df_inds

#now unnest the new dataframes to make sure they worked
rsf_pnts <-  df_inds %>% 
  dplyr::select( id, rsf_pnts ) %>% 
  unnest( cols = rsf_pnts ) 
#check
head( rsf_pnts );dim(rsf_pnts)
plot( rsf_pnts )

#now extract predictor variables for these points###
## Answer: 
#

#######################################################################
##### RSF analyses #################

# Start with amt package  using the NCA as available habitat
mp1 <- df_all %>%  amt::fit_rsf( case_ ~ sage_30m ) %>% 
            summary()

#amt function doesn't take weights so we move to a more flexible #
# package for more accurate analyses
# We focus on comparing models between our two scales:
mp_30m <- glmmTMB( case_ ~ sage_30m,  
                   family = binomial(), data = df_all, 
                weights = weight ) 

summary( mp_30m )

# run analysis for the coarser scale:
mp_300m <- glmmTMB( case_ ~ sage_300m,  
                   family = binomial(), data = df_all, 
                   weights = weight ) 

summary( mp_300m )

#which scale has the most support? How do we choose?

# Interpreting results ###
# since we only have one predictor we don't have to account for
# others in the model and can just exponentiate it to compare the #
# relative intensity or rate of use of two locations that differ by 1 
# SD unit of the explanatory variable but are otherwise equivalent - i.e. #
# they are equally accessible and have identical values for all #
# other explanatory variables. # 
exp( glmmTMB::fixef( mp_30m )$cond[2] )
# this reflects the relative selection strength for choosing sagebrush 
#suggesting that prairie falcons are 1.24 times more likely to choose
# sagebrush with cover that is 1 SD higher. 
# if we want to remind ourselves what the SD for our predictor is
sd( sage_all_30m, na.rm = TRUE )

# we also plot differences in distribution between used and available #
# locations for our predictor of choice. To plot on the real scale we #
# combine unscaled data first:
cbind( cbind( r_all, sage_all_30m, sage_all_300m ) ) %>% 
ggplot( . ) +
  theme_bw( base_size = 15 ) +
  geom_density( aes( x = sage_all_30m, 
                     fill = case_, group = case_ ),
                alpha = 0.5  )

#### can you rerun analysis using the individual-derived points? #####
# Answer:
#

###########################################################
### Save desired results                                  #
# we can save the movement model results
#save workspace if in progress
save.image( 'RSFresults.RData' )
############# end of script  ##################################