##################################################################
# Script developed by Jen Cruz to estimate resource selection     #
# functions. Code adapted from atm vignette:                      #
#https://cran.r-project.org/web/packages/amt/vignettes/p3_rsf.html #
# but modified to remove the need to convert the predictor rasters #
# to a different CRS. Rasters were imported and manipulated using #
# the terra package, while sf was used for spatial objects.       #
# We use landcover data from the National Geospatial Data Asset  #
# https://www.mrlc.gov                                           #
# Habitat predictors include 2018 estimates of sagebrush cover  #
###################################################################

################## prep workspace ###############################

# Clean your workspace to reset your R environment. #
rm( list = ls() )

# load packages relevant to this script:
library( tidyverse ) #easy data manipulation
library( amt )
library( sf )
library( terra ) # for raster manipulation
#####################################################################
## end of package load ###############

###################################################################
#### Load or create data -----------------------------------------

#load the thinned (30min) data for all individuals
trks.thin <- read_rds( "trks.thin" )
#load range for a single individual:
akde_atm <- read_rds( "akde_auto" )
#w_akdes <- read_rds( "weighted_akdes" )
#load range for individuals estimated using atm
load( "../ctmm_akde_w.rda" ) 
#check that it loaded the object
class( akde.w )

#import polygon of the NCA as sf spatial file:
NCA_Shape <- sf::st_read("Z:/Common/QCLData/Habitat/NCA/GIS_NCA_IDARNGpgsSampling/BOPNCA_Boundary.shp")

#import sagebrush raster:
sagebrush <- terra::rast( "Z:/Common/QCLData/Habitat/NLCD/Sage_2007_2018/nlcd_sage_2018_mos_rec_v1.img" )
#view
sagebrush
#plot
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
NCA_buf <- NCA_Shape %>% sf::st_buffer( dist =5e3 )
#create a version that matches coordinates of the predictor raster:
NCA_trans <- sf::st_transform( NCA_buf, st_crs( sagebrush ) ) 
#compare outline of trasformed polygon:
sf::st_bbox( NCA_trans )
#check extent
terra::ext( NCA_trans )
#crop raster to buffered NCA
sage_cropped <- terra::crop( sagebrush, ext( NCA_trans ) )
# Now that we have cropped it to the appropriate area it should be faster #
# to process #

# compare extent with original
terra::ext(sage_cropped)
terra::ext( sagebrush )
#view
sage_cropped
#values greater than 100 are empty so replace with missing
sage_cropped[ sage_cropped >100 ] <- NA
#Plot sagebrush
terra::plot( sage_cropped, main = "Sagebrush (2018)" )
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
# We start with a single individual:
#extract data for individual of interest:
idv.thin <- trks.thin %>% filter( id == 2 )
# Extract available points within range of individual and #
# add used points from track:
r_one <- amt::random_points( akde_atm, presence = idv.thin )
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
r_one_sf <- sf::st_as_sf( r_one, coords = c("x_", "y_"), crs = crstracks )
#now transform to predictor crs:
r_one_trans <- sf::st_transform( r_one_sf, st_crs(sage_cropped) )
#convert to dataframe for use in terra
df_one <- r_one_trans %>%
  dplyr::mutate( x = sf::st_coordinates(.)[,1],
                 y = sf::st_coordinates(.)[,2]) %>% 
  st_drop_geometry() %>% 
  as.data.frame()
#view
head( df_one ); dim( df_one )
#now we can extract values from our predictor rasters while keeping them in #
# their original crs #
# We choose two spatial scales:
#1) original resolution of raster of 30 m grid cells
sage_30m <- terra::extract( x = sage_cropped, y = df_one[ ,c("x","y") ], 
                            method = "simple" )
#check
head( sage_30m ); dim( sage_30m )

#2) doubling it to 60 m
sage_60m <- terra::extract( x = sage_cropped, y = r_one_trans, 
                            method = "bilinear")

# Population-level data management
# We want to determine use within the NCA assuming 10 individuals #
# is a representative sample. When would this be the case? #
# When would it not be the case? #

#use buffer to define available area and the tracks for used points:
# we specify how many available points we want
r_all <- random_points( NCA_buf, n = (dim(trks.thin)[1] * 10 ), 
                        presence = trks.thin )
r_all
plot( r_all )
#extract predictor from used and available points for the population:
rsf_all <- r_all %>% amt::extract_covariates( sage_NCA )

########################################################
##### analyse data  ##########
#We can use fit_rsf, which is just a wrapper around 
#stats::glm with family = binomial(link = "logit").

# starting with single individual:
m1 <- rsf1 %>% fit_rsf( case_ ~ sage_NCA ) %>% 
  summary()

# What are results saying?
# Answer: 
#
# Analyse population habitat use:
m_all <- rsf_all %>%  atm::fit_rsf( case_ ~ sage_NCA ) %>% 
            summary()

# Interpret results:
# Answer:
# 
# How can we define available habitat differently? 
# Answer:
#
# What about standardising data?
# Answer:
#


###########################################################
### Save desired results                                  #
# we can save the movement model results
#save workspace if in progress
save.image( 'RSFresults.RData' )
############# end of script  ##################################