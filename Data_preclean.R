#################################################################
# Script developed by Jen Cruz to clean and format location data #
# Data are Prairie Falcon locations collected during Spring/Summer #
# of 2021 at Morley Nelson Birds of Prey NCA.   

# Vegetation cover  was downloaded from Rangeland Analysis Platform #
# https://rangelands.app/products/ for 2021 and includes #
# % cover for shrub, perennial herbaceous, annual herbaceous #
# tree, litter and bare ground #
# coordinate system is WGS84 EPSG:4326, spatial resolution is 30m #
#################################################################

################## Prep. workspace ###############################
# we will be using new packages:
install.packages( "terra" )
install.packages( "raster" )
install.packages( "rasterVis" )

# load packages relevant to this script:
library( tidyverse ) #easy data manipulation and plotting
# set option to see all columns and more than 10 rows
options( dplyr.width = Inf, dplyr.print_min = 100 )
library( amt ) #creating tracks from location data
library( sf ) #handling spatial data
library( terra ) # for raster manipulation
library( raster )
library( rasterVis ) #for raster visualization (of raster:: objects)

## end of package load ###############
###################################################################
#### Load or create data -----------------------------------------
# Clean your workspace to reset your R environment. #
rm( list = ls() )

#import polygon of the NCA as sf spatial file:
NCA <- sf::st_read( "Data/BOPNCA_Boundary.shp") 

#load RAP cover
#note I don't have RAP data downloaded
#also it doesn't come labeled so below we stack it with the labels
cover <- raster::stack( "Data/vegetation-cover-v3-2021.tif" )

#name stack in the appropriate order
names(cover) <- c( "annual", "bare", "litter", "perennial", 
                 "shrub", "tree" )
#alternatively if you want to use the already cropped version:
cover_NCA <- raster::stack( "Data/RAPcover2021_NCA.img" )

#load 30m steps estimated for all individuals
trks.steps <- read_rds( "Data/trks.steps30" )

#load range for all individuals calculated in amt, which also
#contains tracks data used to estimate the range
akde_all <- read_rds( "Data/akde_all" )
akde_all
##############
#######################################################################
######## preparing raster and polygon data ###############################################
#get coordinates from shapefile
crstracks <- sf::st_crs( NCA )
# checking outline of NCA
sf::st_bbox( NCA )
# We define available habitat as area of NCA with a small buffer #
# around it and draw points from it #
#create a buffer around the NCA using outline of NCA and sf package:
NCA_buf <- NCA %>% sf::st_buffer( dist =1e4 )

#create a version that matches coordinates of the predictor raster:
NCA_trans <- sf::st_transform( NCA_buf, st_crs( cover ) ) 
# this step is crucial. We don't want to match the raster to the #
# polygon. Why?
# Answer:
# 
#create NCA polygon also in the same as raster crs for plotting
NCA_shape <- sf::st_transform( NCA, st_crs( cover ) ) 
#Need to turn sf objects to spatial objects:
NCA_poly <- as(st_geometry(NCA_shape ), Class="Spatial")
#compare outline of trasformed polygon:
sf::st_bbox( NCA_trans )
#check extent
terra::ext( NCA_trans )

# Not all vegetation cover are of interest so we remove 
# layers that we don't care about. This improves efficiency #

# List vegetation layers that you want to keep
keep <- c( "annual", "perennial", "shrub" )

#use list subset raster stack
cover <- cover[[keep ]]

# Now crop the rasters to your study area so that they are manageable
# computation wise. Notice we have to use the polygon that has been 
# transformed to match the CRS of the study area. 
cover_NCA <- raster::crop( cover, NCA_trans )

#select color palete for raster
mapTheme <- rasterTheme( region =  
                           RColorBrewer::brewer.pal(8,"Greens"))
#visualize
rasterVis::levelplot( cover_NCA , margin = FALSE,
                      par.settings = mapTheme )  +
  latticeExtra::layer(sp.lines( NCA_poly, col = "blue",
                                lwd = 3 ) )              

#or with terra:
terra::plot( cover_NCA )
#Note the range of values for each vegetation. 
#What are they?
# Answer: 
#
##############
###############################################################3
######### preparing tracks and step objects  #########################
#####################################################
#We will need to do habitat selection data seperately from step selection, 
#but intergreated step selection will be the same data as step selection
###### for RSFs #########

#For RSF we want to draw random points inside the home range of 
# each individual
#check homerange data
akde_all 

#define how many random points we want to draw per individual:
# as factor that total points will get multiplied by
rn <- 10
#now get random points
df_inds <- akde_all %>% 
  mutate( 
    rsf_pnts =  map( hr_akde_all,
     ~ random_points(., n = nrow(.$data)*rn, presence=.$data) ) )
#view
df_inds

#now unnest the new dataframes to make sure they worked
rsf_pnts <-  df_inds %>% 
  dplyr::select( id, rsf_pnts ) %>% 
  unnest( cols = rsf_pnts ) 
#check
head( rsf_pnts )

#unnest tracks also
trks.thin <- df_inds %>% 
  dplyr::select( id, data ) %>% 
  unnest( cols = data ) 
#view
head( trks.thin)
#check that the correct number of points was extracted
dim(rsf_pnts); dim(trks.thin)

#convert to sf object defining coordinate column
r_all_sf <- sf::st_as_sf( rsf_pnts, coords = c("x_", "y_"), 
                          crs = crstracks )
#now transform to predictor crs:
r_all_trans <- sf::st_transform( r_all_sf, st_crs(cover_NCA) )

#extract predictor at 30m resolution:
cover_tracks_30m <- raster::extract( x = cover_NCA, r_all_trans,
                                 method = "simple" )

#check
cover_tracks_30m
#add resolution to the column labels
colnames(cover_tracks_30m) <- paste( colnames(cover_tracks_30m),
                                     "30m", sep = "_" )
#Now at 300m following process for one individual above
cover_tracks_500m <- raster::extract( x = cover_NCA, r_all_trans,
                                  method = "simple", buffer = 250,
                                  fun  = mean, na.rm = TRUE )

#add resolution to the column labels
colnames(cover_tracks_500m) <- paste( colnames(cover_tracks_500m),
                                     "500m", sep = "_" )

# What proportion of our data are missing values
sum( is.na( cover_tracks_30m ))/ length( cover_tracks_30m )
sum( is.na( cover_tracks_500m ))/ length( cover_tracks_500m )

# We append our predictor estimates to the original spatial tibble
df_all <- cbind( rsf_pnts, cover_tracks_30m, cover_tracks_500m )
#check
head(df_all)

#Now we readd individual attributes that we want to keep
head( trks.thin )
#create a new df with extracted attributes
iddf <- trks.thin %>% 
  dplyr::select( id, territory, sex ) %>% 
  group_by( id ) %>% 
  slice(1)
#view
iddf

#join to your dataframe
df_all <- left_join( df_all, iddf, by = "id" )
#check
tail( df_all)

################
##################################################################
###### Now for SSF and iSSF ###
##################################################################
#######
# Step analysis is at a finer resolution so we don't need the ranges. #
# Instead we rely directly on the steps we created in cleaningtracks.R#
# Here we  use the 30min resolution for computational efficiency.

head(trks.steps)

# start by nesting data using purr to replicate what we did above:
steps_all <- trks.steps %>% nest( data = -"id" )
#view
steps_all

#Draw random steps for each individual
steps_all <- steps_all %>% 
  dplyr::mutate( rnd = lapply( data, function(x){
    amt::random_steps( x ) } ) )
#using the random steps function now

# The default number of random steps drawn per individual is 10.
# By default the random_steps() function fits a tentative #
# gamma distribution to the observed step lengths and a tentative #
# von Mises distribution to the observed turn angles. #
# It then generates random available points by sampling step-lengths #
# and turn angles from these fitted distributions and combining these #
# random steps with the starting locations associated with each observed #
# movement step. #

#now unnest the new dataframes to make sure they worked
stepsdf <- steps_all %>% dplyr::select( id, rnd ) %>% 
  unnest( cols = rnd ) 
#view
head( stepsdf );dim( stepsdf )
# because we are focused on habitat selection, we choose the end 
# of our steps when extracting habitat. 
# We start by turning it to sf object, assigning the correct projection
#Note that we use the second set of coordinates:
steps_sf <- sf::st_as_sf( stepsdf, coords = c("x2_", "y2_"), 
                          crs = crstracks )
#view
steps_sf

# We then transform the crs:
steps_trans <- sf::st_transform( steps_sf, st_crs(cover_NCA) )

#extracting with raster we can used the sf object directly, you also 
# have the choice to use a buffer around each point if you want to increase 
# your resolution:
cover_steps_30m <- raster::extract( x = cover_NCA, steps_trans,
                             method = "simple" )

#check
head(cover_steps_30m)
# What proportion of our data are missing values
sum( is.na( cover_steps_30m ))/ length( cover_steps_30m )

#add resolution to the column labels
colnames(cover_steps_30m) <- paste( colnames(cover_steps_30m),
                                      "30m", sep = "_" )

#extract at 100 m resolution
cover_steps_100m <- raster::extract( x = cover_NCA, 
                              steps_trans,
                              method = "simple", buffer = 50, 
                              fun  = mean, na.rm = TRUE )

#check
head( cover_steps_100m )
# What proportion of our data are missing values
sum( is.na( cover_steps_100m ))/ length( cover_steps_100m )
#add resolution to the column labels
colnames(cover_steps_100m) <- paste( colnames(cover_steps_100m),
                                    "100m", sep = "_" )
# We append our predictor estimates to the original spatial tibble
df_steps <- cbind( stepsdf, cover_steps_30m, cover_steps_100m )
#check
head( df_steps )

#Now we readd individual attributes that we want to keep
#Note that this assumes that the same individuals are #
#found in both your track and step objects:
#join to your dataframe
df_steps <- left_join( df_steps, iddf, by = "id" )
#check
tail( df_steps)
#########
###############################################################
###########    save relevant objects  ########

# If you created it yourself save cropped raster stack
# terra::writeRaster(  cover_NCA, "Data/RAPcover2021_NCA.img",
# overwrite=TRUE )

# Save the tracks and steps dataframe with extracted raster values 
#appended to them
write_rds( df_all, "Data/df_all" )

write_rds( df_steps, "Data/df_steps" )

#save workspace if in progress
save.image( 'HabSelDataCleanWorkspace.RData' )

################### end of script ##############################
