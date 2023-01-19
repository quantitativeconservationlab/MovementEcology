##################################################################
# Script developed by Jen Cruz to estimate SSFs and iSSFs          #
# approach derived from Fieberg et al. 2021 DOI: 10.1111/1365-2656.13441 #
# using code from Appendices B and C                             #
# also Muff et al. 2019 DOI: 10.1111/1365-2656.13087             #
# code here:                                                     #  
# https://conservancy.umn.edu/handle/11299/204737                #
#                                                                #
# We use landcover data from the National Geospatial Data Asset  #
# https://www.mrlc.gov                                           #
# Habitat predictors include 2018 estimates of sagebrush cover   #
###################################################################

################## prep workspace ###############################

# Clean your workspace to reset your R environment. #
rm( list = ls() )

#install new packages
install.packages(c("leaflet", "mapview"))
install.packages("moveVis")
 
# load packages relevant to this script:
library( tidyverse ) #easy data manipulation
# set option to see all columns and more than 10 rows
options( dplyr.width = Inf, dplyr.print_min = 100 )
library( moveVis )
#library( amt )
# library( sf )
# library( terra ) # for raster manipulation
# library( raster )
# library( rasterVis ) #for raster visualization (of raster:: objects)
# library( glmmTMB ) # for analysis

#####################################################################
## end of package load ###############

###################################################################
#### Load or create data -----------------------------------------

#load 30m steps estimated for all individuals
#trks.steps <- read_rds( "trks_all" )
# load points also so that we can combine data
trks <- read_rds( "trks.breed" )
#view
head( trks )

#import polygon of the NCA as sf spatial file:
NCA_Shape <- sf::st_read("Z:/Common/QCLData/Habitat/NCA/GIS_NCA_IDARNGpgsSampling/BOPNCA_Boundary.shp")

################## manipulating data ####################################
# we prepare the predictor data similarly to how we did it for RSFs:
#get coordinates from shapefile
crstracks <- sf::st_crs( NCA_Shape )

#add week and hour columns 
trks <- trks %>% mutate( wk = lubridate::week( t_ ), 
        hr =  lubridate::hour( t_ ) )  
#convert to dataframe
trksdf <- as.data.frame( trks )
#view
head( trksdf )
# now convert to move object
df <- df2move( trksdf, proj = "+proj=utm +zone=11",
               x = "x_", y = "y_", time = "t_",
               track_id = "territory" )

head(df)
# align move_data to a uniform time scale
m <- align_move( df, res = 5, unit = "mins" )

# create spatial frames with a OpenStreetMap watercolour map
frames <- frames_spatial( m, #path_colours = c("red", "green", "blue"),
    map_service = "osm", map_type = "watercolor", alpha = 0.5) %>%
  # add_labels(x = "Longitude", y = "Latitude") %>% # add some customizations, such as axis labels
   add_northarrow() %>%
   add_scalebar() %>%
  # add_timestamps(m, type = "label") %>%
  add_progress()

# return an interactive mapview map
view_spatial( df )

# animate frames
animate_frames(frames, out_file = "moveVis.gif" )
##########################################################################
### Save desired results                                  #
# I save the steps dataframe with extracted raster values so that I 
# don't have to recreate it when estimating issfs 
write_rds( trks_all, "trks_all" )
# I also save the unscaled raster values as a csv
write.csv( sage_30m, "sage_30m_steps.csv", row.names = FALSE )

#save workspace if in progress
save.image( 'SSF_results.RData' )
############# end of script  ##################################