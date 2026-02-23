#################################################################
# Script developed by Jen Cruz to format cleaned GPS points   #
# to even sample rates and to add available points for      #
# habitat selection using 3 methods: RSF, SSF, iSSFs          # 
# We rely heavily on amt getting started vignette here:       #
# https://cran.r-project.org/web/packages/amt/vignettes/p1_getting_started.html#
#                                                               #
# Use locations of Prairie Falcon were collected during Spring/Summer #
# of 2021 at Morley Nelson Birds of Prey NCA.                      #
# Data were collected for multiple individuals and multiple #
# fix rates which we resampled at 5 secs and 30min fixes  #

# Predictors are vegetation cover metrics from Rangeland Analysis Platform #
# https://rangelands.app/products/  #
# you can go to that website to obtain vegetation cover or biomass #
# for western USA from satellite imagery processed at 30 x 30 m resolution #
# Habitat types for the 48 states can also be obtained here:#
# https://www.mrlc.gov/data     #

# we use tmap for spatial data visualization. Learn more at:
# https://r-tmap.github.io/tmap-book/index.html #

# Note that this code takes a lot of time to run #
#For homework choose one section that you need to adapt to your #
# class project, or for your own research, or the one that you #
# are most interested in learning. Submit modified code #
# detailing what you choose and why. Do not include other sections#
#################################################################

################## Prep. workspace ###############################
#install packages
install.packages( "tmap") 
install.packages( "corrplot" )
# load packages relevant to this script:
library( tidyverse ) #easy data manipulation and plotting
# set option to see all columns and more than 10 rows
options( dplyr.width = Inf, dplyr.print_min = 100 )
library( amt ) #handling tracks from location data
library( sf ) #handling vector data
library( raster ) #handle raster data
library( tmap ) #visualize raster and vector data together
library( corrplot )
## end of package load ###############

###################################################################
#### Load or create data -----------------------------------------
# Clean your workspace to reset your R environment. #
rm( list = ls() )

# If this is not the first time working on this script load workspace
# to pick up where you left off
load( "DataCleanRSFs.RData" )

#if you are starting new then load your data:
#import polygon of the NCA as sf spatial file:
NCA_Shape <- sf::st_read(  "Data/BOPNCA_Boundary.shp" )

#import nest locations
nests <- sf::st_read(  "Data/nests.shp" )

#We cropped the vegetation cover raster to our study area so that #
# the image could be shared via github, which has size restrictions #
# Load the cropped raster here:
cover_NCA <- raster::stack( "Data/RAPcover2021_NCA.img" )

#import akde ranges you created which includes  thinned (30min) data
akde_all <- read_rds( "Data/akde_all" )

#########
###################################################################
####  creating available points to match used data          ######
######################################################################

# our GPS trackers only give us used data but we now have to #
# derived available points for our habitat selection analysis #
# We will derive available data at (1) the study area level (1st-order #
# selection ), (2) the home range of each individual (2rd order ) #
#(3) specific for each used point (3rd order) #
# (1) and (2) scales can be analysed using RSFs, while (3) requires #
# SSFs or iSSF approaches #

# For scales (1) and (2) of inferences we do not need high resolution #
# data so we used our data resampled at 30 min intervals
trks.thin <- akde_all %>% 
  dplyr::select( id, data ) %>% 
  unnest( cols = data ) 

#define how many random points we want to draw per individual:
# as factor that total points will get multiplied by
rn <- 5

#For scale 1 We create random points inside study area
sa_pnts <- random_points( NCA_Shape, n = nrow( trks.thin )*rn,
                         type = "random", presence = trks.thin )

#view
head( sa_pnts )

#For scale 2 we create random points inside each individual home range
hr_pnts <- akde_all %>% 
  mutate( 
    rsf_pnts =  map( hr_akde_all,
    ~ random_points(., n = nrow(.$data)*rn, presence=.$data) ) ) %>% 
  #unnest including individual ids
  dplyr::select( id, rsf_pnts ) %>% 
  unnest( rsf_pnts )

#view
head( hr_pnts)

######################
### Extract values of vegetation cover #######################
##########
######## Prepare raster and polygon data ############################

#convert your objects with used and available points to sf objects
# defining the coordinate column to be the same as study area
# since in previous scripts we made them match
#start with study area scale
sa_sf <- sf::st_as_sf( sa_pnts, coords = c("x_", "y_"), 
                       crs = st_crs(NCA_Shape) )

#repeat with home range level scale
hr_sf <- sf::st_as_sf( hr_pnts, coords = c("x_", "y_"), 
                       crs = st_crs(NCA_Shape) )


#To extract raster data at used and available points, we need #
# to turn our vector crs to crs used by the raster #
# We ALWAYS TURN VECTOR PROJECTIONS TO RASTERS NOT THE OTHER WAY AROUND

#now for points at scale 1:
sa_trans  <- sf::st_transform( sa_sf, st_crs( cover_NCA ) ) 

#for points at scale 2:
hr_trans <- sf::st_transform( hr_sf, st_crs( cover_NCA ) ) 

###### Let's check that our raster is suitable for use ###
#view raster attributes
cover_NCA
#note that names of vegetation layers did not save so we add them:
names(cover_NCA) <- c( "annual", "perennial", "shrub" )

#we visualize vegetation rasters
#tm_shape let's you select the object you want to plot
tm_shape( cover_NCA ) +
  #how you want to plot it
  tm_raster( title = "cover (%)" ) + 
  #here you select the NCA polygon to overlay
  tm_shape( NCA_Shape ) +
  #choose to plot the outline in black
  tm_borders( lwd = 3, col = "black" )
#######
#################################################################
# Extract cover around each point. #
#We have to think about scale here again #
# We choose 3 scales 
##########
#the finest resolution is 30 x 30 m cells extracting value at the cell
sa_cover_30m <- raster::extract( x = cover_NCA, sa_trans,
                                 method = "simple" )
#repeat for home range selection:
hr_cover_30m <- raster::extract( x = cover_NCA, hr_trans,
                                 method = "simple"  )

#check
head( hr_cover_30m )
#add resolution to the column labels
colnames(sa_cover_30m) <- colnames(hr_cover_30m) <- paste( colnames(sa_cover_30m),
                                     "30m", sep = "_" )
#check
head( sa_cover_30m );head( hr_cover_30m)

#Now scale to 100m by using a buffer of 50m radius
sa_cover_100m <- raster::extract( x = cover_NCA, sa_trans,
                                  method = "simple", buffer = 50,
                                  fun  = mean, na.rm = TRUE )
#repeat for home range selection:
hr_cover_100m <- raster::extract( x = cover_NCA, hr_trans,
                                  method = "simple", buffer = 50,
                                  fun  = mean, na.rm = TRUE )

#add resolution to the column labels
colnames(sa_cover_100m) <- colnames(hr_cover_100m) <- paste( colnames(sa_cover_100m),
                                "100m", sep = "_" )

#Now scale to 500m by using a buffer of 250m radius
sa_cover_500m <- raster::extract( x = cover_NCA, sa_trans,
                                  method = "simple", buffer = 250,
                                  fun  = mean, na.rm = TRUE )
#repeat for home range selection:
hr_cover_500m <- raster::extract( x = cover_NCA, hr_trans,
                                  method = "simple", buffer = 250,
                                  fun  = mean, na.rm = TRUE )
#add resolution to the column labels
colnames(sa_cover_500m) <- colnames(hr_cover_500m) <- paste( colnames(sa_cover_500m),
                                      "500m", sep = "_" )

# What proportion of our data are missing values
sum( is.na( hr_cover_30m ))/ length( hr_cover_30m )
sum( is.na( hr_cover_100m ))/ length( hr_cover_100m )
sum( is.na( hr_cover_500m ))/ length( hr_cover_500m )

# for 1st-order selection (study area) we combine our 3 scales
df_sa <- cbind( sa_cover_30m, sa_cover_100m, sa_cover_500m ) 
#view
head(df_sa)
# we merge with our gps points from ORIGINAL crs (not the cover raster #
# since we already extracted cover values) 
head( sa_pnts)
dim( df_sa ); dim(sa_pnts)
df_sa <- cbind( sa_pnts, df_sa )
head( df_sa )

# for 2st-order selection (home range) we combine our 3 scales
df_hr <- cbind( hr_cover_30m, hr_cover_100m, hr_cover_500m ) 
#view
head( df_hr )
# we merge with our gps points from ORIGINAL crs (not the cover raster #
# since we already extracted cover values) 
head( hr_pnts)
dim( df_hr ); dim(hr_pnts)
df_hr <- cbind( hr_pnts, df_hr )

# we want to add other individual attributes
#start by extracting them
head( trks.thin )
id_df <- trks.thin %>% 
  dplyr::select( id, territory, sex ) %>% 
  dplyr::group_by( id ) %>% 
  slice(1)
#view
id_df
#combine with our dataframe
df_hr <- left_join( df_hr,id_df, by = "id" )
#view
head( df_hr )

# next we nee to check correlation among predictors 
#create vector of predictors taking advantage of naming commonality 
# to automatically extract them:
prednames <- grep('0m', colnames(df_hr), value = TRUE)
#check for correlation but also whether the values sum to one or close 
# to. See: https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.4256
# for reasons why that is an issue.
corrplot::corrplot( round(cor(df_hr[prednames]),1 ), method = "number" )

y <- rbinom( dim(df_hr)[1], 1, 0.4 )

prednames
#fit model
mod30 <- glm( y ~ annual_30m + perennial_30m + shrub_30m, 
              data = cbind(y,df_hr), family = binomial )
#calculate VIF
sort( car::vif(mod30), decreasing = T ) 
#fit model
mod100 <- glm( y ~ annual_100m + perennial_100m + shrub_100m, 
               data = df_hr, family = binomial)
#calculate VIF
sort( car::vif(mod100), decreasing = T ) 
#fit model
mod500 <- glm( y ~ annual_500m + perennial_500m + shrub_500m, 
               data = df_hr, family = binomial)
#calculate VIF
sort( car::vif(mod500), decreasing = T ) 

#what if we included all scales into the same model
modall <-  glm( y ~ annual_30m + perennial_30m + shrub_30m + 
                  annual_100m + perennial_100m + shrub_100m +
                  annual_500m + perennial_500m + shrub_500m, 
                data = df_hr, family = binomial)
#calculate VIF
sort( car::vif( modall ), decreasing = T ) 

#check for correlation but also whether the values sum to one or close 
# to. See: https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.4256
# for reasons why that is an issue.
#sum covariates at appropriate scale and then plot results:
#for 30m
hist(apply( df_hr[ ,prednames[1:3] ], 1,sum ))
#for 100m
hist(apply( df_hr[ ,prednames[4:6] ], 1,sum ))
#for 500m
hist(apply( df_hr[ ,prednames[7:9] ], 1,sum ))


#############################################################
######### step lengths and turning angles  ##################
#########################

#For scale (3) we want to concentrate on foraging/traveling and remove #
# nest locations and territorial movements #

#create reference dataframe that keeps individual information
iddf <- trks.thin %>% 
  group_by( id ) %>% 
  dplyr::select( id, territory, sex ) %>% 
  slice(1)

iddf
#create a buffer around the nest based on territory estimates 
nest_buffer <- nests %>% 
  dplyr::select( territory= terrtry ) %>% 
  st_buffer(750) %>% 
  right_join( iddf, by = "territory" )

#extract territory IDs
terids <- unique( nest_buffer$territory )

#add row id
trks.thin <- trks.thin %>% 
  dplyr::mutate( rowid = row_number())

#create a points only sf
thin_sf <- trks.thin %>% 
  dplyr::select( territory, rowid, x_, y_ ) %>% 
  amt::as_sf_points()
# calculation needs to be individual specific start with 1 indv
#choose buffer of individual
buf <-  nest_buffer %>% dplyr::filter( territory == "SG" )
#choose points belonging to that individual
bpnts <- thin_sf %>% dplyr::filter( territory == "SG" )
#get points outside the buffer using the st_difference()
b <- st_difference( bpnts, buf )
#create new dataframe to store results
forage.thin <- b
#now loop through the rest of the individuals
for( i in terids[2:length(terids)] ){
  buf <-  nest_buffer %>% dplyr::filter( territory == i )
  bpnts <- thin_sf %>% dplyr::filter( territory == i )
  #this function keeps nonoverlapping points only
  b <- st_difference( bpnts, buf )
  #append updated points to original
  forage.thin <- bind_rows(forage.thin, b )
}

#now filter our original thin dataset using rowids
trks.forage <- trks.thin %>% 
  dplyr::filter( rowid %in% forage.thin$rowid )
#the burst id restarts for each individual
# we need to create indiv specific ones
#add counts 
trks.forage <- trks.forage %>% 
  mutate( burstid = paste(burst_, id, sep = "_")) %>% 
  add_count( burstid )

#check 
head(trks.forage)
#how many points per burst?
table(trks.forage$n)
#note that we have a lot of bursts with a single point
#we remove those
trks.forage <- trks.forage %>% 
  group_by( burstid ) %>% 
  dplyr::filter( n > 1 ) %>% ungroup()

###put them all in the same map to check if it worked
ggplot() +
  theme_bw( base_size = 15 ) + 
  theme( legend.position = "bottom" ) +
  geom_sf( data = nest_buffer, 
           aes( color = as.factor(id) ),
           linewidth = 1.2 ) +
  geom_sf( data = as_sf_points( trks.forage ), 
           aes(color = as.factor(id)), size = 0.5 ) +
  geom_sf(data = NCA_Shape, linewidth = 1.5,
          inherit.aes = FALSE, fill=NA ) 

# We now connect the points into steps, which allow us to derive step lengths and turning angles #
# which are movement parameters that we can incorporate into iSSFs #
# select the data tibble and calculate movement metrics
# The default number of random steps drawn per individual is 10.
# By default the random_steps() function fits a tentative #
# gamma distribution to the observed step lengths and a tentative #
# von Mises distribution to the observed turn angles. #
# It then generates random available points by sampling step-lengths #
# and turn angles from these fitted distributions and combining these #
# random steps with the starting locations associated with each observed #
# movement step. #

# turn tracks into steps
steps_30 <- trks.forage %>% 
      steps_by_burst( keep_cols = 'start' ) 
#view
head(steps_30)

# We can plot step lengths by:
steps_30 %>%   
  ggplot(.) +
  #geom_density( aes( x = sl_, fill = as.factor(burst_)), alpha = 0.4 ) +
  geom_histogram( aes( x = sl_ ) ) +
  xlab("Step length" ) + 
  #ylim( 0, 0.01 ) + xlim(0, 2000 ) +
  theme_bw( base_size = 19 )  +
  theme( legend.position = "none" ) +
  facet_wrap( ~id, scales = 'free' )

#Note there are 0 step lengths we need to remove

# Turning angles:
steps_30 %>%
#removes individual 4 so we can see better for others with less data:
 # dplyr::filter( id != 4 ) %>% 
  ggplot(.) +
  geom_histogram( aes( x = ta_ ) ) +
  coord_polar() +
  ylab("Turning angle") + xlab("") + 
  theme_bw( base_size = 19 ) +
  facet_wrap( ~id) #, scales = 'free_y' )

#remove 0 step lengths and missing turning angles
steps_30 <- steps_30 %>% 
  dplyr::filter( !is.na(ta_) ) %>% 
  dplyr::filter( sl_ > 0 )

#Draw random steps for each individual
steps_30df <- steps_30 %>% amt::nest( data = -"id" ) 
#view
steps_30df
#add random steps
steps_30df <- steps_30df %>% 
  dplyr::mutate( rnd = lapply( data, function(x){
    amt::random_steps( x ) } ) ) 
#view
steps_30df
#unnnest
steps_30df <- steps_30df %>% 
    dplyr::select( id, rnd ) %>% 
    unnest( cols = rnd ) 
#view
head( steps_30df )

#now we can turn to sf object, transform CRS to CRS of cover raster #
#extract cover at 30 m due to fine resolution of analysis
# because we are focused on habitat selection, we choose the end 
# of our steps when extracting habitat. 
# We start by turning it to sf object, assigning the correct projection
#Note that we use the second set of coordinates:
steps30_sf <- sf::st_as_sf( steps_30df, coords = c("x2_", "y2_"), 
                          crs = st_crs(NCA_Shape) )
#view
head(steps30_sf)
# We then transform the crs:
steps30_trans <- sf::st_transform( steps30_sf, st_crs(cover_NCA) )
#view
head(steps30_trans )
#extracting with raster
cover_steps30 <- raster::extract( x = cover_NCA, steps30_trans,
                                    method = "simple" )

head(cover_steps30)
#now combine cover values with original step dataframe
cover_steps30_df <- cbind( steps_30df, cover_steps30 )
#check
head( cover_steps30_df)
#dim( cover_steps30_df )
#colnames( cover_steps30_df)[24:26] <- c( "annual", "perennial", "shrub" )

#############################################################################
# Saving relevant objects and data ---------------------------------
write.csv( df_hr, "Data/df_hr.csv" )
write.csv( df_sa, "Data/df_sa.csv" )
write.csv( cover_steps30_df, "Data/df_steps30.csv" )

#save workspace if in progress
save.image( 'DataCleanRSFs.RData' )


########## end of save #########################
############### END OF SCRIPT ########################################