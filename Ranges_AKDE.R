##################################################################
# Script developed by Jen Cruz to estimate ranges using AKDE     # 
# For this script we rely on Fleming et al.(2015) Ecology 96(5):1182-1188#
# We aim to estimate AKDE using high-resolution data and compare #
# with results from data that has been thinned to remove        #
# autocorrelation                                               #
# For instructions on how to use ctmm directly check out:       #
# https://cran.r-project.org/web/packages/ctmm/vignettes/variogram.html #
# https://cran.r-project.org/web/packages/ctmm/vignettes/akde.html #
###################################################################

################## prep workspace ###############################

# Clean your workspace to reset your R environment. #
rm( list = ls() )
#install.packages( "ctmm" )

# load packages relevant to this script:
library( tidyverse ) #easy data manipulation
library( amt )
library( sf )
#library(ctmm )#for more detailed functionality 
#####################################################################
## end of package load ###############

###################################################################
#### Load or create data -----------------------------------------

#load cleaned data:
#download full data for breeding season monitoring of prairie #
# falcons at the Birds of Prey NCA
trks.breed <- read_rds( "trks.breed" )
#view
trks.breed
#download the thinned (30min) data
trks.thin <- read_rds( "trks.thin" )
#view
trks.thin
#check class
class( trks.thin )
#check that the crs was correctly imported 
get_crs( trks.thin )

###############################################################
##### Estimate ranges using AKDE continuous-time movement model:#
################################################################

# We start by plotting points from both datasets for each individual:
trks.breed %>% # filter( id == 2 ) %>%
#trks.thin %>%   
  ggplot(., aes( x = x_, y = y_, color = speed ) ) +
  theme_bw( base_size = 15 ) + 
  geom_point() +
  # geom_point( data = trks.thin, 
  #   aes( x = x_, y = y_, color = speed ), 
    # size = 3, shape = 8 ) +
  #labs( title = ids[i], fill = "week", x = "lat") +
  facet_wrap( ~territory, scales = "free" )

######### variograms using ctmm ###############
# We can use ctmm to explore the autocorrelation in our data #
# extract ids
ids <- unique( trks.thin$territory )
#loop through each animal
for( i in 1:2){#length(ids) ){ 
  #extract data for each individual for thinned and autocorrelated
  # datasets
  t <- trks.thin %>% filter( id == i )
  a  <- trks.breed %>% filter( id == i )
  #convert to ctmm object
  ctmm.t <- as_telemetry( t )
  ctmm.a <- as_telemetry( a )
  svf.a <- variogram( ctmm.a )
  plot(svf.a, fraction = 1, level = 0.95)
}  
#compare to what it looks like for unthinned irregularly
# sampled data for one individual
#define which individual you want to check 
i <- 2
# filter tracks
a  <- trks.breed %>% filter( id == i )
#convert to ctmm object
ctmm.t <- as_telemetry( t )
ctmm.a <- as_telemetry( a )
svf.a <- variogram( ctmm.a )
plot(svf.a, fraction = 1, level = 0.95)

#extract details for each animal
t <- trks.thin %>% filter( id == i )
#convert to ctmm object
ctmm.t <- as_telemetry( t )
#Use ctmm directly to calculate empirical variograms:
svf.t <- variogram( ctmm.t )
par(mfrow = c(2,1))
#plot it
plot(svf.t, fraction = 1, level = 0.95, main = ids[i])
# automate the process of estimating a suitable movement 
# model for the observed data
m.best <- ctmm.guess( ctmm.t, interactive = FALSE)
# estimate the IID model which assumes no autocorrelation
# in movement
m.iid <- ctmm.fit( ctmm.t )
# now fit the range estimate 
# using the top movement model without weights
akde.uw <- ctmm::akde( ctmm.t, m.best )
# using the top movement model without weights
akde.w <- ctmm::akde( ctmm.t, m.best, weights = TRUE )
#using the IID movement model without weights
kde.iid <- ctmm::akde( ctmm.t, m.iid )
#look at results


#inspect details for the chosen animal
idv.breed
idv.thin
# # you can look at it directly on leaflet: 
inspect( idv.breed )


##############################################################
# Now that we are ready to estimate AKDE, but what movement  #
# model option do we choose?                                 #
# options are "iid": for uncorrelated independent data,      #
#  "bm": Brownian motion, "ou": Ornstein-Uhlenbeck process,  #
# "ouf": Ornstein-Uhlenbeck forage process,                  #
# "auto": uses model selection with AICc to find bets model  #
# These model choices have real consequences to inference    #
##############################################################

# We use amt package (talks to ctmm) to estimate AKDE #
# We start with a single individual to check for computational #
# efficiency#

# define the individual id as an object at the start so you can 
# easily change it and try new ones, without having to alter the #
# rest of the code
ind <- 2
# subset the two datasets accordingly
idv.breed <- trks.breed %>% filter( id == ind )
idv.thin <- trks.thin %>% filter( id == ind )
ctmm.thin <- as_telemetry( idv.thin )
class( ctmm.thin )

#inspect details for the chosen animal
idv.breed
idv.thin
# # you can look at it directly on leaflet: 
inspect( idv.breed )

#Use ctmm directly to calculate empirical variograms:
SVF <- variogram( ctmm.thin ) 
par(mfrow = c(1,1))
plot(SVF, fraction = 1, level = 0.95)

# Calculate an automated model guesstimate:
GUESS1 <- ctmm.guess(animal1_buffalo, interactive = FALSE)

#fit models to data first
#start with thinned dataset and use traditional kde that assumes #
# no autocorrelation in the data:
f.t.iid <- fit_ctmm( idv.thin, "iid" )
#check
summary( f.t.iid )

#now fit the most complicated model version
f.t.ouf <- fit_ctmm( idv.thin, "ouf" )
#summary values
summary( f.t.ouf )
#What is the run time for that individual?
# Answer:
#


#Estimate home range for that individual using the KDE 
# approach that assumes no autocorrelation (traditional KDE)
akde_iid <- idv.thin %>% 
  amt::hr_akde(., model = f.t.iid, #fit_ctmm(., "iid" ),
               levels = 0.95 )

# # if it didn't take too long for one then you can try that 
# # approach for all
# akde_iid_all <- trks.thin %>% dplyr::filter( id == ind ) %>% 
#   amt::hr_akde(., model = fit_ctmm(., "iid" ),
#                levels = 0.95 )


#run the ouf model
akde_ouf <- idv.thin %>% 
  amt::hr_akde(., model = fit_ctmm(., "ouf" ),
               levels = 0.95 )

#Instead of manually choosing, we can do model selection for each 
# individual to find the optimal way of dealing with the observed
# movement behavior and autocorrelation #
#model selection
hr_akde <- idv.thin %>% 
  amt::hr_akde( ., model = fit_ctmm(., "auto" ),
                levels = 0.95 )


#Plot comparisons from the different data choices
ggplot() +
  theme_bw( base_size = 15 ) +
  #extract isopleths for autocorrelated AKDE estimates #
  # as we expect those would be the largest
  geom_sf( data = hr_isopleths( akde_ouf ),
           fill = NA, col = "black", size = 3 ) +
  #extract isopleths for AKDE estimates using thinned data
  #estimated with the model selection approach:
  geom_sf( data = hr_isopleths( akde_idd ), 
           fill = NA, col = "grey", size = 2 ) +
  #extract isopleths for AKDE estimates using single movement #
  # behavior method
  # geom_sf( data = hr_isopleths( akde_ou ), 
  #          fill = NA, col = "blue", size = 1 ) +
  #add points of autocorrelated data
  # #note that we turn them into sf points for plotting
  # geom_sf( data = as_sf_points( subset( trks.breed, id == ind) ) ) +
  #add thinned points
  geom_sf( data = as_sf_points( subset( trks.thin, id == ind)), 
           col = "red" )

# If computation time didn't hinder the process with thinned #
# data then we can try it for the full autocorrelated dataset #
# Here we combine the AKDE and model fit estimation in one step #

#start with the KDE traditional approach:
akde_all <- trks %>%
  mutate(
  akde_iid = map( data, ~hr_akde(., 
            model = fit_ctmm(., "iid" ),
               levels = 0.95 ) ) )

akde_idd_all <- amt::hr_akde( idv.breed, 
          model = fit_ctmm( idv.breed, "iid" ),
                levels = 0.95 )

#now for the more complext model
akde_ouf_all <- idv.breed %>% 
  amt::hr_akde(., model = fit_ctmm(., "ouf" ),
               levels = 0.95 )


#Estimate areas for each method
amt::hr_area( akde_ouf ) 
amt::hr_area( akde_ouf_all ) 
amt::hr_area( hr_akde )
# comment on the results
# Answer:
#

###########################################################
### Save desired results                                  #
#save breeding season data (not thinned)
write_rds( hr_wk, "hr_wk")
#save range area estimates
write_rds( ci_wk, "ci_wk" )
save.image( 'AKDEresults.RData' )
############# end of script  ##################################
