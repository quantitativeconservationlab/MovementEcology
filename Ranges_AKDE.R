##################################################################
# Script developed by Jen Cruz to estimate ranges using AKDE     # 
# For this script we rely on Fleming et al.(2015) Ecology 96(5):1182-1188#
# We use ctmm first, and then use amt                           #
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
library(ctmm )#for more detailed functionality 
#####################################################################
## end of package load ###############

###################################################################
#### Load or create data -----------------------------------------

#load cleaned data:
#download full data for breeding season monitoring of prairie #
# falcons at the Birds of Prey NCA
trks.breed <- read_rds( "Data/trks.breed" )
#view
trks.breed
#download the thinned (30min) data
trks.thin <- read_rds( "Data/trks.thin" )
#view
trks.thin
#check class
class( trks.thin )
#check that the crs was correctly imported 
get_crs( trks.thin )

###############################################################
##### Estimate ranges using AKDE continuous-time movement model:#
################################################################

# We start by plotting points for each individual:
#
#you can choose which dataset by uncommenting and commenting lines 51 and 52 respectively
#trks.breed %>% 
trks.thin %>%   
  ggplot(., aes( x = x_, y = y_, color = speed ) ) +
  theme_bw( base_size = 15 ) + 
  geom_point() +
  facet_wrap( ~territory, scales = "free" )
# note that I color code speed...what else could you color code?

####################################################################
############## estimating ranges in CTMM ###########################
###################################################################

######### variograms using ctmm ###############
# We can use ctmm to explore the remaining autocorrelation in our data #
#Compare variograms for thinned and unthinned data for 1 individual 
#choose an individual id:
i <- 2
# filter tracks to select that individual's data
t <- trks.thin %>% filter( id == i )
a  <- trks.breed %>% filter( id == i )
#convert to ctmm object
ctmm.t <- as_telemetry( t )
ctmm.a <- as_telemetry( a )
#estimate empirical variograms
svf.t <- variogram( ctmm.t )
svf.a <- variogram( ctmm.a )
#now plot them side by side
par(mfrow = c(2,1) )
plot(svf.t, fraction = 1, level = 0.95)
plot(svf.a, fraction = 1, level = 0.95 )
#That second plot took a long time on my computer....
# for efficiency in the class I'm sticking to the thinned #
# data, moving forward. #

##### ALL individuals using ctmm     ###############
###########
#Plot variograms for all individuals
# extract names for individuals first into an object
ids <- unique( trks.thin$territory )
#create objects to store results
svf.t <- list()
ctmm.t <- list()
#xlimz <- c(0,36 %#% "hour" )
#set plot parameters
par( mfrow = c(3,3))
#loop through all individuals 
for( i in 1:length(ids)){
  #print progress
  print( i )
  # extract data for individual i
  t <- trks.thin %>% filter( id == i )
  #convert to ctmm object and add to list
  ctmm.t[[i]] <- as_telemetry( t )
  #Calculate empirical variograms:
  svf.t[[i]] <- variogram( ctmm.t[[i]] )
  #plot variograms for each individual
  plot( svf.t[[i]])#, xlim =  xlimz )
}
# Note how they are unique for each individual. This supports 
# the authors suggestions to estimate unique movement models
# for each individual

# automate the process of estimating a suitable movement     #
# model for the observed data using the empirical            # 
# variogram as a guide.                                     #
# options are "iid": for uncorrelated independent data,      #
#  "bm": Brownian motion, "ou": Ornstein-Uhlenbeck process,  #
# "ouf": Ornstein-Uhlenbeck forage process,                  #
# "auto": uses model selection with AICc to find bets model  #
# These model choices have real consequences to inference    #

#create and object to store results 
m.best <- list()
#loop through each individual
#this won't be fast...remember that we are estimating all #
# possible movement models for each individual and then #
# using AIC to pick a best model from the model choices #
# we also plot the empirical variograms vs the model results #
for( i in 1:2){#length(ids)){
  print( i )
  #use empirical variogram estimated in the previous step 
  # as a way of guiding the choice of movement model
  guess <- ctmm.guess(data = ctmm.t[[i]], variogram = svf.t[[i]],
                      interactive = FALSE )
  #here we actually compare among 6 movement model options 
  # and compare fit using AIC to select the top model
  m.best[[i]] <- ctmm.select( ctmm.t[[i]], guess, verbose = TRUE,
                              trace = 2 )
  #view summary output for model comparison for each individual
  print(summary( m.best[[i]] ))
}
#use individual names to replace those in the list:
names( m.best ) <- ids

#define plotting parameters:
par(mfrow = c(2,2))
#Now compare top model choice against traditional KDE
for( i in 1:2){#length(ids) ){
  #trace progress:
  print(i)
  # add basic IID model to model list
  m.best[[i]]$"IID isotropic" <- ctmm.fit( ctmm.t[[i]],
                                    ctmm(isotropic = TRUE) )
  #extract model name for top model
  an <- rownames(summary( m.best[[i]][1]))
  #plot best model
  ctmm::plot( svf.t[[i]], m.best[[i]][[1]], 
        main = paste( ids[i], an ) )#best model
  #plot worst
  ctmm::plot( svf.t[[i]], m.best[[i]]$"IID isotropic", 
        main = paste( ids[i], "IID" ) ) #worse model
}  
# Now that we have estimated top movement models for each #
# individual we are ready to apply those models to our estimates #
# of ranges. We also have an extra option to choose from #
# we can weight points based on high utilisation to correct the range #
# estimate...refer to the manuscript or the vignettes for more details #

# Here we compare ranges from 3 options: (1) top movement model 
# weighted (2) top movement model without weighing (3) traditional kde no weighing

# we create objects to store output from our 3 options:
akde.uw <- list()
akde.w <- list()
kde.iid <- list()
# We loop through each individual to estimate ranges for each option:
for( i in 1:2){#length(ids) ){
  print(i)
  # using the top movement model without weights
  akde.uw[[i]] <- ctmm::akde( ctmm.t[[i]], m.best[[i]][[1]] )
  # using the top movement model with weights
  akde.w[[i]] <- ctmm::akde( ctmm.t[[i]], m.best[[i]][[1]], 
                        weights = TRUE )
  #using the IID movement model without weights
  kde.iid[[i]] <- ctmm::akde( ctmm.t[[i]], m.best[[i]]$"IID isotropic" )
}

#plot estimate ranges comparing output for each option:
par(mfrow = c(3,2))
for( i in 1:2){#length(ids) ){
  print(i)
  plot( ctmm.t[[i]], akde.w[[i]] )
  title( paste("Weighted best model", ids[i]) )
  plot( ctmm.t[[i]], akde.uw[[i]] )
  title("Unweighted best model")
  plot( ctmm.t[[i]], kde.iid [[i]])
  title("Traditional KDE" )
}

##############################################################
# Estimating AKDE using atm package                          #
##############################################################
################
# We use amt package (talks to ctmm) to estimate AKDE for a #
# single individual to check for computational efficiency #

# define the individual id as an object at the start so you can 
# easily change it and try new ones, without having to alter the #
# rest of the code
i <- 2
#filter data
idv.thin <- trks.thin %>% filter( id == i )
#inspect details for the chosen animal
idv.thin

# start estimating traditional KDE model
f.t.iid <- amt::fit_ctmm( idv.thin, "iid" )
#check
summary( f.t.iid )

#now fit the model chosen as top model by ctmm
f.t.ou <- amt::fit_ctmm( idv.thin, "ou" )
#summary values
summary( f.t.ou )
#What is the run time for that individual?
# Answer:
#
#Estimate ranges for that individual using the KDE 
# approach that assumes no autocorrelation (traditional KDE)
akde_iid <- idv.thin %>% 
  amt::hr_akde(., model = f.t.iid, #fit_ctmm(., "iid" ),
               levels = 0.95 )

#run the ou model
akde_ou <- idv.thin %>% 
  amt::hr_akde(., model = fit_ctmm(., "ou" ),
               levels = 0.95 )

#Instead of manually choosing, we can do model selection for each 
# individual to find the optimal way of dealing with the observed
# movement behavior and autocorrelation #
#model selection
akde_auto <- idv.thin %>% 
  amt::hr_akde( ., model = fit_ctmm(., "auto" ),
                levels = 0.95 )

#Estimate areas for each method
amt::hr_area( akde_ou ) 
amt::hr_area( akde_iid ) 
amt::hr_area( akde_auto )
# comment on the results
# Answer:
#

#Plot comparisons from the different data choices
ggplot() +
  theme_bw( base_size = 15 ) +
  #extract isopleths for AKDE estimates from
  #model selection approach:
  geom_sf( data = hr_isopleths( akde_auto ),
           fill = "blue", col = "blue", size = 1 ) +
  #extract isopleths for ou model
  geom_sf( data = hr_isopleths( akde_ou ),
           fill = NA, col = "black", size = 3 ) +
  #extract isopleths for traditional kde:
  geom_sf( data = hr_isopleths( akde_iid ), 
           fill = NA, col = "grey", size = 2 ) +
  #add points of autocorrelated data
  # #note that we turn them into sf points for plotting
  # geom_sf( data = as_sf_points( subset( trks.breed, id == ind) ) ) +
  geom_sf( data = as_sf_points( subset( trks.thin, id == 2)), 
           col = "red" )


##### What if you want to use atm for all individuals? #######
##### don't try this in class. It will take a long time #
akde_all <- trks.thin %>% nest( data = -"id" ) %>% 
  mutate( hr_akde_all = map( data, ~hr_akde( ., 
              model = fit_ctmm(., "auto" ),
                levels = 0.95 ) ) )
# 

###########################################################
### Save desired results #
# we can save the movement model results
save( m.best,file="../ctmmresults.rda") # save where you want
#save range results for 3 home range options
save( akde.w,file="../ctmm_akde_w.rda")
save( akde.uw,file="../ctmm_akde_uw.rda")
save( kde.iid,file="../ctmm_akde_iid.rda")
#save range for individual estimated using atm
write_rds( akde_auto, "Data/akde_auto") 
#save range for all individuals in atm
write_rds( akde_all, "Data/akde_all" )
#save range for al individuals estimated with ctmm
write_rds( akde.w, "Data/weighted_akdes" )
#save workspace if in progress
save.image( 'AKDEresults.RData' )
############# end of script  ##################################
