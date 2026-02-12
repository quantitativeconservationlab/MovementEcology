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
#uncomment and install the package if you haven't got it
#install.packages( "ctmm" )

# load packages relevant to this script:
library( tidyverse ) #easy data manipulation
# set option to see all columns and more than 10 rows
options( dplyr.width = Inf, dplyr.print_min = 100 )
library( amt )
library( sf )
library(ctmm ) #for more detailed functionality 
#####################################################################
## end of package load ###############

###################################################################
#### Load or create data -----------------------------------------

#load workspace if you have already started working through this script#
load( "AKDEresults.RData")

#if you are starting from scratch load cleaned data:
#download the thinned (30min) data
trks.thin <- read_rds( "Data/trks.thin" )
#view
head( trks.thin )
#download breeding ranges we estimated last week
ranges <- read_rds( "Data/ranges" )
#view import
head( ranges )

#import polygon of the NCA as sf spatial file:
NCA_Shape <- sf::st_read( "Data/BOPNCA_Boundary.shp" )

####################################################################
############## estimating ranges in CTMM ###########################
###################################################################

######### variograms using ctmm ###############
# We use ctmm to explore autocorrelation in our data #
# but using estimates of semi-variance instead of the acf() from 
# last week.

##### We start with a single individual  ######
#choose an individual id:
i <- 2
# filter tracks to select that individual's data
t <- trks.thin %>% filter( id == i )
#convert to ctmm object
ctmm.t <- as_telemetry( t )
#estimate empirical variograms
svf.t <- variogram( ctmm.t )
#now plot them side by side
par(mfrow = c(2,1) )
plot(svf.t, fraction = 1, level = 0.95)
#now zoom in to starting time lags
plot( svf.t, xlim = c(0,2 %#% "day"), 
     fraction = 1, level = 0.95 )

################
##################################################################
##### ALL individuals using ctmm     #############################
###################################################################
#Plot variograms for all individuals
# extract names for individuals first into an object
ids <- sort(unique( trks.thin$id ))
#create objects to store results
svf.t <- list()
ctmm.t <- list()
xlimz <- c(0,36 %#% "hour" )
#set plot parameters
par( mfrow = c(3,3))
#loop through all individuals 
for( i in ids ){
  #print progress
  print( i )
  # extract data for individual i
  t <- trks.thin %>% filter( id == i )
  #convert to ctmm object and add to list
  ctmm.t[[i]] <- as_telemetry( t )
  #Calculate empirical variograms:
  svf.t[[i]] <- variogram( ctmm.t[[i]] )
  #plot variograms for each individual
  plot( svf.t[[i]], xlim =  xlimz )
}
# How are they unique for each individual?
# Answer:
#

#########
##############################################################
# automate the process of estimating a suitable movement     #
# model for the observed data using the empirical            # 
# variogram as a guide.                                     #
# options are "iid": for uncorrelated independent data,      #
#  "ou": Ornstein-Uhlenbeck process,  #
# "ouf": Ornstein-Uhlenbeck forage process,                  #
# "auto": uses model selection with AICc to find bets model  #
################################################################

### we try the model selection method for our class example ###
### YOU DO NOT NEED TO DO THIS SECTION FOR HOMEWORK ######

#create and object to store results 
m.best <- list()
#loop through each individual
#this won't be fast...remember that we are estimating all #
# possible movement models for each individual and then #
# using AIC to pick a best model from the model choices #
# we also plot the empirical variograms vs the model results #
for( i in 1:length(ids)){
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
names( m.best ) <- ids#[1:2]

#define plotting parameters:
par(mfrow = c(2,2))
#Now compare top model choice against traditional KDE
for( i in 1:length(ids) ){
  #trace progress:
  print(i)
  # add basic IID model to model list
  m.best[[i]]$"IID isotropic" <- ctmm.fit( ctmm.t[[i]],
                                     ctmm(isotropic = TRUE) )
  #extract model name for top model
  an <- rownames(summary( m.best[[i]][1]))
  #plot best model
  ctmm::plot( svf.t[[i]], m.best[[i]][[1]], 
              xlim =  xlimz ,
        main = paste( ids[i], an ) )#best model
  #plot two most common models
  a2 <- rownames(summary( m.best[[i]][2]))
  #plot best model
  ctmm::plot( svf.t[[i]], m.best[[i]][[2]], 
              xlim =  xlimz ,
              main = paste( ids[i], a2 ) )#best model
  # #plot against traditional KDE
  # ctmm::plot( svf.t[[i]], m.best[[i]]$"IID isotropic", 
  #             xlim = xlimz,
  #       main = paste( ids[i], "IID isotropic" ) ) 
}  

# Comment on the differences in the variance model assumptions
# between akde and traditional kde ###
# Answer:
#

# How consistent was the top model chosen among individuals?
# For which individuals did it vary most? How?
# Answer:
# 

# Now that we have estimated top movement models for each #
# individual we are ready to apply those models to our estimates #
# of ranges. 

# We also have an extra option to choose from #
# we can weight points based on high utilisation to correct the range #
# estimate.  Weighing may be helpful if there are large data gaps or multiple #
# sampling rates #

# Here we compare ranges from 3 options: (1) top movement model 
# (2) 2nd top movement model  (3) traditional kde no weighing

# we create objects to store output from our 3 options:
akde.ou <- list()
akde.ouf <- list()
kde.iid <- list()
# We loop through each individual to estimate ranges for each option:
for( i in 1:length(ids) ){
  print(i)
  #extract most common movement model
  akde.ou[[i]] <- ctmm::akde( ctmm.t[[i]], m.best[[i]]$"OU anisotropic",
                              weights = FALSE  )
  # extract second most comment movement model
  akde.ouf[[i]] <- ctmm::akde( ctmm.t[[i]], m.best[[i]]$"OUF anisotropic", 
                        weights = FALSE )
  #extract the IID movement model for comparison
  kde.iid[[i]] <- ctmm::akde( ctmm.t[[i]], m.best[[i]]$"IID isotropic" )
}
#plot estimate ranges comparing output for each option:
par(mfrow = c(3,2))
for( i in 1:length(ids) ){# 2
  print(i)
  plot( ctmm.t[[i]], akde.ouf[[i]] )
  title( paste("OUF model", ids[i]) )
  plot( ctmm.t[[i]], akde.ou[[i]] )
  title("OU model")
  plot( ctmm.t[[i]], kde.iid [[i]])
  title("Traditional KDE" )
}

#extract mean HR estimates for weighted and unweighted approaches 
# as sf polygon and combine 
ouf_list <- list()
ou_list <- list()
i_list <- list()
for( i in 1:length(ids) ){
  #extract home range for each animal and turn into sf object
  sf.ouf <- as.sf( akde.ouf[[i]] )
  sf.ou <- as.sf( akde.ou[[i]] )
  sf.i <- as.sf( kde.iid[[i]] )
  # convert crs to study area (otherwise their crs won't match)
  sf.ouf.t <- st_transform( sf.ouf, crs = get_crs( trks.thin ) ) 
  sf.ou.t <- st_transform( sf.ou, crs = get_crs( trks.thin ) ) 
  sf.i.t <- st_transform( sf.i, crs = get_crs( trks.thin ) ) 
  #extract only the point estimate (mean range) and add to list
  ouf_list[[i]] <- sf.ouf.t[2,]
  ou_list[[i]] <- sf.ou.t[2,]
  i_list[[i]] <- sf.i.t[2,]
}

ouf_akdes <-  ouf_list %>%  dplyr::bind_rows()
ou_akdes <-  ou_list %>%  dplyr::bind_rows()
kde_akdes <-  i_list %>%  dplyr::bind_rows()
#re-add attributes for each individual
iddf <- trks.thin %>% 
  group_by( id ) %>% 
  select( id, territory, sex ) %>% 
  slice(1 )
#view
iddf

class( ouf_akdes)
head(ouf_akdes)
ouf_akdes$name <- iddf$territory
ouf_akdes$id <- iddf$id
ouf_akdes$sex <- iddf$sex
ou_akdes$name <- iddf$territory
ou_akdes$id <- iddf$id
ou_akdes$sex <- iddf$sex
kde_akdes$name <- iddf$territory
kde_akdes$id <- iddf$id
kde_akdes$sex <- iddf$sex

head( ouf_akdes)

#Plot comparisons from the different model choices
ggplot() +
  theme_bw( base_size = 15 ) +
#compare against  ouf model using all data from ctmm
geom_sf( data = ouf_akdes,
         fill = NA, col = "purple", linewidth = 2 ) +
#compare against  ou model using all data from ctmm
geom_sf( data = ou_akdes,
         fill = NA, col = "orange", linewidth = 1 ) +
  #compare against kde model using all data from ctmm
  geom_sf( data = kde_akdes,
           fill = NA, col = "black", linewidth = 1 ) +
  facet_wrap( ~id )

## How do the options compare?
#Answer;
#

#####

##############################################################
# Estimating AKDE using atm package                          #
##############################################################
### DO THIS SECTION FOR HOMEWORK   #

# once we have decided which movement model to use we can use 
# amt directly to plot all home ranges at once: 

#make sure individuals are in order so that they can be compared to ctmm results
nested.thin <- trks.thin %>% 
  arrange( id ) %>% 
  nest( data = -"id" ) #nest tibbles

#calculate home range using your chosen movement model:
akde_all <- nested.thin %>% 
  mutate( hr_akde_all = map( data, ~hr_akde( ., 
      model = fit_ctmm(., model = "ou", 
#          uere = uere, ctmm( isotropic = FALSE) 
                ),
                levels = c(0.95) ) ) )
# 
akde_all

for( i in 1:length(ids) ){
#Plot for all individuals against equivalent from ctmm
a <- ggplot() +
  theme_bw( base_size = 8 ) +
  #extract isopleths for ouf model using thinned data from amt
    geom_sf( data = hr_isopleths( akde_all$hr_akde_all[[i]] ),
              col = "black", linewidth = 1, fill=NA ) +
  ### you won't have this object because it was created 
  # on the ctmm model selection section
    geom_sf( data = ou_akdes %>% filter(id == i),
           fill = NA, 
           col = "purple", linewidth = 1 )
  labs( title = ids[i] ) 
print(a)
}

# What is the discrepancy between the two outlines?
# Answer:
# 

# now we plot them with the study area underneath
akde_all %>%
  #choose one home range method at a time
  hr_to_sf( hr_akde_all, id ) %>% 
  #plot with ggplot
  ggplot( . ) +
  theme_bw( base_size = 17 ) + 
  geom_sf( aes( fill = as.factor(id)) ) +  
           #,linewidth = 0.8, alpha = 0.6 ) +
  geom_sf(data = NCA_Shape, inherit.aes = FALSE, fill=NA ) +
  theme( legend.position = "none" ) +
  #plot separate for each individual
  facet_wrap( ~id )

  
###########################################################
### Save desired results #
#save range for your selected individual using your preferred 
# movement model 
#save range for all individuals in atm
write_rds( akde_all, "Data/akde_all" )

# I save range for all individuals estimated with ctmm and 
#converted to sf object. But you don't need to 
# since you did not run that part of the code
# write_rds( ou_akdes, "Data/ou_akdes" )
# write_rds( ouf_akdes, "Data/ouf_akdes" )


#save workspace if in progress
save.image( 'AKDEresults.RData' )
############# end of script  ##################################
