##################################################################
# Script developed by Jen Cruz to estimate ranges using AKDE     # 
# For this script we rely on Fleming et al.(2015) Ecology 96(5):1182-1188#
# We aim to estimate AKDE using high-resolution data and compare #
# with results from data that has been thinned to remove        #
# autocorrelation                                               #
###################################################################

################## prep workspace ###############################

# Clean your workspace to reset your R environment. #
rm( list = ls() )

# load packages relevant to this script:
library( tidyverse ) #easy data manipulation
library( amt )
library( lubridate ) #easy date manipulation

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

###############################################################
##### Estimate ranges using AKDE continuous-time movement model:#
################################################################

# Let's start  by exploring differences in resolution and auto#
# correlation in our two data choices

# We start by plotting points from both datasets for each individual:
trks.breed %>%  #filter( id == 2 ) %>%
  ggplot(., aes( x = x_, y = y_, color = speed ) ) +
  theme_bw( base_size = 15 ) + 
  geom_point() +
  geom_point( data = trks.thin, 
    aes( x = x_, y = y_, color = speed ), 
    size = 3, shape = 8 ) +
  #labs( title = ids[i], fill = "week", x = "lat") +
  facet_wrap( ~territory, scales = "free" )

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
ind <- 6
# subset the two datasets accordingly
idv.breed <- trks.breed %>% filter( id == ind )
idv.thin <- trks.breed %>% filter( id == ind )
#estimate AKDE for thin data
hr_akde <- idv.thin %>% 
        amt::hr_akde( ., model = fit_ctmm(., "auto" ),
                levels = 0.95 )

# did that take too long? If so, then you may want to skip the #
# model selection and choose a realistic model # 
akde_ou <- idv.thin %>% 
  amt::hr_akde(., model = fit_ctmm(., "ou" ),
               levels = 0.95 )

# Did that save computation time? If so then use the same one for #
# your data-rich autocorrelated version: 
akde_ou_all <- idv.breed %>% 
  amt::hr_akde(., model = fit_ctmm(., "ou" ),
               levels = 0.95 )

#Plot comparisons from the different data choices
ggplot() +
  theme_bw( base_size = 15 ) +
  #extract isopleths for autocorrelated AKDE estimates #
  # as we expect those would be the largest
  geom_sf( data = hr_isopleths( akde_ou_all ), 
           fill = NA, col = "black", size = 2 ) +
  #extract isopleths for AKDE estimates using thinned data
  #estimated with the model selection approach:
  geom_sf( data = hr_isopleths( hr_akde ), 
           fill = NA, col = "grey" ) +
  #extract isopleths for AKDE estimates using single movement #
  # behavior method
  geom_sf( data = hr_isopleths( akde_ou ), 
           fill = NA, col = "blue" ) +
  #add points of autocorrelated data
  #note that we turn them into sf points for plotting
  geom_sf( data = as_sf_points( subset( trks.breed, id == 4) ) )+
  #add thinned points
   geom_sf( data = as_sf_points( subset( trks.thin, id == 4)), 
        col = "red" )


#Estimate areas for each method
amt::hr_area( akde_ou ) 
amt::hr_area( akde_ou_all ) 
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
#save.image( 'homerangeresults.RData' )
############# end of script  ##################################
