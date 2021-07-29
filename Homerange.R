#################################################################
# Script developed by Jen Cruz to calculate home ranges        # 
# We rely heavily on amt vignette here:       #
# https://cran.r-project.org/web/packages/amt/vignettes/p2_hr.html #

###################################################################

################## prep workspace ###############################


# load packages relevant to this script:
library( tidyverse ) #easy data manipulation
library( amt )
library( sp )

#####################################################################
## end of package load ###############

###################################################################
#### Load or create data -----------------------------------------
# Clean your workspace to reset your R environment. #
rm( list = ls() )

# # Set working directory. This is the path to your Rstudio folder for this 
# # project. If you are in your correct Rstudio project then it should be:
# getwd()
# # if so then:
# workdir <- getwd()
# 
# # set path to where you can access your data #
# # Note that the path will be different in yours than mine.#
# datapath <- "Z:/Common/PrairieFalcons/"
# 
# #import GPS data
# dataraw <- read.csv( file = paste0( datapath,"ctt_data_export_20210721152221.csv" ),
#                      header = TRUE,
#                      colClasses = c("serial" = "character"))
# #view
# head( dataraw ); dim( dataraw )
# load workspace 
load( "TracksWorkspace.RData" )
##############   Home range          #####################

# estimate home range using a kernel estimator
idv_kde <- amt::hr_kde( tr.idv, levels = c(0.5, 0.9) )
#estimate home range as minimum convex polygon
idv_mcp <- amt::hr_mcp( tr.idv, levels = c(0.5, 0.9) )

#view
par( mfrow = c(1,1))
plot( idv_kde )
plot( idv_akde, add.relocations = FALSE, add = TRUE, lty = 3  )
plot( idv_mcp, add.relocations = FALSE, add = TRUE, lty = 2 )
# Comment on the difference between the two methods
# Answer:
#
#We can also estimate corresponding areas
amt::hr_area( idv_kde ); amt::hr_area( idv_mcp )

#compare methods against reduced data for the same individual:
red_kde <- amt::hr_kde( tr.red, levels = c(0.5, 0.9) )
red_mcp <- amt::hr_mcp( tr.red, levels = c(0.5, 0.9) )
#view
plot( red_kde, add.relocations = FALSE, add = TRUE, 
      lwd = 2, col = 'red' )
plot( red_mcp, add.relocations = TRUE, add = TRUE, 
      col = 'blue', lwd = 2 )
#Estimate areas
amt::hr_area( red_kde ); amt::hr_area( red_mcp )

# What was the influence of removing temporal autocorrelation on 
# results from both methods?
# Answer:
#

###############################################################
##### compare methods for all individuals              ########
##############################################################

# for all data
hrcomp <- trks.comp %>% amt::nest( data = -"id" ) %>% 
  mutate(
    hr_mcp = map(data, ~ hr_mcp(., levels = c(0.5, 0.9)) ),
    hr_kde = map(data, ~ hr_kde(., levels = c(0.5, 0.9)) ),
    n = map_int( data, nrow )
  ) 

#plot homeranges
hrcomp %>%
  #choose one method at a time
  hr_to_sf(  hr_kde, id, n ) %>% 
ggplot( . ) +
  theme_bw( base_size = 17 ) + 
  geom_sf() +
  facet_wrap( ~id )


# for thinned locations:
hrred <- trks.red %>%  amt::nest( data = -"id" ) %>% 
  mutate(
    hr_mcp = map(data, ~ hr_mcp(., levels = c(0.5, 0.9)) ),
    hr_kde = map(data, ~ hr_kde(., levels = c(0.5, 0.9)) ),
    #hr_locoh = map(data, ~ hr_locoh(., n = ceiling(sqrt(nrow(.))))),
    n = map_int( data, nrow )
  )  
#view
hrred 

#plot home ranges
#select tibble 
hrred %>%
  #choose one home range method at a time
  hr_to_sf( hr_kde, id, n ) %>% 
  #plot with ggplot
  ggplot( . ) +
  theme_bw( base_size = 17 ) + 
  geom_sf() +
  #keep separate for each indvidual
  facet_wrap( ~id )

# remove tracking data and convert to long dataframe
hrarea <- hrred %>%  select( -data ) %>% 
  pivot_longer( hr_mcp:hr_kde, names_to = "estimator", 
                values_to = "hr" )

# Calculate area for each method and unnest list
hrarea <- hrarea %>%  mutate( hr_area = map( hr, ~hr_area(.)) ) %>% 
  unnest( cols = hr_area )
#######################################################################
############ Home range overlap ###############################
# amt currently implements methods reviewed by Fieberg & Kochany (2005) #
# hr: proportion of home range of instance i that overlaps with the home #
# range of instance j. This measure does not rely on a UD and is #
# directional (i.e., HRi,j≠HRj,i) and bound between 0 (no overlap) #
# and 1 (complete overlap) #
# phr: Is the probability of instance j being located in the home range #
# of instance i. phr is also directional and bounded between 0 (no #
# overlap) and 1 (complete overlap) #
# vi: The volumetric intersection between two UDs.#
# ba: The Bhattacharyya’s affinity between two UDs. #
# udoi: A UD overlap index. # 
# hd: Hellinger’s distance between two UDs. #





############# end of script  ###########################################