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

# Set working directory. This is the path to your Rstudio folder for this 
# project. If you are in your correct Rstudio project then it should be:
getwd()
# if so then:
workdir <- getwd()

# set path to where you can access your data #
# Note that the path will be different in yours than mine.#
datapath <- "Z:/Common/PrairieFalcons/"

#import GPS data
dataraw <- read.csv( file = paste0( datapath,"ctt_data_export_20210721152221.csv" ),
                     header = TRUE,
                     colClasses = c("serial" = "character"))
#view
head( dataraw ); dim( dataraw )
##############   Home range                      #####################

# use tracking data to create a template raster for the KDE #
trast <- amt::make_trast( tr.slow )

# estimate home range
hr_idv <- amt::hr_kde( tr.slow, trast = trast, levels = c(0.5, 0.9) )
plot( hr_idv)

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