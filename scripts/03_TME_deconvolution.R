
######################## SETTING UP WORKSPACE ########################
# Setting working directory
setwd("C:/Users/julia/OneDrive - Universitaet Bern/PM_Project/DSP_analysis/data") #home
setwd("C:/Users/Julian/OneDrive - Universitaet Bern/PM_Project/DSP_analysis/data")#laptop

# Load Spatial experiment file
spe <- readRDS("spatial_experiment_raw.rds")

######################## Deconvolution of TME ########################

# Extract the TME samples PRE CHEMO only
spe <- spe[, spe$SegmentLabel == 'tme']
spe <- spe[, spe$Pre_Post == 'Pre' & !is.na(spe$Pre_Post == 'Pre')]

