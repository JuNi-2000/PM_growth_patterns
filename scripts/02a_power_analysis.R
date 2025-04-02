library(edgeR)
library(standR)
library(ssizeRNA)

######################## Load Data ########################
setwd("C:/Users/Julian/OneDrive - Universitaet Bern/PM_Project/DSP_analysis/data")#laptop
spe <- readRDS("spatial_experiment_raw.rds")

# Extract the tumour samples PRE CHEMO only for first DE analysis
spe <- spe[, spe$SegmentLabel == 'tumour']
spe <- spe[, spe$Pre_Post == 'Pre' & !is.na(spe$Pre_Post == 'Pre')]

# Convert spatial experiment to DGE list for limma voom pipeline 
dge <- SE2DGEList(spe)

# Extract raw counts from SpatialExperiment object
counts_matrix <- counts(spe)

# Estimate mean counts in the control group
# define non-solid as control group
control <- spe$solid_like == FALSE
mu <- na.omit(rowMeans(counts_matrix[, control]))
# extract sample size
sum(spe$solid_like == TRUE)
sum(spe$solid_like == FALSE)

# define design
design <- model.matrix(~0 + solid_like, data = colData(spe))
# Estimate dispersion
dge <- DGEList(counts = counts_matrix)
dge <- calcNormFactors(dge)
dge <- estimateDisp(dge, design = design)
disp <- dge$common.dispersion  # Use common dispersion

######################## Power analysis ########################
library(ggplot2)
library(dplyr)

# Define parameters
pi0_values <- c(0.8, 0.9, 0.95)
fc_values <- seq(1.2, 2.0, by = 0.1)

# Initialize an empty data frame to store results
power_results_df <- data.frame()

# Loop over each pi0 value
for (pi0_level in pi0_values) {
  
  bh_ave_vector <- numeric(length(fc_values)) # Store bh_ave for this pi0
  
  # Loop over fold change values
  for (i in seq_along(fc_values)) {
    fc_value <- fc_values[i]
    
    power_results <- check.power(
      nGenes = 10436,   # Total genes
      pi0 = pi0_level,  # Proportion of non-DE genes
      m = 13,           # Sample size per group
      mu = mean(mu),    # Mean counts in control group
      disp = disp,      # Dispersion estimate
      fc = fc_value,    # Fold change for DE genes
      up = 0.5,         # 50% of DE genes are upregulated
      fdr = 0.05,       # Adjusted p-value cutoff
      sims = 10         # Number of simulations
    )
    
    # Store bh_ave value
    bh_ave_vector[i] <- power_results$pow_bh_ave  # Adjust if needed
    
    print(paste0("pi0 = ", pi0_level, ", FC = ", fc_value, " done"))
  }
  
  # Create a temporary data frame for this pi0 level
  temp_df <- data.frame(
    FC = fc_values,
    Power = bh_ave_vector,
    pi0 = factor(pi0_level)  # Convert to factor for ggplot grouping
  )
  
  # Append to main results data frame
  power_results_df <- rbind(power_results_df, temp_df)
}

# Plot results using ggplot2
ggplot(power_results_df, aes(x = FC, y = Power, color = pi0, group = pi0)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "Power Analysis Across Different Fold Changes and pi0 Levels",
    x = "Fold Change (FC)",
    y = "Power (bh_ave)",
    color = "Proportion of Non-DE Genes (pi0)"
  ) +
  geom_hline(yintercept = 0.75, linetype = "dashed")+
  theme_minimal()
