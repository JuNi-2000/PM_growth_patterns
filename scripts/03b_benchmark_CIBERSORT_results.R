######################## Setting up workspace ########################
library(tidyverse)
library(pheatmap)
library(DescTools)


setwd("C:/Users/Julian/OneDrive - Universitaet Bern/PM_Project/DSP_analysis/data")
# load data
simulated_bulk <- readRDS("./CIBERSORT/simulate_bulk.rds")
results_CIBERSORT <- read.csv("./CIBERSORT/CIBERSORTx_Job2_Results.csv", row.names = 1)

# Extract the simulated "true" fractions
simulated_fractions <- as.data.frame(simulated_bulk["simulated_frac"])
# Extract the cell composition only
CIBERSORT_matrix <- results_CIBERSORT[, 1:(ncol(results_CIBERSORT)-3)]

# Match the column names of the two data frames
colnames(CIBERSORT_matrix)
colnames(simulated_fractions)
clean_names <- as.vector(colnames(simulated_fractions) |>
  str_remove("^simulated_frac\\."))
clean_names <- gsub("\\.", "_", clean_names)
# Manually curate last ones
clean_names[clean_names == "Exhausted_T_cells_CD8_"] <- "Exhausted_T_cells_CD8"
clean_names[clean_names == "Cytotoxic_T_cells_CD8_"] <- "Cytotoxic_T_cells_CD8"
# Check
all(clean_names %in% colnames(CIBERSORT_matrix))
# Apply clean names
colnames(simulated_fractions) <- clean_names
# Add 0 entries for mast cells in the simulation (there were none in the dataset)
simulated_fractions <- cbind(Mast_cells = as.numeric(10e-8), simulated_fractions)

#reorder the dataframes so they match
CIBERSORT_matrix <- as.matrix(CIBERSORT_matrix[, colnames(simulated_fractions)])
simulated_fractions <- as.matrix(simulated_fractions)


# Absolute difference matrix
diff_matrix <- abs(CIBERSORT_matrix - simulated_fractions) / simulated_fractions
# Cap max value at 5
diff_matrix_capped <- pmin(diff_matrix, 1)

# Heatmap of differences
pheatmap(diff_matrix_capped, main = "Relative Difference per Cell (absolute difference /\n true cell fractions)")

ccc <- sapply(1:ncol(simulated_fractions), function(i) {
  CCC(simulated_fractions[, i], CIBERSORT_matrix[, i])$rho.c
})
ccc <- unlist(ccc[1, ])

# Plot the Concordance Correlation coefficient
par(mar=c(12,4, 4, 4))
names(ccc) <- colnames(CIBERSORT_matrix)
barplot(ccc, las = 2, col = "skyblue", 
        main = "CIBERSORT Concordance per Cell Type", 
        ylab = "CCC (Concordance Correlation Coefficient)",
        cex.names = 0.8,
        ylim = c(0, 1))
legend("topleft", print(paste0("Mean: ", round(mean(na.omit(ccc)), 2), 
                               "\nMedian: ", round(median(na.omit(ccc)), 2), 
                               "\nSD: ", round(sd(na.omit(ccc)), 2))), 
       bty = "n")


# Plot the correlation per celltype
cor_per_celltype <- sapply(1:ncol(CIBERSORT_matrix), function(j) {
  cor(CIBERSORT_matrix[, j], simulated_fractions[, j], method = "pearson")
})
names(cor_per_celltype) <- colnames(CIBERSORT_matrix)
barplot(cor_per_celltype, las = 2, main = "CIBERSORT Per Cell Type Correlation", 
        ylim = c(0, 1), 
        col = "lightgreen")
legend("topleft", print(paste0("Mean: ", round(mean(na.omit(cor_per_celltype)), 2),
                               "\nMedian: ", round(median(na.omit(cor_per_celltype)), 2), 
                               "\nSD: ", round(sd(na.omit(cor_per_celltype)), 2))), 
                        bty = "n")




# Plot the correlation per sample
cor_per_sample <- sapply(1:nrow(CIBERSORT_matrix), function(i) {
  cor(CIBERSORT_matrix[i,], simulated_fractions[i,], method = "pearson", 
      use = "complete.obs")
})
names(cor_per_sample) <- rownames(CIBERSORT_matrix)
barplot(cor_per_sample, las = 2, main = "CIBERSORT Per Cell Type Correlation", 
        ylim = c(0, 1))

par(mfrow = c(3, 6), mar = c(4, 4, 2, 1))  # mar tweaks margin size
for (j in 1:ncol(CIBERSORT_matrix)) {
  plot(simulated_fractions[, j], CIBERSORT_matrix[, j],
       xlab = "True", ylab = "Estimated",
       main = colnames(CIBERSORT_matrix)[j],
       pch = 16, col = "steelblue")
  abline(0, 1, col = "red", lty = 2)
}










