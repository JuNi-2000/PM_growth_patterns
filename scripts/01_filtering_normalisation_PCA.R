library(tidyverse)
library(readxl)
library(standR)
library(SpatialExperiment)

######################## SETTING UP WORKSPACE ########################
# Setting working directory
setwd("C:/Users/julia/OneDrive - Universitaet Bern/PM_Project/DSP_analysis/data")
setwd("C:/Users/Julian/OneDrive - Universitaet Bern/PM_Project/DSP_analysis/data") #laptop

# Load data
count_matrix <- as.data.frame(read_excel("./DSP_22_09_meso_RawData.xlsx", sheet=3))
sample_anno <- as.data.frame(read_excel("./DSP_22_09_meso_RawData.xlsx", sheet=1))
feature_anno <- as.data.frame(read_excel("./DSP_22_09_meso_RawData.xlsx", sheet=2))
morphologies <- read_delim("./Meso_morphology_TMA_22_09.csv", delim = ";")
metadata <- read_delim("./metadata_merged_anon.tsv", delim = "\t")


############### DATA CLEANING ####################
# Extract the legend for easy reference
morphology_legend <- names(morphologies)

# Creating clean data structure for morphology file
names(morphologies) <- c("SlideName", "core.x", "core.y", "growthPattern", "secondaryGrowthPattern")

# Secondary growthpattern is NA so remove it
morphologies <- morphologies[, -ncol(morphologies)]

# Make sure the growth patterns are numeric entries only
morphologies$growthPattern <- as.numeric(morphologies$growthPattern)

# Remove NA entries for primary growth patterns if present
morphologies <- drop_na(morphologies, c("growthPattern"))

# Fill NAs in secondary growth patterns with 0 to prevent working with NA entries
#morphologies$secondaryGrowthPattern[is.na(morphologies$secondaryGrowthPattern)] <- 0

# Merge morphologies based on prognosis: Solid, sarcomatoid and micropapillary have bad prognosis
# tubulopapillary and trabecular have more favourable prognosis
morphologies <- morphologies %>%
  mutate(solid_like = growthPattern %in% c(1, 4, 5))
           
# Merging morphologies with sample_anno file
sample_anno_merged <- merge(sample_anno, morphologies, by = c("core.x", "core.y", "SlideName"), all = TRUE)

# Merging pre / post chemo and donor block ID information with sample anno file
relevant_metadata <- metadata[, c("SegmentDisplayName", "Pre_Post", "PID_anon", "Tumour_region", "Donor.Block.ID_anon")]
sample_anno_merged <- merge(sample_anno_merged, relevant_metadata, by = "SegmentDisplayName", all= TRUE)

######################## LOQ based Pre-filtering  ########################

# Check if the names of the anno file match the count matrix order (excl. the first col with the gene names)
all(names(count_matrix)[-1] == sample_anno$SegmentDisplayName)

# Extract the LOQ values from the annotation file
LOQ <- as.data.frame(t(sample_anno$`ExpressionFilteringThreshold (Human NGS Whole Transcriptome Atlas RNA_1.0)`))
names(LOQ) <- sample_anno$SegmentDisplayName

# Merge the LOQ values into the count matrix as a new row
LOQ_count_matrix <- rbind(count_matrix[-1], 
                          sample_anno$`ExpressionFilteringThreshold (Human NGS Whole Transcriptome Atlas RNA_1.0)`)

# Create empty matrix to store a binary mask
Mask <- matrix(0, nrow=nrow(LOQ_count_matrix)-1, ncol = ncol(LOQ_count_matrix))

# Create binary mask indicating if LOQ is reached for each gene
for(i in seq(ncol(LOQ_count_matrix))){
  col <- LOQ_count_matrix[, i]
  LOQ <- col[length(col)]
  above_LOQ <- sapply(col, function(x) as.numeric(x > LOQ))
  Mask[, i] <- above_LOQ[1:length(above_LOQ)-1]
}

# Sample level QC: Decide for each sample how many genes are above LOQ
threshold_ROI <- 1.5 # if less than 1.5% of genes in a sample are above LOQ, its likely bad quality
# Create vector with all the elements to keep
keep_ROI <- vector(length = ncol(Mask))

# Decide for each sample if enough genes are above LOQ
for(i in seq(ncol(Mask))){
  col <- Mask[, i]
  nr_genes <- length(col)
  colsum <- sum(col)
  above_LOQ <- ((colsum / nr_genes)*100) >= threshold_ROI
  keep_ROI[i] <- above_LOQ
}
# Inspect results
table(keep_ROI)
# Apply ROI filtering to Mask 
Mask <- Mask[, keep_ROI]

# Gene wise filtering: Genes above LOQ for at least 6% of samples are kept
threshold_genes <- 6 # 6 percent
keep_genes <- vector(length = nrow(Mask))
for(i in seq(nrow(Mask))){
  gene_row <- Mask[i,]
  nr_samples <- length(gene_row)
  rowsum <- sum(gene_row)
  above_LOQ <- ((rowsum / nr_samples)*100) >= threshold_genes
  keep_genes[i] <- above_LOQ
}
table(keep_genes)


######################## Starting standR workflow ########################
# Create the spatial experiment object
spatial_experiment <- readGeoMx(count_matrix, 
                                sample_anno_merged,
                                feature_anno, 
                                rmNegProbe = FALSE)

# Apply the LOQ based pre filtering to spatial experiment object
spatial_experiment <- spatial_experiment[, keep_ROI]
spatial_experiment <- spatial_experiment[keep_genes, ]

# Remove NA entries in the growth patterns
spatial_experiment <- spatial_experiment[, !is.na(colData(spatial_experiment)$growthPattern)]


######################## standR QUALITY CONTROL ########################
# Checking QC flags
colData(spatial_experiment)$QCFlags
# Most samples show low NTC or low neg. probe count

# Sample level QC
plotSampleInfo(spatial_experiment, column2plot = c("SlideName", "growthPattern"))

# Add statistics to spatial experiment
# Removes genes with low count in more than 90% of the samples
# Removes genes with less than a minimum count of 5
spatial_experiment <- addPerROIQC(spatial_experiment, rm_genes = TRUE, sample_fraction = 0.8)
dim(spatial_experiment) # Check if any ROIs were removed

# Investigate percentage of lowly expressed genes in each sample
plotGeneQC(spatial_experiment)

# ROI-level QC: Plot library size and cell numbers. Cutoff threshold usually recommended between
# 100 and 200, depending on the analysis method
plotROIQC(spatial_experiment, x_threshold = 150, color= SlideName)

qc <- colData(spatial_experiment)$AOINucleiCount > 150
table(qc) #20 ROIs get filtered out with a cell number threshold of 150

# Apply the threshold to the data
spatial_experiment <- spatial_experiment[, qc]

# QC for the AreaSize
plotROIQC(spatial_experiment,  x_axis = "AOISurfaceArea", x_lab = "AreaSize", 
          y_axis = "lib_size", y_lab = "Library size", col = SlideName, 
          x_threshold = 1.6e4) # No area filtering necessary for recommended threshold of 16000

# Visualize the relative log expression (identify technical variation in the data)
plotRLExpr(spatial_experiment)
plotRLExpr(spatial_experiment, ordannots = "growthPattern", assay = 2, color = growthPattern)
plotRLExpr(spatial_experiment, ordannots = "SlideName", assay = 2, color = SlideName)
plotRLExpr(spatial_experiment, ordannots = "SegmentLabel", assay = 2, color = SegmentLabel)
# Variation is more or less uniform / random across these conditions => no evidence of systematic bias

# Plot PCA to further investigate systematic differences
spatial_experiment$growthPattern <- as.factor(spatial_experiment$growthPattern)
spatial_experiment$PID_anon <- as.factor(spatial_experiment$PID_anon)
drawPCA(spatial_experiment, assay = 2, color = solid) #Clear separation between tumour and TME


# Split tumour and TME to see if patterns show clear signature within these two groups
tme <- colData(spatial_experiment)$SegmentLabel == 'tme'
spe_TME_subset <- spatial_experiment[, tme]
drawPCA(spe_TME_subset, assay = 2, color = solid) # No clear separation visible between groups

tumour <- colData(spatial_experiment)$SegmentLabel == 'tumour'
spe_tumour_subset <- spatial_experiment[, tumour]
spe_tumour_subset <- spe_tumour_subset[, !is.na(spe_tumour_subset$Pre_Post)]
spe_tumour_subset <- spe_tumour_subset[, spe_tumour_subset$Pre_Post == "Pre"]
drawPCA(spe_tumour_subset, assay = 2, color = PID_anon) 

## Investigate other principal components
#spe_tumour_subset <- scater::runPCA(spe_tumour_subset)
#pca_results <- reducedDim(spe_tumour_subset, "PCA")
#plotPairPCA(spe_tumour_subset, col = Pre_Post, precomputed = pca_results, 
#            n_dimension = 4, title = "PCA plots for tumour compartment only")
#

spatial_experiment <- scater::runUMAP(spatial_experiment, dimred = "UMAP")
plotDR(spatial_experiment, dimred = "UMAP", col = growthPattern)

# Save the file for subsequent analysis (DE analysis requires non normalized data)
saveRDS(spatial_experiment, file = "spatial_experiment_raw.rds")

######################## Aggregating count matrix ########################
col_data <- colData(spatial_experiment)

# Convert colData to a tibble for easier manipulation
col_data_tbl <- as_tibble(col_data)

# Perform aggregation
aggregated_data <- col_data_tbl %>%
  group_by(Donor.Block.ID_anon, Pre_Post, Tumour_region) %>%
  summarise(across(everything(), mean, na.rm = TRUE), .groups = 'drop')


# Example: Creating a new SpatialExperiment object with aggregated data
new_spe <- spatial_experiment[, match(aggregated_data$Donor.Block.ID_anon, col_data$Donor.Block.ID_anon)]
colData(new_spe) <- DataFrame(aggregated_data)



###################### Normalisation ##############################

# Normalize using TMM
spa_experiment_TMM <- geomxNorm(spatial_experiment, method = "TMM")

# Investigate the result of the normalisation 
plotRLExpr(spa_experiment_TMM, assay = 2, color = SlideName) + ggtitle("TMM normalized")


###################### Removing batch effect #########################
# Since there is still a small batch effect visible in the normalized expression data we want to remove it
spa_experiment_TMM<- findNCGs(spa_experiment_TMM, batch_name = "SlideName", top_n = 300)
metadata(spa_experiment_TMM) |> names()

# Run the geomxBatchCorrection on the normalized data
# ADD PRE AND POST THERAPY AS WELL TO BIOLOGICAL FACTORS??
for(i in seq(5)){
  spe_ruv <- geomxBatchCorrection(spa_experiment_TMM, factors = c("SegmentLabel","growthPattern"), 
                                  NCGs = metadata(spa_experiment_TMM)$NCGs, k = i)
  
  print(plotPairPCA(spe_ruv, assay = 2, n_dimension = 4, color = SegmentLabel, title = paste0("k = ", i)))
  
}

# Using k=1 for this analysis
spe_ruv <- geomxBatchCorrection(spa_experiment_TMM, factors = c("SegmentLabel","growthPattern"),  
                                NCGs = metadata(spa_experiment_TMM)$NCGs, k = 1)
# Assess the expression after batch correction
plotRLExpr(spe_ruv, assay = 2, color = SlideName) + ggtitle("RUV4 Removed batch effect and TMM normalized")


# Compare batch correction to the limma removeBatchEffect
spe_lrb <- geomxBatchCorrection(spa_experiment_TMM,
                                batch = colData(spa_experiment_TMM)$SlideName, method = "Limma",
                                design = model.matrix(~growthPattern + SegmentLabel, data = colData(spa_experiment_TMM)))

plotRLExpr(spe_lrb, assay = 2, color = SlideName) + ggtitle("LIMMA Removed batch effect and TMM normalized")
# RUV4 looks better in this case => less deviation of the sample means from 0

# There are some statistics to formally validate which method perfoms better on this data
spe_list <- list(spa_experiment_TMM, spe_ruv, spe_lrb)
plotClusterEvalStats(spe_list = spe_list,
                     bio_feature_name = c("SegmentLabel"),
                     batch_feature_name = "SlideName",
                     data_names = c("Raw","RUV4","Limma"))


# Save the normalized batch corrected data
saveRDS(spe_ruv, file = "spatial_experiment_normalized_batch_corrected.rds")


######################## Aggregating technical replicates ########################
# Load the spatial experiment file
spe_agg <- readRDS("./spatial_experiment_normalized_batch_corrected.rds")
# Extract the counts and features data
counts_t <- as.data.frame(t(spe_agg@assays@data@listData$counts))
features <- as.data.frame(colData(spe_agg))
features <- features[match(rownames(counts_t), rownames(features)), ]

# Extract relevant entries for the aggregation of technical replicates
rel_features <- features[, c("Pre_Post", "Tumour_region", "Donor.Block.ID_anon", 
                             "SegmentLabel", "solid_like", "growthPattern", 
                             "PID_anon")]

# Bind counts to features for easy aggregation
rel_features <- cbind(rel_features, counts_t)

# Keep treatment naive samples only
rel_features <- drop_na(rel_features[rel_features$Pre_Post == "Pre", ])

# Check if growth patterns are the same with, in triplets (technical replicates)
# If not the same, aggregation might not make sense
cols_to_check <- c("Donor.Block.ID_anon", "Tumour_region", "Pre_Post", "SegmentLabel")
df_with_multiple_growth <- rel_features %>%
  group_by(across(all_of(cols_to_check))) %>%
  summarise(n_growth_patterns = n_distinct(growthPattern), .groups = 'keep') %>%
  filter(n_growth_patterns > 1)
print(paste0(nrow(df_with_multiple_growth), " out of ", 
             nrow(rel_features), " samples have ambiguous growth patterns"))

# Remove these ambiguous samples


# Collapse technical replicates by taking mean over the gene counts
aggregated_df <- rel_features %>%
  group_by(Donor.Block.ID_anon, Tumour_region, Pre_Post, SegmentLabel) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  ungroup()

#


# Extract the aggregated counts matrix
counts_agg <- as.matrix(t(aggregated_df[, 5:nrow(aggregated_df)]))

######################## GSVA ########################
library(GSVA)
library(GSEABase)
library(limma)
library(pheatmap)

# Start with hallmark gene set
hallmark_gene_sets <- getGmt("https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5.1/h.all.v7.5.1.symbols.gmt")
# Want to look at pre treatment only
counts_pre_treatment <- aggregated_df[aggregated_df$Pre_Post == "Pre", ]
# Extract tumour and tme compartment
counts_pre_treatment_tumor <- na.omit(counts_pre_treatment[counts_pre_treatment$SegmentLabel == "tumour", ])
counts_pre_treatment_tme <- na.omit(counts_pre_treatment[counts_pre_treatment$SegmentLabel == "tme", ])

# Extract the aggregated counts matrix
counts_pre_treatment_tumor <- as.matrix(t(counts_pre_treatment_tumor[, 5:nrow(counts_pre_treatment_tumor)]))
counts_pre_treatment_tme <- as.matrix(t(counts_pre_treatment_tme[, 5:nrow(counts_pre_treatment_tme)]))

# Build the input objects for GSVA
gsva_input_pre_tumour <- gsvaParam(counts_pre_treatment_tumor, hallmark_gene_sets)
gsva_input_pre_tme <- gsvaParam(counts_pre_treatment_tme, hallmark_gene_sets)

# Run GSVA
gsva_pre_tumour <- gsva(gsva_input_pre_tumour, verbose = TRUE)
gsva_pre_tme <- gsva(gsva_input_pre_tme, verbose = TRUE)


# Visualize results: Tumour compartment
pheatmap(gsva_pre_tumour, scale = "row", cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = TRUE, show_colnames = TRUE, 
         main = "Hallmark GSVA pre-chemotherapy on \n aggregated tumour compartment")

# Visualize results: TME compartment
pheatmap(gsva_pre_tme, scale = "row", cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = TRUE, show_colnames = TRUE, 
         main = "Hallmark GSVA pre-chemotherapy on \n aggregated TME compartment")



### Oncogenic signatures ###
oncogenic_signatures_genes <- getGmt("https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5.1/c6.all.v7.5.1.symbols.gmt")
gsva_input_pre_tumour <- gsvaParam(counts_pre_treatment_tumor, oncogenic_signatures_genes)
gsva_pre_tumour <- gsva(gsva_input_pre_tumour, verbose = TRUE)
pheatmap(gsva_pre_tumour, scale = "row", cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = TRUE, show_colnames = TRUE, 
         main = "C6 Oncogenic signatures pre-chemotherapy on \n aggregated tumour compartment")





