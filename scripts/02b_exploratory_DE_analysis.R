library(ggplot2)
library(ggalluvial)
library(edgeR)
library(limma)
library(tidyverse)


######################## SETTING UP WORKSPACE ########################
# Setting working directory
setwd("C:/Users/julia/OneDrive - Universitaet Bern/PM_Project/DSP_analysis/data")

# Load Spatial experiment file
spe <- readRDS("spatial_experiment_raw.rds")


######################## Differential gene expression ########################

# Extract the tumour samples PRE CHEMO only for first DE analysis
spe <- spe[, spe$SegmentLabel == 'tumour']
#spe <- spe[, spe$Pre_Post == 'Pre' & !is.na(spe$Pre_Post == 'Pre')]

# Convert spatial experiment to DGE list for limma voom pipeline 
dge <- SE2DGEList(spe)
dim(dge)

dge <- dge[, !is.na(dge$samples$Pre_Post)]


# Construct the design matrix for DE analysis
design <- model.matrix(~0 + growthPattern + Pre_Post, data = colData(spe))


# Filter low expression genes
keep <- filterByExpr(dge, design)
table(keep) # None of the genes are highlighted => no filtering necessary


# Calculate normalisation factors for TMM 
dge <- calcNormFactors(dge)

# Check ratio of library sizes
max(spe$lib_size) / min (spe$lib_size) # Ratio is big!
# Non-consistant library size => use voom approach!!

# Apply voom transformation to the normalized DGElist object
v <- voom(dge, design, plot = TRUE, normalize.method = "quantile")

# Estimate the correlation between measurements from same patients
corfit <- duplicateCorrelation(v, design, block=interaction(dge$samples$PID_anon, dge$samples$SlideName))

# Fit the linear model
fit <- lmFit(v, design, block=interaction(dge$samples$PID_anon, dge$samples$SlideName), correlation = corfit$consensus)

# Create contrast matrix for comparisons between growth patterns
contrast_matrix <- makeContrasts(
  solid_vs_papillary = growthPattern1 - growthPattern2,
  solid_vs_glandular = growthPattern1 - growthPattern3,
  papillary_vs_glandular = growthPattern2 - growthPattern3, 
  levels = design)


fit_contrast <- contrasts.fit(fit, contrasts = contrast_matrix)
efit <- eBayes(fit_contrast, robust = TRUE)

results_efit <- decideTests(efit, adjust.method = "fdr", p=0.05) 
summary(results_efit)

##### CHANGE FILTERING THRESHOLD??


# BGV check
dge <- estimateDisp(dge, design = design, robust = TRUE) # error: no samples belong to category 4. Remove cat 4 in this case
design <- design[, c(1:3)]
# Try again
dge <- estimateDisp(dge, design = design, robust = TRUE) 

# Plot BGV 
plotBCV(dge, legend.position = "topleft", ylim = c(0, 1.3))
bcv_df <- data.frame(
  'BCV' = sqrt(dge_all$tagwise.dispersion),
  'AveLogCPM' = dge_all$AveLogCPM,
  'gene_id' = rownames(dge_all)
)

highbcv <- bcv_df$BCV > 0.8
highbcv_df <- bcv_df[highbcv, ]
points(highbcv_df$AveLogCPM, highbcv_df$BCV, col = "red")
text(highbcv_df$AveLogCPM, highbcv_df$BCV, labels = highbcv_df$gene_id, pos = 4)



# Get results 
library(ggrepel)
library(tidyverse)

de_results_BvT <- topTable(efit, coef = 2, sort.by = "P", n = Inf)

de_genes_toptable_BvT <- topTable(efit, coef = 1, sort.by = "P", n = Inf, p.value = 0.05)


de_results_BvT <- mutate(de_results_BvT, DE = ifelse(logFC > 0 & adj.P.Val <0.05, "UP", ifelse(logFC <0 & adj.P.Val<0.05, "DOWN", "NOT DE")))

ggplot(de_results_BvT, aes(AveExpr, logFC, col = DE)) + 
  geom_point(shape = 2, size = 1) + 
  geom_text_repel(data = de_genes_toptable_BvT %>% 
                    mutate(DE = ifelse(logFC > 0 & adj.P.Val <0.05, "UP", 
                                       ifelse(logFC <0 & adj.P.Val<0.05, "DOWN", "NOT DE"))) %>%
                    rownames_to_column(), aes(label = rowname)) +
  theme_bw() +
  xlab("Average log-expression") +
  ylab("Log-fold-change") +
  ggtitle("B cell zone vs. T cell zone in Lymph node (limma-voom)") +
  scale_color_manual(values = c("blue","gray","red")) +
  theme(text = element_text(size=15))



######################## GO enrichment analysis ########################








