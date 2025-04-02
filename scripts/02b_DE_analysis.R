library(ggplot2)
library(ggalluvial)
library(edgeR)
library(limma)
library(tidyverse)


######################## SETTING UP WORKSPACE ########################
# Setting working directory
setwd("C:/Users/julia/OneDrive - Universitaet Bern/PM_Project/DSP_analysis/data") #home
setwd("C:/Users/Julian/OneDrive - Universitaet Bern/PM_Project/DSP_analysis/data")#laptop

# Load Spatial experiment file
spe <- readRDS("spatial_experiment_raw.rds")

######################## Differential gene expression ########################

# Extract the tumour samples PRE CHEMO only for first DE analysis
spe <- spe[, spe$SegmentLabel == 'tumour']
spe <- spe[, spe$Pre_Post == 'Pre' & !is.na(spe$Pre_Post == 'Pre')]

# Convert spatial experiment to DGE list for limma voom pipeline 
dge <- SE2DGEList(spe)
dim(dge)

# Construct the design matrix for DE analysis
design <- model.matrix(~0 + solid_like, data = colData(spe))


# Filter low expression genes
keep <- filterByExpr(dge, design)
table(keep) # None of the genes are highlighted => no filtering necessary
dge <- dge[keep, ]

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
  non_solid_vs_solid = solid_likeTRUE - solid_likeFALSE, 
  levels = design)

fit_contrast <- contrasts.fit(fit, contrasts = contrast_matrix)
efit <- eBayes(fit_contrast)

results_efit <- decideTests(efit, p=0.1)
summary(results_efit)

topTable(efit, sort.by = "logFC")
# How many genes are significant before multiple testing correction
sum(efit$p.value <= 0.05)
######################## Biological variation ########################
# BGV check
dge <- estimateDisp(dge, design = design, robust = TRUE) # error: no samples belong to category 4. Remove cat 4 in this case
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

de_results_BvT <- topTable(efit, sort.by = "P", n = Inf)

#top 50 genes logFold change
sort(abs(de_results_BvT$logFC), decreasing = TRUE)[1:50]


de_genes_toptable_BvT <- topTable(efit, coef = 1, sort.by = "P", n = Inf, p.value = 0.05)


de_results_BvT <- mutate(de_results_BvT, DE = ifelse(logFC > 0 & adj.P.Val <0.05, "UP", ifelse(logFC <0 & adj.P.Val<0.05, "DOWN", "NOT DE")))

ggplot(de_results_BvT, aes(AveExpr, logFC, col = DE)) + 
  geom_point(shape = 1, size = 1) + 
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



######################## GSVA ########################
library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)


# Get the names of the genes
gene_list_logFC <- topTable(efit, sort.by = "logFC", n = Inf)[1]
gene_list_logFC <- rownames(gene_list_logFC)
symbol_to_Entrez <- bitr(gene_list_logFC, fromType = "SYMBOL", toType = "ENTREZID", OrgDb="org.Hs.eg.db") 

# Keep only converted entries
gene_list_logFC <- symbol_to_Entrez[, 2]
gene_list_logFC <- as.vector(gene_list_logFC)
# Extract the logFC values and create the gene list for ReactomePA
logFC_values <- topTable(efit, sort.by = "logFC", n = Inf)[, c("logFC")]
names(logFC_values) <- gene_list_logFC

# Check the gene list
head(logFC_values)
logFC_values <- sort(logFC_values, decreasing = TRUE)

#remove NA entries
logFC_values <- logFC_values[!is.na(names(logFC_values))]

# Check enriched Reactome
reactome <- gsePathway(geneList = na.omit(logFC_values), 
                       organism = 'human',
                       pvalueCutoff = 0.05)

gene_ontology_results <- gseGO(geneList = logFC_values, 
                               OrgDb = org.Hs.eg.db, 
                               )


dotplot(reactome, showCategory = 15) +
  ggtitle("Top Enriched Reactome Pathways")


ridgeplot(reactome, showCategory = 15) +
  ggtitle("Reactome solid vs. non-solid ")

heatplot(reactome, showCategory = 15)
cnetplot(reactome, showCategory = 10, foldChange = logFC_values)

# GO results
ridgeplot(gene_ontology_results, showCategory = 15) +
  ggtitle("Reactome Enrichment Results")


dotplot(gene_ontology_results, showCategory = 15, font.size = 9) +
  ggtitle("GSEA of GO Pathways (Upregulated and Downregulated)")


######################## Hallmarks GSEA ########################
library("msigdbr")
msig_h <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene) %>%
  dplyr::rename(ont = gs_name, gene = entrez_gene)

# Run GSEA on Hallmarks
msig_gsea <- GSEA(logFC_values, TERM2GENE = msig_h, scoreType = "pos") %>%
  as_tibble
msig_gsea

# Select top enriched pathways (adjust this as needed)
top_n_pathways <- 20  # Change to the number of pathways you want to show
msig_gsea <- msig_gsea %>% arrange(p.adjust) %>% head(top_n_pathways)

# Plot
ggplot(msig_gsea, aes(x = NES, y = reorder(ID, NES), size = -log10(p.adjust), color = NES)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "Normalized Enrichment Score (NES)", 
       y = "Pathway",
       size = "-log10(Adjusted P-value)",
       color = "NES",
       title = "GSEA Dot Plot Visualization") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        legend.position = "right")
