library(tidyverse)
library(BayesPrism)
library(pheatmap)
library(DescTools)


# Bayesprism needs 4 input files to work:
# RAW count matrix of bulkRNA seq data
# RAW single cell count matrix
# Cell type labels containing cell type info for each single cell
# Cell state labels containing fine grained cell type info (subclusters)

# Load dataset and setup workspace
setwd("C:/Users/julian/OneDrive - Universitaet Bern/PM_Project/DSP_analysis/data/PRISM")
sc_count_matrix <- read.delim("./count_matrix_PRISM_raw.tsv")
pseudo_bulk <- read.delim("./pseudo_bulk_PRISM_raw.tsv")
cell_type_info <- read.delim("./cell_type_info.tsv", sep = " ")

# Change cycling mesothelial cells to state
cell_type_info$cell.type[cell_type_info$cell.type == "Cycling mesothelial cells"] <- "Mesothelial cells"

# Transpose the matrices (format needed for Prism)
bk.dat <- t(pseudo_bulk)
sc.dat <- t(sc_count_matrix)
# Extract the cell state and cell type labels
cell.type.labels <- cell_type_info[, 1]
names(cell.type.labels) <- rownames(cell_type_info)
cell.state.labels <- cell_type_info[, 2]
names(cell.state.labels) <- rownames(cell_type_info)


#
plot.cor.phi(input=sc.dat,
              input.labels=cell.state.labels,
              title="cell state correlation",
              #specify pdf.prefix if need to output to pdf
              #pdf.prefix="gbm.cor.cs", 
              cexRow=0.2, cexCol=0.2,
              margins=c(2,2))

# Visualize outliers in single cell RNA data
sc.stat <- plot.scRNA.outlier(
  input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID 
  cell.type.labels=cell.type.labels,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE #return the data used for plotting. 
  #pdf.prefix="gbm.sc.stat" specify pdf.prefix if need to output to pdf
)

# Visualize outliers in bulk
bk.stat <- plot.bulk.outlier(
  bulk.input=bk.dat,#make sure the colnames are gene symbol or ENSMEBL ID 
  sc.input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID 
  cell.type.labels=cell.type.labels,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE
  #pdf.prefix="gbm.bk.stat" specify pdf.prefix if need to output to pdf
)
# Both show highly expressed ribsomal genes which have very low celltype specificity

# Filter out these non-relevant genes
sc.dat.filtered <- cleanup.genes(input=sc.dat,
                                  input.type="count.matrix",
                                  species="hs", 
                                  gene.group=c("Rb","Mrp","other_Rb","chrM","MALAT1") ,
                                  exp.cells=5)



# Check concordance for different types of genes(only makes sense in non-synthetic RNA dataset!)
plot.bulk.vs.sc (sc.input = sc.dat.filtered,
                 bulk.input = bk.dat
                 #pdf.prefix="gbm.bk.vs.sc" specify pdf.prefix if need to output to pdf
)

# Filter protein coding genes only 
sc.dat.filtered.pc <-  select.gene.type (sc.dat.filtered,
                                         gene.type = "protein_coding")

# Select marker genes (optional)
diff.exp.stat <- get.exp.stat(sc.dat=sc.dat[,colSums(sc.dat>0)>3],# filter genes to reduce memory use
                              cell.type.labels=cell.type.labels,
                              cell.state.labels=cell.state.labels,
                              pseudo.count=0.1, #a numeric value used for log2 transformation. =0.1 for 10x data, =10 for smart-seq. Default=0.1.
                              cell.count.cutoff=50, # a numeric value to exclude cell state with number of cells fewer than this value for t test. Default=50.
                              n.cores=1 #number of threads
)
# Apply the marker gene filtering
# Want at least 50 marker genes per cell type! 
# In this case: Had to lower p cutoff value to 0.2 to get >50 for all cells!
# Default is p = 0.01
sc.dat.filtered.pc.sig <- select.marker (sc.dat=sc.dat.filtered.pc,
                                         stat=diff.exp.stat,
                                         pval.max=0.05, 
                                         lfc.min=0.10)

# Construct prism object
myPrism <- new.prism(
  reference=sc.dat.filtered.pc.sig, 
  mixture=bk.dat,
  input.type="count.matrix", 
  cell.type.labels = cell.state.labels, 
  cell.state.labels = cell.state.labels,
  key=NULL,
  outlier.cut=0.01,
  outlier.fraction=0.1,
)


# Run prism
bp.res <- run.prism(prism = myPrism, n.cores=10)
theta <- as.matrix(get.fraction(bp=bp.res,
                       which.theta="final",
                       state.or.type="type"))

bp.res

#saveRDS(bp.res, "results_PRISM.rds")
######################## Check results ########################
setwd("C:/Users/julian/OneDrive - Universitaet Bern/PM_Project/DSP_analysis/data/PRISM")
bp.res <- readRDS("./results_PRISM.rds")
simulated_bulk <- readRDS("./simulate_bulk.rds")

bulk_expression_matrix <- simulated_bulk$simulated_frac
theta <- as.matrix(get.fraction(bp=bp.res,
                                which.theta="first",
                                state.or.type="state"))


colnames(theta) <- gsub(" ", "_", colnames(theta))
# Keep same colnames only for simplified analysis
common_columns <- intersect(colnames(bulk_expression_matrix), colnames(theta))
bulk_expression_matrix <- bulk_expression_matrix[, common_columns]
theta <- theta[, common_columns]


# Plot the correlation per celltype
cor_per_celltype <- sapply(1:ncol(theta), function(j) {
  cor(theta[, j], bulk_expression_matrix[, j], method = "pearson")
})
names(cor_per_celltype) <- colnames(theta)

par(mar=c(11, 4, 4, 4))
barplot(cor_per_celltype,
        las = 2,              # Rotate x-axis labels vertically
        cex.names = 0.8,      # Slightly reduce label font size
        col = "lightgreen",    # Nicer bar color
        ylim = c(0, 1),
        main = "BayesPrism Per Cell Type Correlation",
        ylab = "Pearson Correlation")
legend("topleft", print(paste0("Mean: ", round(median(na.omit(cor_per_celltype)), 2), 
                               "\nMedian: ", round(median(na.omit(cor_per_celltype)), 2), 
                               "\nSD: ", round(sd(na.omit(cor_per_celltype)), 2))), 
       bty = "n")

# Get the concordance correlation coefficient
ccc <- sapply(1:ncol(bulk_expression_matrix), function(i) {
  CCC(bulk_expression_matrix[, i], theta[, i])$rho.c
})
ccc <- unlist(ccc[1, ])

# Plot the Concordance Correlation coefficient
names(ccc) <- colnames(theta)
barplot(ccc, las = 2, col = "skyblue", 
        main = "BayesPrism Concordance per Cell Type", 
        ylab = "CCC (Concordance Correlation Coefficient)", 
        ylim = c(0, 1), 
        cex.names = 0.8)
legend("topleft", print(paste0("Mean: ", round(mean(na.omit(ccc)), 2), 
                               "\nMedian: ", round(median(na.omit(ccc)), 2), 
                               "\nSD: ", round(sd(na.omit(ccc)), 2))), 
       bty = "n")


S# Absolute difference matrix
diff_matrix <- abs(theta - bulk_expression_matrix) / bulk_expression_matrix
# Cap max value at 5
diff_matrix_capped <- pmin(diff_matrix, 1)

# Heatmap of differences
pheatmap(diff_matrix_capped, main = "Relative Difference per Cell (absolute difference /\n true cell fractions)")


