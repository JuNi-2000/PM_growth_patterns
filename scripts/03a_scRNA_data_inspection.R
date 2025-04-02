######################## Setting up workspace ########################

setwd("C:/Users/Julian/OneDrive - Universitaet Bern/PM_Project/DSP_analysis/data")
setwd("C:/Users/julia/OneDrive - Universitaet Bern/PM_Project/DSP_analysis/data") #home

# Load the datasets
#sc_mesothelial <- readRDS("./single_cell/mesothelial_cells.rds")
#dummy <- read.table("./CIBERSORT/DUMMY_dataset_CIBERSORT.txt", sep="\t", header= TRUE) # Example dataset from Cibersort to check required formatting
sc_integrated_annotated <- readRDS("./single_cell/sc_integrated_annotated.rds")
sc_tcells_annotated <- readRDS("./single_cell/tcells_annotated.rds")

######################## Prepare construction of signature matrix ########################
# Filter dataset for treatment naive samples only and exclude pleural effusion samples
sc_integrated_annotated <- sc_integrated_annotated[, sc_integrated_annotated$treatment == "naive"]
sc_integrated_annotated <- sc_integrated_annotated[, sc_integrated_annotated$sample_type == "tumor"]
sc_tcells_annotated <- sc_tcells_annotated[, sc_tcells_annotated$treatment == "naive"]
sc_tcells_annotated <- sc_tcells_annotated[, sc_tcells_annotated$sample_type == "tumor"]

# Check results
table(sc_integrated_annotated$treatment)

# Extract the count matrix
# The raw data is saved in "counts", the normalized data is saved as "data" in the Seurat object
# CIBERSORT requires non-log transformed data but Seurat uses log1p transformation
# Normalize raw data using CPM which is accepted by CIBERSORT
count_matrix <- sc_integrated_annotated[["RNA"]]$counts
#Normalize using counts per million
cpm_sig <- apply(count_matrix, 2, function(x) (x / sum(x)) * 1e6)
# Filter out genes with 0 entries across all samples to reduce file size
cpm_sig <- cpm_sig[rowSums(cpm_sig) > 0, ]

# Extract the cell type information and replace T cell classification by detailed one
celltype_fine <- sc_integrated_annotated$celltype_fine
tcell_fine <- sc_tcells_annotated$tcell_annotation
intersection <- intersect(names(celltype_fine), names(tcell_fine))
celltype_fine[intersection] <- tcell_fine[intersection]
table(celltype_fine)

# Remove spaces (CIBERSORT falsely merges cell phenotypes when spaces are present!!!)
celltype_fine <- gsub(" ", "_", celltype_fine)
table(celltype_fine)

# Unname samples to use cell types as new col names (format required for CIBERSORT)
colnames(cpm_sig) <- unname(celltype_fine)

# Convert the sparse matrix to matrix
cpm_sig <- as.matrix(cpm_sig)

# The sample is too big to upload on CIBERSORT
# Idea: Sample randomly from the dataset and create a subset with less samples
set.seed(123)
nr_samples <- ncol(cpm_sig)
random_cell_samples <- sample(1:nr_samples, (nr_samples * 0.5), replace = FALSE) # Keep 50% of the cell samples
cpm_sig_small <- cpm_sig[, random_cell_samples]
table(colnames(cpm_sig_small))

# Manual adjustments: Full dataset only contains 3 mast cells and 8 pericytes
# Want to include them in the signature matrix!
pericytes <- which(colnames(cpm_sig) == "Pericytes")
mast_cells <- which(colnames(cpm_sig) == "Mast_cells")
random_cell_samples <- c(random_cell_samples, pericytes[1:4], mast_cells)

# Kick out duplicates if present after manual curation
random_cell_samples <- unique(random_cell_samples)
cpm_sig_small <- cpm_sig[, random_cell_samples]
table(colnames(cpm_sig_small))

# Convert to data frame
cpm_sig_small <- as.data.frame(cpm_sig_small)
# Populate the first col with the gene names (required formatting for CIBERSORT)
cpm_sig_small <- data.frame(Gene = rownames(cpm_sig_small),
                                 cpm_sig_small,
                                 row.names = NULL)

# Export the matrix as tsv file - Important: row.names = FALSE or CIBERSORT gives an error
#write.table(cpm_sig_small, file = "./CIBERSORT/count_matrix.tsv", sep = "\t", 
#            quote = FALSE, row.names = FALSE, col.names = TRUE)

######################## Create synthetic test scRNA dataset for validation ########################
# Load necessary library
#library(devtools)
#install_github("humengying0907/deconvBenchmarking")
library(deconvBenchmarking)

# DeconvBenchmark requires a metadata file with sample names as rows specifying
# the cell type and the sample ID for each single cell sample.

# Get the complementary samples which are NOT used for the creation of the matrix (test data)
complement_cell_samples <- setdiff(1:nr_samples, random_cell_samples)
# Check (should be 0)
intersect(complement_cell_samples, random_cell_samples)

# Create complementary single cell dataset
sc_complementary <- sc_integrated_annotated[, complement_cell_samples]
# Create vector with sample names
sample_names <- colnames(sc_complementary[["RNA"]]$data)
# Extract the cell type information
celltype_complementary <- sc_complementary$celltype_fine

# Exchange T cell classification with more fine grained classification
tcell_fine <- sc_tcells_annotated$tcell_annotation
intersection <- intersect(names(tcell_fine), names(celltype_complementary))
celltype_complementary[intersection] <- tcell_fine[intersection]
celltype_complementary <- gsub(" ", "_", celltype_complementary)
# Quick check (All must be TRUE!)
length(celltype_complementary) + length(random_cell_samples) == ncol(sc_integrated_annotated)
all(names(tcell_fine %in% names(celltype_complementary)))
identical(colnames(sc_complementary), names(celltype_complementary))

# Build the metadata dataframe
scMeta <- data.frame(cell_type = celltype_complementary)
sample_ID <- substr(rownames(scMeta), 1, 4)
scMeta <- data.frame(scMeta, sample_ID)

# Get the count matrix for these samples
rm(count_matrix) # Free up RAM
scExpr <- as.matrix(sc_complementary[["RNA"]]$counts)
# Apply CPM transformation
scExpr <- apply(scExpr, 2, function(x) (x / sum(x)) * 1e6)


######################## Generate random cell fraction matrix ########################

# Set seed for reproducibility
set.seed(123)

# Define sample and cell type names
num_samples <- 20
num_cell_types <-length(unique(scMeta$cell_type))

samples <- paste0("Sample", 1:num_samples)
cell_types <- unique(scMeta$cell_type)

# Generate a random fraction matrix
simulated_frac <- matrix(runif(num_samples * num_cell_types), 
                         nrow = num_samples, ncol = num_cell_types)

# Normalize each row to sum to 1
simulated_frac <- simulated_frac / rowSums(simulated_frac)

# Assign row and column names
rownames(simulated_frac) <- samples
colnames(simulated_frac) <- cell_types

# Print first few rows
print(head(simulated_frac, 5))
all(rowSums(simulated_frac)) == 1

# Run the pseudobulk simulation
simulate_bulk <- bulkSimulator_heter(
            scExpr,
            scMeta,
            colnames_of_cellType = "cell_type",
            colnames_of_sample = "sample_ID",
            simulated_frac = simulated_frac,
            min_chunkSize = min(table(scMeta$cell_type)), #chunk size cannot be bigger than celltype with least counts
            use_chunk = "all",
            export_cellUsage = F,
            n.core = 2
)


pseudo_bulk_dataset <- as.data.frame(simulate_bulk["simulated_bulk"])
colnames(pseudo_bulk_dataset) <- paste0("Sample", 1:num_samples)
pseudo_bulk_dataset <- data.frame(Gene = rownames(pseudo_bulk_dataset), 
                                  pseudo_bulk_dataset, 
                                  row.names = NULL)

# Export the pseudo bulk dataset for CIBERSORT deconvolution
#write.table(pseudo_bulk_dataset, file = "./CIBERSORT/pseudo_bulk.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
#saveRDS(simulate_bulk, file="./CIBERSORT/simulate_bulk.rds")
