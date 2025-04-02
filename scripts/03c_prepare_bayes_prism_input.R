######################## Setting up workspace ########################
library(dplyr)

# Bayesprism needs 4 input files to work:
# RAW count matrix of bulkRNA seq data
# RAW single cell count matrix
# Cell type labels containing cell type info for each single cell
# Cell state labels containing fine grained cell type info (subclusters)
# The filtering and synthetic dataset creation is done in the same way as for CIBERSORT
# So the two methods can be compared!

setwd("C:/Users/Julian/OneDrive - Universitaet Bern/PM_Project/DSP_analysis/data")
setwd("C:/Users/julia/OneDrive - Universitaet Bern/PM_Project/DSP_analysis/data") #home

# Load the datasets
sc_integrated_annotated <- readRDS("./single_cell/sc_integrated_annotated.rds")
sc_tcells_annotated <- readRDS("./single_cell/tcells_annotated.rds")

######################## Setup raw count matrix for prism ########################
# Filter dataset for treatment naive samples only and exclude pleural effusion samples
sc_integrated_annotated <- subset(sc_integrated_annotated, 
                                  subset = treatment == "naive" & sample_type == "tumor")
sc_tcells_annotated <- subset(sc_tcells_annotated, 
                              subset = treatment == "naive" & sample_type == "tumor")
# Check results
table(sc_integrated_annotated$treatment)

# Merge fine celltype information
celltype_fine <- sc_integrated_annotated$celltype_fine
t_celltype_fine <- sc_tcells_annotated$tcell_annotation
intersection <- intersect(names(celltype_fine), names(t_celltype_fine))
# Create copy of old celltypes and merge fine cell type info
sc_integrated_annotated$celltype_fine_old <- celltype_fine
sc_integrated_annotated$celltype_fine[intersection] <- t_celltype_fine[intersection]

# Sample randomly from the dataset and create a subset with less samples
set.seed(123)
nr_samples <- ncol(sc_integrated_annotated)
random_cell_samples <- sample(1:nr_samples, (nr_samples * 0.5), replace = FALSE) # Keep 50% of the cell samples
sc_subset <- sc_integrated_annotated[, random_cell_samples]
# Check the cell distrubtion in the subsample
table(sc_subset$celltype_fine)

# Manual adjustments: Subset only contains 1 mast cells and 4 pericytes!
# Want to include more of them in the signature matrix!
mast_cells <- which(sc_integrated_annotated$celltype_fine == "Mast cells")
pericytes <- which(sc_integrated_annotated$celltype_fine == "Pericytes")
random_cell_samples <- unique(c(random_cell_samples, pericytes[1:4], mast_cells))
# Redefine the subset with the added cell samples
sc_subset <- sc_integrated_annotated[, random_cell_samples]
table(sc_subset$celltype_fine)

# Extract the count matrix with raw counts matrix
count_matrix_small <- as.matrix(sc_subset[["RNA"]]$counts)

## Export the raw count matrix as tsv file for Bayes Pris
write.table(count_matrix_small, file = "./PRISM/count_matrix_PRISM_raw.tsv", sep = "\t", 
            quote = FALSE, row.names = TRUE, col.names = TRUE)

######################## Create the cell type annotation files ########################

# Use the old cell type definition as the broad cell type definition
# And the fine grained definition containing immune cells as cell states
celltype_info_prism <- data.frame(sc_subset$celltype_fine_old, 
                                  sc_subset$celltype_fine)
names(celltype_info_prism) <- c("cell.type", "cell.state")

# Check presence of ambiguous cell state mappings 
celltype_info_prism %>%
  distinct(cell.state, cell.type) %>%
  group_by(cell.state) %>%
  summarise(n_types = n()) %>%
  filter(n_types > 1)

# Some cell states map to multiple cell types!! (NOT ALLOWED BY BAYES PRISM)
# Fine grained classification does not always agree with broad classification!!!
# Since these disagreements are rare will use fine classification where they do not agree

# Define mapping table for correct classification
state_type_map <- c(
  "Cytotoxic T cells CD8+" = "T cells",
  "Exhausted T cells CD8+" = "T cells",
  "Naive T cells" = "T cells",
  "NK cells" = "NK cells"
)

# Fix cell.type wherever cell.state is known to belong to a specific type
celltype_info_prism <- celltype_info_prism %>%
  mutate(cell.type = ifelse(
    cell.state %in% names(state_type_map) & cell.type != state_type_map[cell.state],
    state_type_map[cell.state],
    cell.type
  ))

# Change other to T_other and NK_other to keep details
celltype_info_prism <- celltype_info_prism %>%
  mutate(cell.state = ifelse(cell.state == "other",
                             paste0(gsub(" cells", "", cell.type), "_other"),
                             cell.state))

# Check again if no ambiguous states exist anymore
celltype_info_prism %>%
  distinct(cell.state, cell.type) %>%
  group_by(cell.state) %>%
  summarise(n_types = n()) %>%
  filter(n_types > 1)

# Export the dataframe
write.table(celltype_info_prism, "./PRISM/cell_type_info.tsv")

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
#colnames(pseudo_bulk_dataset) <- paste0("Sample", 1:num_samples)
#pseudo_bulk_dataset <- data.frame(Gene = rownames(pseudo_bulk_dataset), 
#                                  pseudo_bulk_dataset, 
#                                  row.names = NULL)

# Export the pseudo bulk dataset for CIBERSORT deconvolution
write.table(pseudo_bulk_dataset, file = "./PRISM/pseudo_bulk_PRISM_raw.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
saveRDS(simulate_bulk, file="./PRISM/simulate_bulk.rds")
