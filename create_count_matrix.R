## -------------------------------------
## Function to construct a gene count matrix from multiple .tabular files
## Each file is expected to have gene identifiers in column 1 and counts 
## in column 2
## -------------------------------------
create_gene_count_matrix <- function(files) {
  ## Read the first file
  matrix <- read.table(files[1], header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  
  genes <- matrix[,1] ## Extract gene names
  counts <- matrix[,2, drop = FALSE] ## Extract counts as a data frame
  
  ## Loop through remaining files and append counts
  for (i in 2:length(files)) {
    table <- read.table(files[i], header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    counts <- cbind(counts, table[,2])
  }
  
  # Asign sample names (file names without extensions)
  colnames(counts) <- tools::file_path_sans_ext(basename(files))
  rownames(counts) <- genes
  counts = as.data.frame(counts)[-1,]
  
  return(as.data.frame(counts))
}

## -------------------------------------
## Define paths to counts files for each condition
## -------------------------------------
## (Conditions 1 to 14; includes control and treatment samples for each experiment)
## Example for condition 1 (HeLa cell infection)
condition_1 = c("counts_control/control_1.tabular",
                "counts_control/control_2.tabular",
                "counts_control/control_3.tabular",
                "counts_control/COND_1/ERR6412208.tabular",
                "counts_control/COND_1/ERR6412209.tabular",
                "counts_control/COND_1/ERR6412210.tabular")
## Define all other conditions anagously
## (condition 2 to condition 14 omitted for brevity)
                

## -------------------------------------
## Build gene counts matrices for each condition
## -------------------------------------
gene_count_matrix_1 = create_gene_count_matrix(condition_1)
## Repeat for condition 2 to condition 14

## -------------------------------------
## Rename columns for clarity: control and treatment replicates
## -------------------------------------
colnames(gene_count_matrix_1) = c("control_1", "control_2", "control_3",
                                  "hela_1", "hela_2", "hela_3")
## Repeat for gene_count_matrix_2 to gene_count_matrix_14 accordingly

## -------------------------------------
## Save each count matrix to CSV
## -------------------------------------
write.csv(gene_count_matrix_1, file = "gene_count_matrix/gcm_1.csv", row.names = TRUE)
## Repear for gcm_2 to gcm_14