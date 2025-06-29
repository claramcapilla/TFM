# Comparative Analysis of Gene Expression in RNA-Seq Studies of Salmonella enterica: Focus on the Type VI Secretion System

### metadata_cleaning.R
This script processes and cleans RNA-Seq sample metadata from the NCBI SRA for Salmonella enterica. It performs the following tasks:
 - Loads raw metadata from an Excel spreadsheet.
 - Removes metadata columns with only missing values.
 - Filters for transcriptomic experiments of type RNA-Seq.
 - Selects and organises relevant experimental variables (e.g., strain, growth condition, sequencing platform).
 - Standardises inconsistent or misspelled entries in key variables such as Serovar and strain.
 - Exports a cleaned metadata table (metadatos.xlsx) for downstream analysis.
This script supports the creation of a reproducible transcriptomic dataset and facilitates filtering by experimental condition, strain, or host model.

### functions_DE.R
This script provides a collection of modular R functions designed to facilitate differential gene expression (DGE) analysis of RNA-seq datasets in Salmonella enterica. It supports workflows for both experiments with biological replicates (using DESeq2) and those without replicates (using NOISeq). The resulting outputs include tables of significantly activated and repressed genes, as well as comprehensive log2 fold-change values. Additionally, the script allows automated export of results to neatly formatted Excel files, supporting analyses involving one to four treatment conditions.
The script includes functions for:
- Performing DGE analysis using DESeq2, which is applied when two or more biological replicates are available.
- Performing DGE analysis using NOISeq, suitable for single-replicate datasets.
- Extracting gene lists for activated and repressed genes based on user-defined thresholds for fold change and statistical significance.
- Exporting results to Excel workbooks using the openxlsx package, allowing for convenient review and downstream use.

The primary inputs include:
- A gene count matrix (gene_count_matrix) containing raw read counts, with genes in rows and samples in columns.
- A design matrix (design) indicating the treatment (treat) and control (control) groups via a group column.
- Threshold parameters such as adjusted p-value (pv), log2 fold-change (lfc), or probability (pb) depending on the method used.

For each comparison, the functions return:
- Two data frames listing genes that are significantly activated or repressed.
- Optionally, a complete list of log2 fold-change values, with or without NA filtering.
- Excel files with separate worksheets for each condition, supporting comparisons with one to four treatments.
