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

### create_count_matrix.R
This R script is designed to generate structured gene count matrices from multiple .tabular files output by tools such as featureCounts (e.g., in Galaxy workflows). Each .tabular file contains two columns: gene identifiers and raw read counts. The script performs the following steps:
- Reads a series of raw count files (control and treatment replicates) corresponding to different experimental conditions.
- Assembles each condition’s counts into a unified matrix, assigning appropriate sample names.
- Stores the matrices as .csv files for downstream differential expression analysis using packages such as DESeq2 or NOISeq.

### functions_DE.R
This script provides a collection of modular R functions designed to facilitate differential gene expression (DGE) analysis of RNA-seq datasets. It supports workflows for both experiments with biological replicates (using DESeq2) and those without replicates (using NOISeq). The resulting outputs include tables of significantly activated and repressed genes, as well as comprehensive log2 fold-change values. Additionally, the script allows automated export of results to neatly formatted Excel files, supporting analyses involving one to four treatment conditions.
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

### DE_tables.R
This R script processes RNA-seq gene count data from multiple experimental conditions and executes differential gene expression (DGE) analyses using DESeq2 (for replicated designs) and NOISeq (for single-replicate datasets). The script:
- Loads gene count matrices exported from Galaxy workflows.
- Defines experimental designs for 14 distinct conditions (e.g., HeLa cells, lettuce, starvation, pH stress, antibiotics, temperature shifts).
- Runs DGE analysis using custom wrapper functions.
- Exports results to structured Excel files, with activated/repressed gene lists per condition.
- Compiles complete and reduced log2 fold-change (log2FC) matrices across all genes and selected conditions.
- Extracts log2FC data for genes in the Type VI Secretion System (T6SS) operon.
- Prepares gene lists for use in ClueGO (e.g., for enrichment analysis in Cytoscape).

### T6SS_genes.R
This R script processes differential expression (DE) results obtained under 14 environmental and stress-related conditions to identify genes involved in the Type VI Secretion System (T6SS). What the script does:
- Loads DE results from multiple Excel files, each representing a condition with one or more comparisons (replicates or treatments).
- Loads and parses a GCA feature table to extract genomic annotation, including gene IDs, names, strand, and coordinates.
- Defines a curated list of T6SS-associated genes, mapped by locus_tag and manually annotated with functional names.
- Extracts activation and repression events for T6SS genes across all conditions, based on their presence in DE result tables.
- Stores results in a structured list (t6ss_matches) for further analysis or export (e.g., visualization or functional enrichment).

### representation.R
This script generates visual representations and performs principal component analysis (PCA) on the expression of genes within the T6SS operon of Salmonella Typhimurium 14028S. It builds on previously computed log₂ fold-change (log₂FC) values from RNA-seq differential expression analysis under diverse environmental and stress conditions. Key functionalities:
- Gene Map Plotting: Visualizes the genomic organization of the T6SS operon using gggenes, highlighting activated genes
- Heatmap Generation: Displays a heatmap of log₂FC values for all T6SS genes across selected experimental conditions.
- Gene Map for Heatmap Reference: Creates a static map of all T6SS genes for consistent reference alongside the heatmap.
- PCA Analysis: Performs PCA on log₂FC data to identify global expression trends of the T6SS operon across conditions. Genes and conditions are represented in PC space for intuitive biological interpretation.

### salmonella_T6SS.sql
This script contains the complete definition and data of the salmonella_T6SS database, developed as part of the project focused on gene expression analysis in the T6SS of Salmonella. Script contents:
- Instructions to create and configure the database (CREATE DATABASE, ALTER DATABASE)
- Definition of all relational tables (CREATE TABLE)
- Setup of primary and foreign key constraints (PRIMARY KEY, FOREIGN KEY)
- Insertion of all records into each table (INSERT INTO ... VALUES (...))
All relationships are normalized and ensure referential integrity. How to use:
- Open SQL Server Management Studio (SSMS).
- Create a new query window and load the salmonella_T6SS.sql file.
- Execute the full script (F5).
The salmonella_T6SS database will be created with all data ready for querying.

