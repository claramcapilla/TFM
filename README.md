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
