## Load required libraries
library(DESeq2) ## Differential expression analysis with replicates
library(NOISeq) ## Differential expression analysis without replicates
library(openxlsx) ## Writing Excel files

## --------------------------------
## DESeq2 function (with replicates)
## Returns activated and repressed genes based on padj and log2FC cutoffs
## --------------------------------
DE_treat_vs_wt <- function(gene_count_matrix, design, pv, lfc) {
  dds <- DESeqDataSetFromMatrix(countData = gene_count_matrix, 
                                colData = design, design = ~ group)
  dds <- DESeq(dds) 
  
  res <- results(dds, contrast = c("group", "treat", "control")) 
  res <- na.omit(res)
  
  gene_ids_DE <- rownames(res) 
  log_FC <- res$log2FoldChange 
  adjusted_pval <- res$padj
  names(log_FC) <- gene_ids_DE
  names(adjusted_pval) <- gene_ids_DE
  
  activated_genes <- gene_ids_DE[log_FC > lfc & adjusted_pval < pv]
  df_activated <- data.frame(gene_id = activated_genes,
                             padj = adjusted_pval[activated_genes],
                             log2FC = log_FC[activated_genes],
                             row.names = NULL)
  
  repressed_genes <- gene_ids_DE[log_FC < -lfc & adjusted_pval < pv]
  df_repressed <- data.frame(gene_id = repressed_genes,
                             padj = adjusted_pval[repressed_genes],
                             log2FC = log_FC[repressed_genes],
                             row.names = NULL)
  
  return(list(act = df_activated, rep = df_repressed))
}


## --------------------------------
## NOISeq function (no replicates)
## Returns activated and repressed genes based on probability and log2FC cutoffs
## --------------------------------
NOISeq_treat_vs_wt <- function(gene_count_matrix, pb, lfc) {
  conds <- data.frame(condition = c("control", "treat"))
  data <- readData(data = gene_count_matrix, factors = conds)
  res <- noiseq(data, factor = "condition",norm = "tmm", 
                replicates = "no", k = 0.5) 
  table = na.omit(res@results[[1]])
  gene_id <- rownames(table)
  table = cbind(gene_id, table)
  
  df_activated = table[table$prob > pb & table$M > lfc, ][, c("gene_id", "prob", "M")]
  
  df_repressed = table[table$prob > pb & table$M < -lfc, ][, c("gene_id", "prob", "M")]
  
  return(list(act = df_activated, rep = df_repressed))
}

## --------------------------------
## DESeq2 function returning log2FC values (NA removed)
## --------------------------------
DE_logfc <- function(gene_count_matrix, design) {
  dds <- DESeqDataSetFromMatrix(countData = gene_count_matrix, 
                                colData = design, design = ~group)
  dds <- DESeq(dds) 
  
  res <- results(dds, contrast = c("group", "treat", "control")) 
  res <- na.omit(res)
  
  gene_ids_DE <- rownames(res) 
  log_FC <- res$log2FoldChange 
  
  names(log_FC) <- gene_ids_DE
  
  df_logfc <- data.frame(gene_id = gene_ids_DE,
                         log2FC = log_FC,
                         row.names = NULL)
  
  return(df_logfc)
}

## --------------------------------
## NOISeq function returning log2FC and probability (NA removed)
## --------------------------------
NOISeq_logfc <- function(gene_count_matrix) {
  conds <- data.frame(condition = c("control", "treat"))
  data <- readData(data = gene_count_matrix, factors = conds)
  res <- noiseq(data, factor = "condition",
                norm = "tmm", 
                replicates = "no", 
                k = 0.5) 
  table = na.omit(res@results[[1]])
  gene_id <- rownames(table)
  table = cbind(gene_id, table)
  
  return(table)
}

## --------------------------------
## DESeq2 function returning log2FC values (includes NA)
## --------------------------------
DE_logfc_NA <- function(gene_count_matrix, design) {
  dds <- DESeqDataSetFromMatrix(countData = gene_count_matrix, 
                                colData = design, design = ~ group)
  dds <- DESeq(dds) 
  
  res <- results(dds, contrast = c("group", "treat", "control")) 
  
  gene_ids_DE <- rownames(res) 
  log_FC <- res$log2FoldChange 
  
  names(log_FC) <- gene_ids_DE
  
  df_logfc <- data.frame(gene_id = gene_ids_DE,
                         log2FC = log_FC,
                         row.names = NULL)
  
  return(df_logfc)
}

## --------------------------------
## NOISeq function returning all results (includes NA)
## --------------------------------
NOISeq_logfc_NA <- function(gene_count_matrix) {
  conds <- data.frame(condition = c("control", "treat"))
  data <- readData(data = gene_count_matrix, factors = conds)
  res <- noiseq(data, factor = "condition",
                norm = "tmm", 
                replicates = "no", 
                k = 0.5) 
  table = res@results[[1]]
  gene_id <- rownames(table)
  table = cbind(gene_id, table)
  
  return(table)
}

## --------------------------------
## Save DE results (activated & repressed) to Excel file
## One experimental condition
## --------------------------------
save_excel_1cond <- function(activated, repressed, filename) {
  wb <- createWorkbook()
  
  addWorksheet(wb, "Activated genes")
  writeData(wb, sheet = "Activated genes", activated)
  
  addWorksheet(wb, "Repressed genes")
  writeData(wb, sheet = "Repressed genes", repressed)
  
  saveWorkbook(wb, file = filename, overwrite = TRUE)
}

## --------------------------------
## Save DE results for 2 conditions
## --------------------------------
save_excel_2cond <- function(activated1, repressed1, activated2, repressed2, filename) {
  wb <- createWorkbook()
  
  addWorksheet(wb, "Activated genes - Treat 1")
  writeData(wb, sheet = "Activated genes - Treat 1", activated1)
  
  addWorksheet(wb, "Repressed genes - Treat 1")
  writeData(wb, sheet = "Repressed genes - Treat 1", repressed1)
  
  addWorksheet(wb, "Activated genes - Treat 2")
  writeData(wb, sheet = "Activated genes - Treat 2", activated2)
  
  addWorksheet(wb, "Repressed genes - Treat 2")
  writeData(wb, sheet = "Repressed genes - Treat 2", repressed2)
  
  saveWorkbook(wb, file = filename, overwrite = TRUE)
}

## --------------------------------
## Save DE results for 3 conditions
## --------------------------------
save_excel_3cond <- function(activated1, repressed1, 
                             activated2, repressed2, 
                             activated3, repressed3, filename) {
  wb <- createWorkbook()
  
  addWorksheet(wb, "Activated genes - Treat 1")
  writeData(wb, sheet = "Activated genes - Treat 1", activated1)
  
  addWorksheet(wb, "Repressed genes - Treat 1")
  writeData(wb, sheet = "Repressed genes - Treat 1", repressed1)
  
  addWorksheet(wb, "Activated genes - Treat 2")
  writeData(wb, sheet = "Activated genes - Treat 2", activated2)
  
  addWorksheet(wb, "Repressed genes - Treat 2")
  writeData(wb, sheet = "Repressed genes - Treat 2", repressed2)
  
  addWorksheet(wb, "Activated genes - Treat 3")
  writeData(wb, sheet = "Activated genes - Treat 3", activated3)
  
  addWorksheet(wb, "Repressed genes - Treat 3")
  writeData(wb, sheet = "Repressed genes - Treat 3", repressed3)
  
  saveWorkbook(wb, file = filename, overwrite = TRUE)
}

## --------------------------------
## Save DE results for 4 conditions
## --------------------------------
save_excel_4cond <- function(activated1, repressed1, 
                             activated2, repressed2, 
                             activated3, repressed3, 
                             activated4, repressed4, filename) {
  wb <- createWorkbook()
  
  addWorksheet(wb, "Activated genes - Treat 1")
  writeData(wb, sheet = "Activated genes - Treat 1", activated1)
  
  addWorksheet(wb, "Repressed genes - Treat 1")
  writeData(wb, sheet = "Repressed genes - Treat 1", repressed1)
  
  addWorksheet(wb, "Activated genes - Treat 2")
  writeData(wb, sheet = "Activated genes - Treat 2", activated2)
  
  addWorksheet(wb, "Repressed genes - Treat 2")
  writeData(wb, sheet = "Repressed genes - Treat 2", repressed2)
  
  addWorksheet(wb, "Activated genes - Treat 3")
  writeData(wb, sheet = "Activated genes - Treat 3", activated3)
  
  addWorksheet(wb, "Repressed genes - Treat 3")
  writeData(wb, sheet = "Repressed genes - Treat 3", repressed3)
  
  addWorksheet(wb, "Activated genes - Treat 4")
  writeData(wb, sheet = "Activated genes - Treat 4", activated4)
  
  addWorksheet(wb, "Repressed genes - Treat 4")
  writeData(wb, sheet = "Repressed genes - Treat 4", repressed4)
  
  saveWorkbook(wb, file = filename, overwrite = TRUE)
}

