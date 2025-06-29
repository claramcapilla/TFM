## Load dependencias
source("functions_DE.R") ## custom DE functions
library(writexl)
library(readxl)

set.seed(123)

## Helper to load gene count matrix
load_gcm <- function(path) {
  df <- read.csv(path, header = TRUE, sep = ",")
  rownames(df) <- df[[1]]
  return(df[,-1])
}

## -------------------------------
## CONDITION 1: HeLa cells
## -------------------------------
gcm_1 = load_gcm("gene_count_matrix/gcm_1.csv")
design_1 = data.frame(sample=colnames(gcm_1),group=c(rep("control",3),rep("treat",3)))
DE_1 = DE_treat_vs_wt(gcm_1, design_1, 0.05, 2)
save_excel_1cond(DE_1$act, DE_1$rep, "DE_results/DE_cond1.xlsx")

## -------------------------------
## CONDITION 2: Lettuce & root exudate
## -------------------------------
gcm_2 = load_gcm("gene_count_matrix/gcm_2.csv")
design_2_treat1 = data.frame(sample = colnames(gcm_2[, 1:6]),
                             group = c(rep("control", 3), rep("treat", 3)))
design_2_treat2 = data.frame(sample = colnames(gcm_2[, c(1,2,3,7,8,9)]),
                             group = c(rep("control", 3), rep("treat", 3)))
DE_2_treat1 = DE_treat_vs_wt(gcm_2[, 1:6], design_2_treat1, 0.05, 2)
DE_2_treat2 = DE_treat_vs_wt(gcm_2[, c(1,2,3,7,8,9)], design_2_treat2, 0.05, 2)
save_excel_2cond(DE_2_treat1$act, DE_2_treat1$rep,
                 DE_2_treat2$act, DE_2_treat2$rep,
                 "DE_results/DE_cond2.xlsx")

## -------------------------------
## CONDITION 1: Starvation
## -------------------------------
gcm_3 = load_gcm("gene_count_matrix/gcm_3.csv")
design_3_treat1 = data.frame(sample = colnames(gcm_3[, c(1,2,7,8)]),
                             group = c(rep("control", 2), rep("treat", 2)))
design_3_treat2 = data.frame(sample = colnames(gcm_3[, c(1,2,9,10)]),
                             group = c(rep("control", 2), rep("treat", 2)))
design_3_treat3 = data.frame(sample = colnames(gcm_3[, 1:4]),
                             group = c(rep("control", 2), rep("treat", 2)))
design_3_treat4 = data.frame(sample = colnames(gcm_3[, c(1,2,5,6)]),
                             group = c(rep("control", 2), rep("treat", 2)))
DE_3_treat1 = DE_treat_vs_wt(gcm_3[, c(1,2,7,8)], design_3_treat1, 0.05, 2)
DE_3_treat2 = DE_treat_vs_wt(gcm_3[, c(1,2,9,10)], design_3_treat2, 0.05, 2)
DE_3_treat3 = DE_treat_vs_wt(gcm_3[, 1:4], design_3_treat3, 0.05, 2)
DE_3_treat4 = DE_treat_vs_wt(gcm_3[, c(1,2,5,6)], design_3_treat4, 0.05, 2)
save_excel_4cond(DE_3_treat1$act, DE_3_treat1$rep,
                 DE_3_treat2$act, DE_3_treat2$rep,
                 DE_3_treat3$act, DE_3_treat3$rep,
                 DE_3_treat4$act, DE_3_treat4$rep,
                 "DE_results/DE_cond3.xlsx")

## -------------------------------
## CONDITION 1: pH 3
## -------------------------------
gcm_4 = load_gcm("gene_count_matrix/gcm_4.csv")
DE_4 = NOISeq_treat_vs_wt(gcm_4, 0.9, 1)
save_excel_1cond(DE_4$act, DE_4$rep, "DE_results/DE_cond4.xlsx")

## -------------------------------
## CONDITION 5: Novobiocin
## -------------------------------
gcm_5 = load_gcm("gene_count_matrix/gcm_5.csv")
design_5_treat1 = data.frame(sample = colnames(gcm_5[, c(1:3, 7:9)]),
                             group = c(rep("control", 3), rep("treat", 3)))
design_5_treat2 = data.frame(sample = colnames(gcm_5[, c(1:3, 10:12)]),
                             group = c(rep("control", 3), rep("treat", 3)))
design_5_treat3 = data.frame(sample = colnames(gcm_5[, c(1:3, 13:15)]),
                             group = c(rep("control", 3), rep("treat", 3)))
DE_5_treat1 = DE_treat_vs_wt(gcm_5[, c(1:3, 7:9)], design_5_treat1, 0.05, 2)
DE_5_treat2 = DE_treat_vs_wt(gcm_5[, c(1:3, 10:12)], design_5_treat2, 0.05, 2)
DE_5_treat3 = DE_treat_vs_wt(gcm_5[, c(1:3, 13:15)], design_5_treat3, 0.05, 2)
save_excel_3cond(DE_5_treat1$act, DE_5_treat1$rep,
                 DE_5_treat2$act, DE_5_treat2$rep,
                 DE_5_treat3$act, DE_5_treat3$rep,
                 "DE_results/DE_cond5.xlsx")

## -------------------------------
## CONDITION 6: Mitomycin C
## -------------------------------
gcm_6 = load_gcm("gene_count_matrix/gcm_6.csv")
design_6 = data.frame(sample=colnames(gcm_6),group=c(rep("control",3),rep("treat",3)))
DE_6 = DE_treat_vs_wt(gcm_6, design_6, 0.05, 2)
save_excel_1cond(DE_6$act, DE_6$rep, "DE_results/DE_cond6.xlsx")

## -------------------------------
## CONDITION 7: MgCl2
## -------------------------------
gcm_7= load_gcm("gene_count_matrix/gcm_7.csv")
design_7 = data.frame(sample=colnames(gcm_7),group=c(rep("control",2),rep("treat",2)))
DE_7 = DE_treat_vs_wt(gcm_7, design_7, 0.05, 2)
save_excel_1cond(DE_7$act, DE_7$rep, "DE_results/DE_cond7.xlsx")

## -------------------------------
## CONDITION 8: Oxidative stress
## -------------------------------
gcm_8 = load_gcm("gene_count_matrix/gcm_8.csv")
design_8 = data.frame(sample=colnames(gcm_8),group=c(rep("control", 4),rep("treat", 4)))
DE_8 = DE_treat_vs_wt(gcm_8, design_8, 0.05, 2)
save_excel_1cond(DE_8$act, DE_8$rep, "DE_results/DE_cond8.xlsx")

## -------------------------------
## CONDITION 9: Coculture with A. castellanii
## -------------------------------
gcm_9 = load_gcm("gene_count_matrix/gcm_9.csv")
design_9 = data.frame(sample=colnames(gcm_9),group=c(rep("control", 3),rep("treat", 3)))
DE_9 = DE_treat_vs_wt(gcm_9, design_9, 0.05, 2)
save_excel_1cond(DE_9$act, DE_9$rep, "DE_results/DE_cond9.xlsx")

## -------------------------------
## CONDITION 10: MOPS, Mg, pH 5.8 and 7.7
## -------------------------------
gcm_10 = load_gcm("gene_count_matrix/gcm_10.csv")
DE_10_treat1 = NOISeq_treat_vs_wt(gcm_10[, 1:2], 0.9, 1)
DE_10_treat2 = NOISeq_treat_vs_wt(gcm_10[, c(1,3)], 0.9, 1)
DE_10_treat3 = NOISeq_treat_vs_wt(gcm_10[, c(1,4)], 0.9, 1)
save_excel_3cond(DE_10_treat1$act, DE_10_treat1$rep,
                 DE_10_treat2$act, DE_10_treat2$rep,
                 DE_10_treat3$act, DE_10_treat3$rep,
                 "DE_results/DE_cond10.xlsx")

## -------------------------------
## CONDITION 11: L-arabinose
## -------------------------------
gcm_11 = load_gcm("gene_count_matrix/gcm_11.csv")
design_11 = data.frame(sample=colnames(gcm_11),group=c(rep("control",2),rep("treat", 2)))
DE_11 = DE_treat_vs_wt(gcm_11, design_11, 0.05, 2)
save_excel_1cond(DE_11$act, DE_11$rep, "DE_results/DE_cond11.xlsx")

## -------------------------------
## CONDITION 12: DS soil
## -------------------------------
gcm_12 = load_gcm("gene_count_matrix/gcm_12.csv")
design_12_treat1 = data.frame(sample = colnames(gcm_12[, 1:6]),
                             group = c(rep("control", 3), rep("treat", 3)))
design_12_treat2 = data.frame(sample = colnames(gcm_12[, c(1,2,3,7,8,9)]),
                             group = c(rep("control", 3), rep("treat", 3)))
DE_12_treat1 = DE_treat_vs_wt(gcm_12[, 1:6], design_12_treat1, 0.05, 2)
DE_12_treat2 = DE_treat_vs_wt(gcm_12[, c(1,2,3,7,8,9)], design_12_treat2, 0.05, 2)
save_excel_2cond(DE_12_treat1$act, DE_12_treat1$rep,
                 DE_12_treat2$act, DE_12_treat2$rep,
                 "DE_results/DE_cond12.xlsx")

## -------------------------------
## CONDITION 13: 42ºC
## -------------------------------
gcm_13 = load_gcm("gene_count_matrix/gcm_13.csv")
design_13 = data.frame(sample=colnames(gcm_13),group=c(rep("control", 2), rep("treat", 2)))
DE_13 = DE_treat_vs_wt(gcm_13, design_13, 0.05, 2)
save_excel_1cond(DE_13$act, DE_13$rep, "DE_results/DE_cond13.xlsx")

## -------------------------------
## CONDITION 14: 15ºC
## -------------------------------
gcm_14 = load_gcm("gene_count_matrix/gcm_14.csv")
design_14 = data.frame(sample=colnames(gcm_14),group=c(rep("control", 4), rep("treat", 4)))
DE_14 = DE_treat_vs_wt(gcm_14, design_14, 0.05, 2)
save_excel_1cond(DE_14$act, DE_14$rep, "DE_results/DE_cond14.xlsx")

## -----Export DE logfc CSVs (NA removed)-----
export_logc <- function(func, gcm, design, out_csv) {
  res <- func(gcm, design)
  write.csv(res, out_csv, row.names = FALSE)
}
export_logfc(DE_logfc, gcm_1, design_1, "DE_result_csv/DE_cond1.csv")
## Repeat for condition 1 to condition 14

## -----Export DE logfc CSVs (NA included)-----
export_logfc(DE_logfc_NA, gcm_1, design_1, "DE_result_csv2/DE_cond1.csv")
## Repeat for condition 1 to condition 14

## -----Export NOISeq logFC CSVs-----
cond4 <- NOISeq_logfc(load_gcm("gene_count_matrix/gcm.csv"))
write.csv(cond4, "DE_results_csv/DE_cond4.csv", row.names = FALSE)
cond4_na <- NOISeq_logfc_NA(load_gcm("gene_count_matrix/gcm.csv"))
write.csv(cond4_na, "DE_results_csv2/DE_cond4.csv", row.names = FALSE)
## etc

## -----Build full logFC tables across all conditions-----
logfc_list_all <- list(
  "HeLa cells" = DE_1_NA$log2FC,
  "Lettuce" = DE_2_t1_NA$log2FC,
  "Root exudate" = DE_2_t2_NA$log2FC,
  "Starvation low CFU 24h" = DE_3_t1_NA$log2FC,
  "Starvation low CFU 3d" = DE_3_t2_NA$log2FC,
  "Starvation high CFU 4h" = DE_3_t3_NA$log2FC,
  "Starvation high CFU 4d" = DE_3_t4_NA$log2FC,
  "pH 3" = DE_4_NA$log2FC,
  "Novobiocin 10 min" = DE_5_t1_NA$log2FC,
  "Novobiocin 20 min" = DE_5_t2_NA$log2FC,
  "Novobiocin 60 min" = DE_5_t3_NA$log2FC,
  "Mitomycin C" = DE_6_NA$log2FC,
  "MgCl2" = DE_7_NA$log2FC,
  "Oxidative stress" = DE_8_NA$log2FC,
  "Coculture A. castellanii" = DE_9_NA$log2FC,
  "MOPS" = DE_10_t1_NA$log2FC,
  "Mg + pH 5.8" = DE_10_t2_NA$log2FC,
  "Mg + pH 7.7" = DE_10_t3_NA$log2FC,
  "L-arabinose" = DE_11_NA$log2FC,
  "DS soil" = DE_12_t1_NA$log2FC,
  "Autoclaved DS soil" = DE_12_t2_NA$log2FC,
  "42ºC" = DE_13_NA$log2FC,
  "15ºC" = DE_14_NA$log2FC
)

table_logfc_all <- as.data.frame(logfc_list_all)
rownames(table_logfc_all) <- load_gcm("gene_count_matrix/gcm_1.csv") |> rownames()

write_xlsx(table_logfc_all, "logfc_all.xlsx")

##----Build reduced log2FC table for selected conditions-----
selected_conditions <- c(
  "HeLa cells", "Lettuce", "Starvation low CFU 3d", "Starvation high CFU 4d",
  "pH 3", "Novobiocin 60 min", "Mitomycin C", "MgCl2", "Oxidative stress",
  "Coculture A. castellanii", "MOPS", "Mg + pH 5.8", "Mg + pH 7.7",
  "L-arabinose", "Autoclaved DS soil", "42ºC", "15ºC"
)

table_logfc_reduced <- table_logfc_all[, selected_conditions, drop = FALSE]
write_xlsx(table_logfc_reduced, "logfc_reduced.xlsx")

##-----Build T6SS gene-specific tables-----
t6ss_genes <- sprintf("STM14_%04d", 313:343)

extract_t6ss <- function(tbl, genes) {
  subset(tbl, rownames(tbl) %in% genes)
}
logfc_t6ss_all <- extract_t6ss(table_logfc_all, t6ss_genes)
write_xlsx(logfc_t6ss_all, "logfc_t6ss_all.xlsx")

logfc_t6ss_reduced <- extract_t6ss(table_logfc_reduced, t6ss_genes)
write_xlsx(logfc_t6ss_reduced, "logfc_t6ss_reduced.xlsx")

##-----Write gene lists to ClueGO input-----
export_cluego_list <- function(excel_file, sheet, out_txt) {
  df <- read_excel(excel_file, sheet = sheet, col_names = TRUE)
  writeLines(df$gene_id, out_txt)
}
export_cluego_list("DE_results/DE_cond1.xlsx", 1, "clueGO/cond1_expr.txt")
## Repeat for condition 2 to condition 14
