# Load required libraries
library(readxl)
library(ggplot2)
library(gggenes)
library(dplyr)
library(scales)
library(pheatmap)
library(FactoMineR)
library(factoextra)

##------------------------------
## T6SS Gene Map and Expression Analysis
##------------------------------
## Load annotated T6SS feature table
feature_table_t6ss <- read_excel("feature_t6ss.xlsx",col_names = TRUE)

## Build gene table for T6SS cluster
genes_t6ss <- data.frame(
  gene = c("tssA                 ","tssG","tssF","tssE","tagJ","tagK","clpV","tssB","tssC",
           "STM14_0322","STM14_0323","hcp1","tae4","tai4","hcp2","tssJ",
           "tssK","tssL","STM14_0331","sciR","tssM","tagF","tldi1",
           "tlde1","STM14_0337","vgrG","eagR","PAAR","STM14_0341",
           "Ntax47","STM14_0343"),
  start = feature_table_t6ss$start,
  end = feature_table_t6ss$end,
  strand = c(rep(FALSE,6), rep(TRUE,25)),
  act = c(rep(TRUE,6),"",rep(TRUE,4),rep("",3),TRUE,rep("",4),rep(TRUE,7),
          "",rep(TRUE,4)))
## Set gene factor levels to preserve orden in plotting
genes_t6ss$gene = factor(genes_t6ss$gene, levels = genes_t6ss$gene)
## Coordinates for marking active genes
activacion <- genes_t6ss %>%
  filter(!is.na(start)) %>%
  mutate(pos = (start + end)/2)

## Plot gene map with activation markers
x_min <- min(genes_t6ss$start)
x_max <- max(genes_t6ss$end)

ggplot(genes_t6ss, aes(xmin = start, xmax = end, y = 1, fill = gene, forward = strand)) +
  geom_gene_arrow(arrow_body_height = unit(6, "mm"), 
                  arrowhead_height = unit(8, "mm")) + 
  scale_x_continuous(limits = c(x_min-60, x_max+60), expand = c(0, 0)) +
  theme_genes() +
  theme(
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.box = "horizontal",
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 12)
  ) +
  guides(fill = guide_legend(nrow = 7, byrow = F)) +
  geom_point(data = activacion %>% filter(act != ""), 
             aes(x = pos, y = 1), shape = 21, size = 2, fill = "yellow") 

##------------------------------
## Heatmap of log2FC
##------------------------------
logfc_t6ss = read_excel("logfc_t6ss_reduced.xlsx")
rownames(logfc_t6ss) <- c("tssA        ","tssG        ","   tssF        ",
                          "   tssE         ","tagJ         ","    tagK         ",
                          "    clpV         ","    tssB         ","    tssC         ",
                          "STM14_0322","STM14_0323","    hcp1         ",
                          "    tae4         ","    tai4         ","    hcp2         ",
                          "    tssJ         ","    tssK         ","    tssL         ",
                          "STM14_0331","    sciR         ","    tssM         ",
                          "    tagF         ","    tldi1        ",
                          "    tlde1         ","STM14_0337","    vgrG         ",
                          "    eagR         ",
                          "    PAAR         ","STM14_0341","    Ntax47       ",
                          "STM14_0343")
## Create heatmap
pheatmap(t(logfc_t6ss),
         cluster_rows = TRUE, 
         cluster_cols = FALSE, 
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-3, 3, length.out = 101), 
         main = "Differential expression in T6SS",
         angle_col = 90,
         fontsize_row = 10,
         fontsize_col = 9,
         border_color = NA)

##------------------------------
## Auxiliary Gene Map for Heatmap
##------------------------------
## Generate synthetic coordinates for visual map
genes_t6ss_heatmap <- data.frame(
  gene = c("STM14_0313 (tssA)","STM14_0314 (tssG)","STM14_0315 (tssF)",
           "STM14_0316 (tssE)","STM14_0317 (tagJ)","STM14_0318 (tagK)",
           "STM14_0319 (tssH)","STM14_0320 (tssB)","STM14_0321 (tssC)",
           "STM14_0322","STM14_0323","STM14_0324","STM14_0325",
           "STM14_0326","STM14_0327 (tssD)","STM14_0328 (tssJ)",
           "STM14_0329 (tssK)","STM14_0330 (tssL)","STM14_0331",
           "STM14_0332 (sciR)","STM14_0333 (tssM)","STM14_0334","STM14_0335",
           "STM14_0336","STM14_0337","STM14_0338 (tssI)","STM14_0339",
           "STM14_0340 (PAAR)","STM14_0341", "STM14_0342 (Ntax47)",
           "STM14_0343"),
  start = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140,
            150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270,
            280,290,300),
  end = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140,
            150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270,
            280, 290, 300, 310),
  strand = c(rep(FALSE,6), rep(TRUE,25)))

genes_t6ss_heatmap$gene = factor(genes_t6ss_heatmap$gene, levels = genes_t6ss_heatmap$gene)

## Plot map
ggplot(genes_t6ss_heatmap, aes(xmin = start, xmax = end, y = 1, forward = strand)) +
  geom_gene_arrow(fill = "gray90", 
                  arrow_body_height = unit(6, "mm"), 
                  arrowhead_height = unit(8, "mm")) +
  geom_gene_label(aes(label = gene), size = 2.5, angle = 45, align = "left") +
  theme_genes() +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank())

##------------------------------
## PCA of T6SS log2FC across conditions
##------------------------------
## Load reduced log2FC matrix
logfc_red = read_excel("logfc_t6ss_reduced.xlsx",col_names = TRUE)
colnames(logfc_red) = c("HeLa cells","Lettuce", "Starvation at \n low CFU",
                        "Starvation at high CFU", "pH 3", "Novobiocin",
                        "Mitomycin C", "MgCl2", "Oxidative stress",
                        "Coculture A. castellanii", "MOPS", "Mg + pH 5.8",
                        "Mg + pH 7.7", "L-arabinose", "Autoclaved \n DS soil",
                        "42ºC", "15ºC")
logfc_red=t(logfc_red)
logfc_red[is.na(logfc_red) | logfc_red == ""] <- 0

## Run PCA
pca_result <- prcomp(logfc_red, scale = FALSE)
pca_scores <- as.data.frame(pca_result$x)
pca_scores$Condition <- rownames(pca_scores)

pca_loadings <- as.data.frame(pca_result$rotation)
pca_loadings$Gene <- rownames(pca_loadings)

# Scale arrows for better visualization
scale_factor <- max(abs(pca_scores[, 1:2])) / max(abs(pca_loadings[, 1:2])) * 0.7
pca_loadings_scaled <- pca_loadings
pca_loadings_scaled[, 1:2] <- pca_loadings_scaled[, 1:2] * scale_factor

# Variance explained
total_var <- sum(pca_result$sdev^2)
var_explained <- (pca_result$sdev^2 / total_var) * 100
x_lab <- paste0("PC1 (", round(var_explained[1], 2), "%)")
y_lab <- paste0("PC2 (", round(var_explained[2], 2), "%)")

# Plot PCA
ggplot(pca_scores, aes(x = PC1, y = PC2)) +
  geom_point(size = 5, color = "steelblue") +
  geom_text(aes(label = Condition), vjust = -1, size = 3.5, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  labs(x = x_lab, y = y_lab) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
  )