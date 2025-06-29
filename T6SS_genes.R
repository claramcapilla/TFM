## Load required libraries
library(readxl)
library(writexl)

##-------------------------------
## Load DE results from Excel
##-------------------------------
# Helper function to load all sheets for one condition
load_de_sheets <- function(file_path, num_sheets) {
  lapply(1:num_sheets, function(s) read_excel(file_path, sheet = s))
}

# Define how many sheets each DE_results file has
de_sheets_map <- list(
  DE_cond1 = 2, DE_cond2 = 4, DE_cond3 = 8, DE_cond4 = 2, DE_cond5 = 6,
  DE_cond6 = 2, DE_cond7 = 2, DE_cond8 = 2, DE_cond9 = 2, DE_cond10 = 6,
  DE_cond11 = 2, DE_cond12 = 4, DE_cond13 = 2, DE_cond14 = 2
)

# Load all DE results into a named list
de_results <- lapply(names(de_sheets_map), function(cond) {
  path <- file.path("DE_results", paste0(cond, ".xlsx"))
  load_de_sheets(path, de_sheets_map[[cond]])
})
names(de_results) <- names(de_sheets_map)

##-------------------------------
## Load annotation table
##-------------------------------
feature_table <- read.delim("GCA_feature_table.txt",header = TRUE, sep = "\t", 
                            comment.char = "#", quote = "", 
                            stringsAsFactors = FALSE)
genes <- feature_table[feature_table$gene == "gene", ]
locus_info <- genes[,c(17,15,8,9,10)]
colnames(locus_info) <- c("gene_id","name","start","end","strand")

##-------------------------------
## Define T6SS gene list
##-------------------------------
t6ss_genes <- sprintf("STM14_%04d", 313:343)
t6ss_info <- locus_info[locus_info$gene_id %in% t6ss_genes, ]

## Assign curated gene names 
t6ss_info$name = c("tssA","tssG","tssF","tssE","tagJ","tagK","clpV","tssB",
                   "tssC","STM14_0322","STM14_0323","hcp1","tae4",
                   "tai4","hcp2","tssJ","tssK","tssL","STM14_0331",
                   "sciR","tssM","tagF","tldi1","tlde1",
                   "STM14_0337","vgrG","eagR","PAAR","STM14_0341","Ntax47",
                   "STM14_0343")
write_xlsx(t6ss_info, "feature_t6ss.xlsx")

##-------------------------------
## Extract T6SS activation/repression info
##-------------------------------
# Helper function to extract matching t6ss_info for a given result table
extract_t6ss_matches <- function(t6ss_df, result_df) {
  t6ss_df[t6ss_df$gene_id %in% result_df$gene_id, ]
}

# Extract all activation/repression matches for all DE result sheets
t6ss_matches <- list()

for (cond in names(de_results)) {
  for (i in seq_along(de_results[[cond]])) {
    result_name <- paste(cond, paste0("sheet", i), sep = "_")
    t6ss_matches[[result_name]] <- extract_t6ss_matches(t6ss_info, 
                                                        de_results[[cond]][[i]])
  }
}