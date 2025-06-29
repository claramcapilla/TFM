## Load required libraries
library(readxl)
library(openxlsx)

## Read metadata Excel file
res <- read_excel("SRA - S. enterica RNA.xlsx", col_names = TRUE) 

## Remove columns that contain only NA values
datos <- res[,colSums(is.na(res))<nrow(res)]
ncol(res)-ncol(datos) ## Number or removed columns
setdiff(colnames(res),colnames(datos)) ## Names of removed columns
rm(res) ## Clean up workspace

## Filter for RNA-Seq experiments only
datos=datos[datos$`Assay Type`=="RNA-Seq",]
datos=datos[datos$LibrarySource=="TRANSCRIPTOMIC",]

## Function to summarise null values in any variable
resumen=function(var) {
  print(table(var))
  print(any(is.na(var)))
}

attach(datos)
## Select only relevant metadata columns
metadatos=as.data.frame(cbind(Run,Organism,Serovar,strain,genotype,
                              `genotype/variation`,growth_condition,growth_phase,
                              source_type,source_name,medium,Culture_condition,
                              isolate,isol_growth_condt,isolation_source,HOST,
                              Host_disease,cell_type,replicate,number_of_cells,
                              Instrument,Platform,LibraryLayout,LibrarySelection,
                              LibrarySource,`DATASTORE filetype`,treatment,
                              sample_type,temp,stimulus,antibiotic_treatment,plasmid))

## Standarise values in the 'Serovar' column
metadatos$Serovar=ifelse(metadatos$Serovar=="Typhimuriun","Typhimurium",metadatos$Serovar)
metadatos$Serovar=ifelse(metadatos$Serovar=="Enteriditis","Enteritidis",metadatos$Serovar)

## Standarise values in the 'strain' column
metadatos$strain=ifelse(metadatos$strain=="Salmonella enterica subsp. enterica serovar typhimurium ATCC 14028","ATCC 14028",metadatos$strain)
metadatos$strain=ifelse(metadatos$strain=="14028S","14028s",metadatos$strain)
metadatos$strain=ifelse(metadatos$strain=="ATCC14028s","ATCC 14028s",metadatos$strain)
metadatos$strain=ifelse(metadatos$strain=="Salmonella enterica subsp. enterica serovar Typhimurium str. ATCC 14028","ATCC 14028",metadatos$strain)
metadatos$strain=ifelse(metadatos$strain=="strain 14028s","14028s",metadatos$strain)
metadatos$strain=ifelse(metadatos$strain=="Typhimurium 14028S","14028s",metadatos$strain)
metadatos$strain=ifelse(metadatos$strain=="wild type","WT",metadatos$strain)
metadatos$strain=ifelse(metadatos$strain=="wildtype","WT",metadatos$strain)

## Export cleaned metadata
write.xlsx(metadatos, file = "metadata.xlsx")
