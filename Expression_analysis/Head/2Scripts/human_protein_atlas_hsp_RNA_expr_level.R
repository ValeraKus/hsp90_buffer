rm(list=ls(all=TRUE))

hpa <- read.table('../../Body/1_Raw/normal_tissue.tsv', sep = '\t', header = T)

hsp_hpa <- read.table("../../Body/1_Raw/ENSG00000096384.tsv", sep = '\t', header = T)

tissue <- colnames(hsp_hpa)[grepl("Tissue", colnames(hsp_hpa))]
expression_level <- as.numeric(hsp_hpa[tissue])
names(expression_level) <- NULL

tissue <- gsub("..NX.", "", gsub("Tissue.RNA...", "", tissue))
tissue <- gsub("_1", "", gsub("\\.", "_", tissue))


hsp_expr <- data.frame(tissue, expression_level)


write.table(hsp_expr, "../../Body/2_Derived/human_protein_atlas_hsp_expression_level.csv", sep = ",")

hsp90 <- hpa[hpa$Gene == "ENSG00000096384",]
write.table(hsp90, "../../Body/2_Derived/human_protein_atlas_hsp_protein_level.csv", sep = ",")


hsc70 <- read.table('../../Body/1_Raw/HSPA8.tsv', sep = '\t', header = T)
hsc70 <- hsc70[hsc70$Ensembl == 'ENSG00000109971',]
tissue <- colnames(hsc70)[grepl("Tissue", colnames(hsc70))]
expression_level <- as.numeric(hsc70[tissue])
names(expression_level) <- NULL
tissue <- gsub("..NX.", "", gsub("Tissue.RNA...", "", tissue))
tissue <- gsub("_1", "", gsub("\\.", "_", tissue))

hsp_expr <- data.frame(tissue, expression_level)
write.table(hsp_expr, "../../Body/2_Derived/human_protein_atlas_hsc70_expression_level.csv", sep = ",")
