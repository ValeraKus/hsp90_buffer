rm(list=ls(all=TRUE))



List = list.files("../../Body/1_Raw/GTEx_Analysis_v7_eQTL/")
All_best_var = c()
for (i in 1:length(List))
{ # i = 2
  infile = paste('../../Body/1_Raw/GTEx_Analysis_v7_eQTL/',List[i], sep = '')
  if (grepl('v7.egenes.txt.gz',infile))
  {tissue = gsub('.v7.egenes.txt.gz','',List[i])
  QTL = read.table(infile, sep = '\t', header = TRUE)
  QTL$gene_id = gsub("\\.(.*)",'',QTL$gene_id)
  QTL$Tissue = tissue
  All_best_var = rbind(All_best_var,QTL)
  }}


# to leave cis-eQTL that is outside of the gene
All_best_var <- All_best_var[All_best_var$gene_start > All_best_var$pos | All_best_var$pos > All_best_var$gene_end, ]

write.table(All_best_var, file = '../../Body/2_Derived/All.best.vat.all.tissue.variants.outside.genes.txt')
All_best_var <- read.table('../../Body/2_Derived/All.best.vat.all.tissue.variants.outside.genes.txt')

All_best_var$cis_eQTL_id <- paste(All_best_var$variant_id, All_best_var$gene_id, sep='_')

All_best_var$Assessed_Allele_Freq <- rep(0, dim(All_best_var)[1])  
All_best_var[All_best_var$ref_factor == 1,]$Assessed_Allele_Freq <- All_best_var[All_best_var$ref_factor == 1,]$maf
All_best_var[All_best_var$ref_factor == -1,]$Assessed_Allele_Freq <- 1 - All_best_var[All_best_var$ref_factor == -1,]$maf
All_best_var <- All_best_var[,c('cis_eQTL_id', 'Assessed_Allele_Freq', 'slope', 'Tissue')]
write.table(All_best_var, file = '../../Body/2_Derived/All.best.vat.all.tissue.variants.outside.genes.important.columns.txt')
