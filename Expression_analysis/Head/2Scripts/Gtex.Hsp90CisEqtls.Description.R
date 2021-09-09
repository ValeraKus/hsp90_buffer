rm(list=ls(all=TRUE))

List = list.files("../../Body/1_Raw/GTEx_Analysis_v7_eQTL/")
# v7.egenes.txt.gz; *.egenes.txt.gz files contain data for all genes tested; to obtain the list of eGenes, select the rows with 'qval' ≤ 0.05.
# v7.signif_variant_gene_pairs.txt.gz; all significant variant-gene associations based on permutations
AllSignifVarGenePairs = c()
TheBestVarPerGene = c()
for (i in 1:length(List))
{ # i = 2
  infile = paste('../../Body/1_Raw/GTEx_Analysis_v7_eQTL/',List[i], sep = '')
  if (grepl('v7.signif_variant_gene_pairs.txt.gz',infile))
  {
  tissue = gsub('.v7.signif_variant_gene_pairs.txt.gz','',List[i])  
  QTL = read.table(infile, sep = '\t', header = TRUE)
  QTL$gene_id = gsub("\\.(.*)",'',QTL$gene_id)
  QTL = QTL[QTL$gene_id == 'ENSG00000096384',]
  if (nrow(QTL) > 0) {QTL$Tissue = tissue; AllSignifVarGenePairs = rbind(AllSignifVarGenePairs,QTL)}
  }
  if (grepl('v7.egenes.txt.gz',infile))
   {
    tissue = gsub('.v7.egenes.txt.gz','',List[i])  
    QTL = read.table(infile, sep = '\t', header = TRUE)
    QTL$gene_id = gsub("\\.(.*)",'',QTL$gene_id)
    QTL = QTL[QTL$gene_id == 'ENSG00000096384',]
    if (nrow(QTL) > 0) {QTL$Tissue = tissue; TheBestVarPerGene = rbind(TheBestVarPerGene,QTL)}
   }
}

write.table(TheBestVarPerGene,file = "../../Body/2_Derived/Gtex.Hsp90.TheBestVarPerGene.txt")
write.table(AllSignifVarGenePairs,file = "../../Body/2_Derived/Gtex.Hsp90.AllSignifVarGenePairs.txt")


plot(AllSignifVarGenePairs$slope,AllSignifVarGenePairs$maf)                      # which allele is minor allele here (assessed or not?)  
plot(TheBestVarPerGene$slope*TheBestVarPerGene$ref_factor,TheBestVarPerGene$maf) # minor allele effect versus minor allele frequency . Is it correct?

