rm(list=ls(all=TRUE))

Genes = read.table("../../Body/1_Raw/gencode.v25.annotation.gtf.Genes.Shet.pLI.FIS.RVIS.GHIS.KnKs.GC.BrainSpecificRanking.Branch", header = TRUE)
names(Genes)
Genes = Genes[,c(1,15)]
Genes = Genes[!is.na(Genes$GenomeWideHaploinsufficiencyScore),]
length(unique(Genes$EnsemblId))
Genes = unique(Genes)

List = list.files("../../Body/1_Raw/GTEx_Analysis_v7_eQTL/")
All_sign_pairs = c()
for (i in 1:length(List))
{ # i = 2
  infile = paste('../../Body/1_Raw/GTEx_Analysis_v7_eQTL/',List[i], sep = '')
  if (grepl('v7.signif_variant_gene_pairs.txt.gz',infile))
  {tissue = gsub('.v7.signif_variant_gene_pairs.txt.gz','',List[i])
  QTL = read.table(infile, sep = '\t', header = TRUE)
  QTL$gene_id = gsub("\\.(.*)",'',QTL$gene_id)
  QTL$varian_geneid_tss <- paste(QTL$variant_id, QTL$gene_id, QTL$tss_distance, sep = "_")
  QTL$Tissue = tissue
  QTL <- QTL[, c('varian_geneid_tss', "maf", "pval_nominal", "slope", "slope_se", "Tissue")]
  ExtractEnsemblId <- function(x){ unlist(strsplit(x,'_'))[6]}
  QTL$EnsemblId = apply(as.matrix(QTL$varian_geneid_tss),1,FUN = ExtractEnsemblId)
  QTL = QTL[QTL$EnsemblId %in% Genes$EnsemblId,]
  All_sign_pairs = rbind(All_sign_pairs, QTL)
  }
}


Agg1 = aggregate(list(All_sign_pairs$slope,All_sign_pairs$maf,All_sign_pairs$pval_nominal),by = list(All_sign_pairs$EnsemblId), FUN = median)
names(Agg1)=c('EnsemblId','MedianSlope','MedianMaf','MedianP')

All_sign_pairs$one = 1
Agg2 = aggregate(All_sign_pairs$one, by = list(All_sign_pairs$EnsemblId), FUN = sum)
names(Agg2)=c('EnsemblId','NumberEqtl')

Agg3 = aggregate(All_sign_pairs[All_sign_pairs$slope < 0,]$one, by = list(All_sign_pairs[All_sign_pairs$slope < 0,]$EnsemblId), FUN = sum)
names(Agg3)=c('EnsemblId','NumberLoeEqtl')

Agg4 = aggregate(All_sign_pairs[All_sign_pairs$slope > 0,]$one, by = list(All_sign_pairs[All_sign_pairs$slope > 0,]$EnsemblId), FUN = sum)
names(Agg4)=c('EnsemblId','NumberGoeEqtl')

Agg=merge(Agg1,Agg2)
Agg=merge(Agg,Agg3, all = TRUE)
Agg=merge(Agg,Agg4, all = TRUE)
length(unique(Agg$EnsemblId))
Agg=merge(Agg,Genes,all.x = TRUE)
length(unique(Agg$EnsemblId))
Agg$MedianP = -log10(Agg$MedianP)
summary(Agg$MedianP)

write.table(Agg,"../../Body/2_Derived/Gtex.All.sign.pairs.all.tissues.Haploins.txt")

