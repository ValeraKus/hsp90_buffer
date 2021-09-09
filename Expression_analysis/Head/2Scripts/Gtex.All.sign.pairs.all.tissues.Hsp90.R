rm(list=ls(all=TRUE))
library(ggplot2)

pdf("../../Body/4_Figures/GTEx.All.sign.pairs.all.tissues.Hsp90")
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
  QTL = QTL[grepl("ENSG00000096384",QTL$varian_geneid_tss),]
  All_sign_pairs = rbind(All_sign_pairs, QTL)
  }
}

write.table(All_sign_pairs,"../../Body/2_Derived/GTEx.All.sign.pairs.all.tissues.Hsp90.txt")

All_sign_pairs <- read.table('../../Body/2_Derived/GTEx.All.sign.pairs.all.tissues.Hsp90.txt')


plot(All_sign_pairs$slope,All_sign_pairs$maf)
abline(v=0,col='red')
summary(All_sign_pairs$slope)
wilcox.test(All_sign_pairs$slope, mu = 0) # hehe!!

theme_set(
  theme_minimal() +
    theme(legend.position = "right", axis.text.x = element_text(angle = 0))
)
ggplot(data = All_sign_pairs, aes(slope, maf, colour = Tissue))+
  geom_point()+
  ggtitle('cis-eQTLs of HSP90AB1')

par(pin = c(6,3))


#counts <- table(All_sign_pairs$slope >0, All_sign_pairs$Tissue)
#barplot(counts, main="LOE and GOE cis-eQTLs of HSP90AB1",
#        col=c("darkblue","red"), xlab = "Tissue", ylab = 'Number of cis-eQTLs',
#        legend = c("LOE", "GOE"), beside=TRUE)

All_sign_pairs$GOE <- All_sign_pairs$slope > 0
theme_set(
  theme_minimal() +
    theme(legend.position = "right", axis.text.x = element_text(angle = -20, hjust = 0))
)

ggplot(data=All_sign_pairs, aes(x=Tissue, fill= GOE)) +
  geom_histogram(stat="count", position=position_dodge())+
  ggtitle("HSP90AB1")

for (tissue in unique(All_sign_pairs$Tissue)){
  plot(All_sign_pairs[All_sign_pairs$Tissue == tissue,]$slope,All_sign_pairs[All_sign_pairs$Tissue == tissue,]$maf,
       main = paste('cis-eQTLs of HSP90AB1 in ', tissue))
  
}


dev.off()







