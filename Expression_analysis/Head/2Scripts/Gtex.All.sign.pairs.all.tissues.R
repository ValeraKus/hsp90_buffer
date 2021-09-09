rm(list=ls(all=TRUE))


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
  All_sign_pairs = rbind(All_sign_pairs, QTL)
  }}


f <- read.table('../../Body/1_Raw/GTEx_Analysis_v7_eQTL/Adipose_Subcutaneous.v7.signif_variant_gene_pairs.txt.gz', sep = '\t', header = TRUE)

# to leave cis-eQTL that is outside of the gene
#All_sign_pairs <- All_best_var[All_best_var$gene_start > All_best_var$pos | All_best_var$pos > All_best_var$gene_end, ]

gz1 <- gzfile("../../Body/2_Derived/All.sign.pairs.all.tissues.txt.gz", "w")

write.csv(All_sign_pairs, gz1)
close(gz1)


All_sign_pairs <- read.table('../../Body/2_Derived/All.sign.pairs.all.tissues.txt.gz')

pdf("../../Body/4_Figures/All.sign.pairs.all.tissues.hsp90ab1.pdf")

plot(All_sign_pairs[grepl("ENSG00000096384", All_sign_pairs$varian_geneid_tss), "slope"], All_sign_pairs[grepl("ENSG00000096384", All_sign_pairs$varian_geneid_tss), "maf"],
     xlab = "slope", ylab = "maf", main = 'HSP90AB1')

x <- barplot(table(All_sign_pairs[grepl("ENSG00000096384", All_sign_pairs$varian_geneid_tss), "Tissue"]), xaxt="n",
        ylab = 'number of cis-eQTLs', main = "HSP90AB1")
labs <- paste(names(table(All_sign_pairs[grepl("ENSG00000096384", All_sign_pairs$varian_geneid_tss), "Tissue"])), "Tissue")
text(cex=1, x=x-.25, y=-1.25, labs, xpd=TRUE, srt=45)

dev.off()

All_sign_pairs[All]

