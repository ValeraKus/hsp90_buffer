rm(list = ls(all = T))

tissue_num <- read.table('../../Body/2_Derived/Gtex.LOE.GOE.tissues.numbers.frec.more.5.perc.txt')
gene <- gsub('.*b37_','', tissue_num$cis_eQTL_id)
tissue_num <- cbind(tissue_num, gene)

LOE_tissues <- c()
LOE_cis_eQTLs <- c()
GOE_tissues <- c()
GOE_cis_eQTLs <- c()
unique_genes <- unique(gene)
for (g in unique_genes){
  #count common amount of tissues for each gene in which there is cis-eQTLs
  loe_tis <- sum(tissue_num[tissue_num$gene == g,'LOE_tissues'])
  goe_tis <- sum(tissue_num[tissue_num$gene == g,'GOE_tissues'])
  LOE_tissues <- c(LOE_tissues,loe_tis)
  GOE_tissues <- c(GOE_tissues,goe_tis)
  #count common amount of unique cis-eQTLs for each gene
  loe_eqtl <- dim(tissue_num[(tissue_num$gene == g) & (tissue_num$LOE_tissues > 0),])[1]
  goe_eqtl <- dim(tissue_num[(tissue_num$gene == g) & (tissue_num$GOE_tissues > 0),])[1]
  LOE_cis_eQTLs <- c(LOE_cis_eQTLs, loe_eqtl)
  GOE_cis_eQTLs <- c(GOE_cis_eQTLs, goe_eqtl)
}
##usually tissues and cis_eQTLs are equal. But if the same cis-eQTL is presented in several tissues tissue number will be higher 


#saving results
df <- data.frame(cbind(unique_genes, as.numeric(as.vector(LOE_tissues)), as.numeric(as.vector(GOE_tissues)), as.numeric(as.vector(LOE_cis_eQTLs)), as.numeric(as.vector(GOE_cis_eQTLs))) )
colnames(df) <- c('gene', 'LOE_tissues', 'GOE_tissues', 'LOE_cis_eQTLs', 'GOE_cis_eQTLs')
write.table(df, '../../Body/2_Derived/Gtex.LOE.GOE.tissue.ciseQTL.numbers.per.gene.like.hsp90.txt')


#drawing plots
  pdf('../../Body/4_Figures/Gtex.LOE.GOE.tissue.ciseQTL.numbers.genes.like.hsp.pdf')
  par(mfcol =c(1,1))
  
  plot(LOE_tissues,GOE_tissues, asp = 1, xlab = 'LOE cis-eQTL sum in all tissue', ylab = 'GOE cis-eQTL sum in all tissue')
  points(LOE_tissues[which(unique_genes == 'ENSG00000096384')], GOE_tissues[which(unique_genes == 'ENSG00000096384')], col = 'red', pch = 19)
  legend(25, 40, legend = 'hsp90', col = 'red', pch=19)
  title('LOE vs GOE in all tissues')
  
  plot(LOE_cis_eQTLs, GOE_cis_eQTLs, xlab = 'unique LOE cis-eQTLs number', ylab ='unique GOE cis-eQTLs numbe')
  points(LOE_cis_eQTLs[which(unique_genes == 'ENSG00000096384')], GOE_cis_eQTLs[which(unique_genes == 'ENSG00000096384')], col = 'red', pch = 19)
  legend(25, 30, legend = 'hsp90', col = 'red', pch=19)
  title('LOE vs GOE cis-eQTL number per gene')
  
  
  #pdf('../../Body/4_Figures/Gtex.LOE.GOE.tissue.numbers.boxplot.hsp.like.genes.pdf')
  
  #in how much tissues there is some particular cis-eQTL in tach gene
  boxplot(tissue_num[tissue_num$gene == 'ENSG00000096384','GOE_tissues'], tissue_num[tissue_num$gene != 'ENSG00000096384','GOE_tissues'], names = c('hsp90', 'other genes'))
  title('GOE cis-eQTL tissue number')
  boxplot(tissue_num[tissue_num$gene == 'ENSG00000096384','LOE_tissues'], tissue_num[tissue_num$gene != 'ENSG00000096384','LOE_tissues'], names = c('hsp90', 'other genes'))
  title('LOE cis-eQTL tissue number')
  
  
  
  boxplot(LOE_cis_eQTLs, GOE_cis_eQTLs, names = c('LOE', 'GOE'))
  title('cis-eQTL number per gene')
  
  
  boxplot(GOE_cis_eQTLs/LOE_cis_eQTLs)
  points(GOE_cis_eQTLs[which(unique_genes == 'ENSG00000096384')]/LOE_cis_eQTLs[which(unique_genes == 'ENSG00000096384')], col = 'red', pch = 19)
  title('GOE/LOE cis-eQTL number')
  legend(1.2, 11, legend = 'hsp90', col = 'red', pch = 19)
  
  boxplot(GOE_tissues/LOE_tissues)
  points(GOE_tissues[which(unique_genes == 'ENSG00000096384')]/LOE_tissues[which(unique_genes == 'ENSG00000096384')], col = 'red', pch = 19)
  title('GOE/LOE tissue number\nAssesed allele frequency > 5%')
  legend('topright', legend = 'hsp90', col = 'red', pch = 19)
  dev.off()






