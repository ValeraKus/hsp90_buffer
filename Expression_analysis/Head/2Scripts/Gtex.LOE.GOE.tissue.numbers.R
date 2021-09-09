rm(list=ls(all=TRUE))


All_best_var <- read.table('../../Body/2_Derived/All.best.vat.all.tissue.variants.outside.genes.important.columns.txt')
All_best_var[,'genes'] <- gsub('[0-9, A-Z]+_[0-9]+_[A-Z]+_[A-Z]+_b37_', '', All_best_var$cis_eQTL_id)
All_best_var['hsp90'] <- All_best_var$genes == 'ENSG00000096384'
All_best_var['GOE'] <- All_best_var$slope > 0

genes <- read.table('../../Body/1_Raw/gencode.v25.annotation.gtf.Genes.Shet.pLI.FIS.RVIS.GHIS.KnKs.GC.BrainSpecificRanking.Branch', header = 1)
cds <-  as.vector(genes[genes$GeneType == 'protein_coding', 'EnsemblId'])


genes.like.hsp <- as.vector(read.table('../../Body/2_Derived/hsp.like.genes.pca.genes.ranged.by.distance.to.hsp.txt')[1:300,1])
All_best_var_genes_like_hsp <- All_best_var[All_best_var$genes %in% genes.like.hsp,] #choose hsp-like genes only

All_best_var
library(stringr)

All_best_var_protein_coding <- All_best_var[cds %in% All_best_var$genes,]


only_ovary <- All_best_var_genes_like_hsp[All_best_var_genes_like_hsp$Tissue == 'Ovary',]
plot(only_ovary$Assessed_Allele_Freq, only_ovary$slope, col = 'grey')
points(only_ovary[only_ovary$hsp90, 'Assessed_Allele_Freq'], only_ovary[only_ovary$hsp90, 'slope'], col = 'red', cex = 1.5, pch = 19)
legend('topright', legend = 'HSP90AB1', col = 'red', pch=19)

pdf('../../Body/4_Figures/Gtex.LOE.GOE.tissue.specific.slopes.assessed.allele.freq.pdf')
for (tissue in unique(All_best_var_protein_coding$Tissue)){
  only_tis <- All_best_var_protein_coding[All_best_var_protein_coding$Tissue == tissue,]
  plot(only_tis$Assessed_Allele_Freq, only_tis$slope, col = 'grey')
  points(only_tis[only_tis$hsp90, 'Assessed_Allele_Freq'], only_tis[only_tis$hsp90, 'slope'], col = 'red', cex = 1.5, pch = 19)
  legend('topright', legend = 'HSP90AB1', col = 'red', pch=19)
  title(sprintf('Assessed allele frequency vs slope \n of cis-eQTLs in %s', tissue))
}
dev.off()

pdf('../../Body/4_Figures/Gtex.LOE.GOE.tissue.specific.slopes.assessed.allele.freq.hsp.like.genes.pdf')
for (tissue in unique(All_best_var_genes_like_hsp$Tissue)){
  only_tis <- All_best_var_genes_like_hsp[All_best_var_genes_like_hsp$Tissue == tissue,]
  plot(only_tis$Assessed_Allele_Freq, only_tis$slope, col = 'grey')
  points(only_tis[only_tis$hsp90, 'Assessed_Allele_Freq'], only_tis[only_tis$hsp90, 'slope'], col = 'red', cex = 1.5, pch = 19)
  legend('topright', legend = 'HSP90AB1', col = 'red', pch=19)
  title(sprintf('Assessed allele frequency vs slope \n of cis-eQTLs in %s', tissue))
}
dev.off()

pdf('../../Body/4_Figures/Gtex.LOE.GOE.tissue.numbers.slopes.assessed.allele.freq.pdf')
#for all genes




boxplot(All_best_var$Assessed_Allele_Freq ~ All_best_var$hsp90, names = c('all genes', 'HSP90AB1'))
title('Assesed Allele Frequency')


boxplot(abs(All_best_var[All_best_var$slope < 50 & All_best_var$slope > -50, 'slope']) ~ 
          All_best_var[All_best_var$slope < 50 & All_best_var$slope > -50, 'hsp90'] * 
          All_best_var[All_best_var$slope < 50 & All_best_var$slope > -50, 'GOE_or_LOE'],
        names = c('all genes LOE', 'HSP90AB1 LOE', 'all genes GOE', 'HSP90AB1 GOE'),
        main = 'Slopes of cis-eQTLs')


boxplot(All_best_var[All_best_var$slope < 50 & All_best_var$slope > -50, 'Assessed_Allele_Freq'] ~ 
          All_best_var[All_best_var$slope < 50 & All_best_var$slope > -50, 'hsp90'] * 
          All_best_var[All_best_var$slope < 50 & All_best_var$slope > -50, 'GOE_or_LOE'],
        names = c('all genes LOE', 'HSP90AB1 LOE', 'all genes GOE', 'HSP90AB1 GOE'),
        main = 'Assesed Allele Frequency of cis-eQTLs')


#for genes like hsp90 from PCA


plot(All_best_var_genes_like_hsp$Assessed_Allele_Freq, All_best_var_genes_like_hsp$slope, col = 'grey')
points(All_best_var_genes_like_hsp[All_best_var_genes_like_hsp$hsp90, 'Assessed_Allele_Freq'], 
       All_best_var_genes_like_hsp[All_best_var_genes_like_hsp$hsp90, 'slope'], col = 'red', pch = 19, cex=1.5)
legend('topright', legend = 'HSP90AB1', col = 'red', pch=19)


plot(log10(All_best_var_genes_like_hsp$Assessed_Allele_Freq), All_best_var_genes_like_hsp$slope, col = 'grey')
points(log10(All_best_var_genes_like_hsp[All_best_var_genes_like_hsp$hsp90, 'Assessed_Allele_Freq']), 
       All_best_var_genes_like_hsp[All_best_var_genes_like_hsp$hsp90, 'slope'], col = 'red', pch = 19, cex=1.5)
legend('topright', legend = 'HSP90AB1', col = 'red', pch=19)


boxplot(All_best_var_genes_like_hsp$Assessed_Allele_Freq ~ All_best_var_genes_like_hsp$hsp90, names = c('hsp90 like genes', 'HSP90AB1'))
title('Assesed Allele Frequency in genes like hsp90')


boxplot(abs(All_best_var_genes_like_hsp[, 'slope']) ~ 
          All_best_var_genes_like_hsp[, 'hsp90'] * 
          All_best_var_genes_like_hsp[, 'GOE_or_LOE'],
        names = c('hsp90 like genes LOE', 'HSP90AB1 LOE', 'hsp90 like genes GOE', 'HSP90AB1 GOE'),
        main = 'Slopes of cis-eQTLs in genes like hsp90')


boxplot(All_best_var_genes_like_hsp[, 'Assessed_Allele_Freq'] ~ 
          All_best_var_genes_like_hsp[, 'hsp90'] * 
          All_best_var_genes_like_hsp[, 'GOE_or_LOE'],
        names = c('hsp90 like genes LOE', 'HSP90AB1 LOE', 'hsp90 like genes GOE', 'HSP90AB1 GOE'),
        main = 'Assesed Allele Frequency of cis-eQTLs in genes like hsp90')



#All_best_var_genes_like_hsp_freq <- (All_best_var_genes_like_hsp[All_best_var_genes_like_hsp$Assessed_Allele_Freq >= 0.05,])
#look at hsp90
All_best_var_hsp <- All_best_var_genes_like_hsp[grepl('ENSG00000096384', All_best_var_genes_like_hsp$cis_eQTL_id),]
length(All_best_var_hsp$cis_eQTL_id) == length(unique(All_best_var_hsp$cis_eQTL_id)) #TRUE 47
length(All_best_var_hsp$Tissue) == length(unique(All_best_var_hsp$Tissue)) #TRUE 47
# it means that hsp90 has one unique cis-eQTL in one unique tissue

dim(All_best_var_genes_like_hsp_freq[grepl('ENSG00000096384', All_best_var_genes_like_hsp_freq$cis_eQTL_id),]) #17

hsp_LOE <- All_best_var_hsp[All_best_var_hsp$slope < 0, 'Assessed_Allele_Freq']
hsp_GOE <- All_best_var_hsp[All_best_var_hsp$slope > 0, 'Assessed_Allele_Freq']

x <- c(hsp_GOE, hsp_LOE)
gp = c(rep("hsp90 GOE Assessed allele freq", length(hsp_GOE)),rep("hsp90 LOE Assessed allele freq",length(hsp_LOE)))
boxplot(x ~ gp, xlab = '', ylab = 'Assesed allele frequency', mail = 'Assessed Allele Frequency HSP90AB1')


t.test(hsp_GOE, hsp_LOE) #no difference p-value = 0.7528
wilcox.test(hsp_GOE, hsp_LOE) #p-value = 0.2658
LOE <- All_best_var_genes_like_hsp[All_best_var_genes_like_hsp$slope < 0, 'Assessed_Allele_Freq']
GOE <- All_best_var_genes_like_hsp[All_best_var_genes_like_hsp$slope > 0, 'Assessed_Allele_Freq']
boxplot(LOE, GOE, names=c('LOE', 'GOE'), ylab = 'Assesed allele frequency')
title('cis-eQTL assesed allele frequency')

dev.off()

# count tissues in which each unique cis-eQTL is occured
cis_eQTL_id <- as.array(unique(All_best_var_genes_like_hsp$cis_eQTL_id))
LOE_ciseQTLs <- c()
GOE_ciseQTLs <- c()

#the loop takes a while
for (QTL in cis_eQTL_id){
  LOE_count <- 0
  LOE <- All_best_var_genes_like_hsp[(All_best_var_genes_like_hsp$cis_eQTL_id == QTL) & (All_best_var_genes_like_hsp$slope < 0),]$Tissue
  LOE_count <- length(unique(LOE))
  LOE_ciseQTLs <- c(LOE_ciseQTLs, LOE_count)
  GOE_count <- 0  
  GOE <- All_best_var_genes_like_hsp[(All_best_var_genes_like_hsp$cis_eQTL_id == QTL) & (All_best_var_genes_like_hsp$slope > 0),]$Tissue
  GOE_count <- length(unique(GOE))
  GOE_ciseQTLs <- c(GOE_ciseQTLs, GOE_count)
}

#saving results
cis_eQTL_id <- as.vector(cis_eQTL_id)
tissue_numbers <- data.frame(cbind(cis_eQTL_id, LOE_ciseQTLs, GOE_ciseQTLs))
write.table(tissue_numbers, '../../Body/2_Derived/Gtex.LOE.GOE.tissues.numbers.txt')





# count tissues in which each unique cis-eQTL is occured
cis_eQTL_id <- as.array(unique(All_best_var_genes_like_hsp$cis_eQTL_id))
LOE_tissues <- c()
GOE_tissues <- c()

#the loop takes a while
for (QTL in cis_eQTL_id){
  LOE_count <- 0
  LOE <- All_best_var_genes_like_hsp_freq[(All_best_var_genes_like_hsp_freq$cis_eQTL_id == QTL) & (All_best_var_genes_like_hsp_freq$slope < 0),]$Tissue
  LOE_count <- length(unique(LOE))
  LOE_tissues <- c(LOE_tissues, LOE_count)
  GOE_count <- 0  
  GOE <- All_best_var_genes_like_hsp_freq[(All_best_var_genes_like_hsp_freq$cis_eQTL_id == QTL) & (All_best_var_genes_like_hsp_freq$slope > 0),]$Tissue
  GOE_count <- length(unique(GOE))
  GOE_tissues <- c(GOE_tissues, GOE_count)
}

#saving results
cis_eQTL_id <- as.vector(cis_eQTL_id)
tissue_numbers <- data.frame(cbind(cis_eQTL_id, LOE_tissues, GOE_tissues))
boxplot(as.numeric(as.vector(tissue_numbers$LOE_tissues)),as.numeric(as.vector(tissue_numbers$GOE_tissues)))

write.table(tissue_numbers, '../../Body/2_Derived/Gtex.LOE.GOE.tissues.numbers.frec.more.5.perc.txt')
