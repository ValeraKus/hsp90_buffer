rm(list = ls(all=T))
setwd('~/Desktop/Work/hsp90_buffer/Human.genetic.variation/head/2scripts/')

#BiocManager::install("rtracklayer")
library(rtracklayer)

pgs.hsp90 <- read.table('../../body/2derived/pgs.phenotypes.hsp90.association.summary.txt')

gtex.hsp90 <- read.table('../../body/2derived/Gtex.signif.variant.gene.pairs.hsp90.txt')
gtex.hsp90$pos <- gsub('chr6_|_[A-Z].*', '', gtex.hsp90$variant_id)

phenotypes.pgs <- pgs.hsp90$phenotypeID

pgs.hsp90$cis_eQTLs_total <- 0
pgs.hsp90$LOE <- 0
pgs.hsp90$GOE <- 0
pgs.hsp90$meanSlope <- 0
pgs.hsp90$meanMAF <- 0
pgs.hsp90$LOE_meanSlope <- 0
pgs.hsp90$GOE_meanSlope <- 0
pgs.hsp90$LOE_meanMAF <- 0
pgs.hsp90$GOE_meanMAF <- 0
pgs.hsp90$cis_eQTLs_positive_association <- 0
pgs.hsp90$LOE_positive_association <- 0




for (phen in phenotypes.pgs){
  
  file <- paste('../../body/2derived/hsp90.PGS.phenotypes/ENSG00000096384.', phen, '.txt.snps.txt', sep = '')
  pgs.vars <- read.table(file, header = F, sep = '\t')
  
  #####change coordinated of pgs variants from hg19 to hg38
  chain <- import.chain("../../body/1raw/hg19ToHg38.over.chain")
  x <- pgs.vars[c(2,3)]
  x$V2 <- paste('chr', x$V2, sep='')
  x$V4 = x$V3
  colnames(x) <- c('chromosome', 'start', 'end')
  x <- GRanges(x)
  y <- liftOver(x, chain)
  pgs.vars.38 <- data.frame(y)
  pgs.vars$pos38 <- pgs.vars.38$start
  
  intersection <- merge(pgs.vars, gtex.hsp90, by.y = 'pos', by.x='pos38')
  
  pgs.hsp90[pgs.hsp90$phenotypeID == phen,]$cis_eQTLs_total <- nrow(intersection)
  pgs.hsp90[pgs.hsp90$phenotypeID == phen,]$LOE <- nrow(intersection[intersection$slope < 0,])
  pgs.hsp90[pgs.hsp90$phenotypeID == phen,]$GOE<- nrow(intersection[intersection$slope > 0,])
  pgs.hsp90[pgs.hsp90$phenotypeID == phen,]$meanSlope <- mean(intersection$slope)
  pgs.hsp90[pgs.hsp90$phenotypeID == phen,]$meanMAF <- mean(intersection$maf)
  pgs.hsp90[pgs.hsp90$phenotypeID == phen,]$LOE_meanSlope <- mean(intersection[intersection$slope < 0,'slope'])
  pgs.hsp90[pgs.hsp90$phenotypeID == phen,]$GOE_meanSlope <- mean(intersection[intersection$slope > 0,'slope'])
  pgs.hsp90[pgs.hsp90$phenotypeID == phen,]$LOE_meanMAF <- mean(intersection[intersection$slope < 0,'maf'])
  pgs.hsp90[pgs.hsp90$phenotypeID == phen,]$GOE_meanMAF <- mean(intersection[intersection$slope > 0,'maf'])
  pgs.hsp90[pgs.hsp90$phenotypeID == phen,]$cis_eQTLs_positive_association <- nrow(intersection[intersection$V6 > 0,])
  pgs.hsp90[pgs.hsp90$phenotypeID == phen,]$LOE_positive_association <- nrow(intersection[(intersection$V6 > 0) & (intersection$slope > 0),])
  
  
}


phenotypes_decode <- read.csv('../../body/1raw/pgs_scores_data.csv')
phenotypes_decode <- phenotypes_decode[,c(1,3)]
phenotypes_decode$Polygenic.Score.ID...Name <- gsub('\\(.*\\)', '', phenotypes_decode$Polygenic.Score.ID...Name)
colnames(phenotypes_decode) <- c('phenotypeID', 'Trait')

pgs.hsp90 <- merge(pgs.hsp90, phenotypes_decode, by = 'phenotypeID', all.x = T)
pgs.hsp90 <- pgs.hsp90[order(pgs.hsp90$num_of_snps, decreasing = T),]
write.table(pgs.hsp90, '../../body/2derived/hsp90.pgs.GTEx.intersection.summary.table.txt')


library(ggplot2)

ggplot(na.omit(pgs.hsp90), aes(num_of_snps, cis_eQTLs_total))+
  geom_point()


ggplot(na.omit(pgs.hsp90), aes(num_of_snps, num_of_pos_snsps))+
  geom_point()

ggplot(na.omit(pgs.hsp90), aes(x = num_of_snps, y = proportion_of_pos_snps))+
  geom_bar(stat="identity")+
  ylim(0,1)


pdf('../../body/4figures/hsp90.pgs.GTEx.intersection.summary.table.pdf', width = 14)

colors <- c("All PGS catalog SNPs" = "gray75", "Positively associated with phenotype" = "darkolivegreen2", "GTEx cis-eQTLs" = "mediumpurple2", 'GTEx cis-eQTLs positively associated with phenotype' = 'midnightblue')

ggplot(na.omit(pgs.hsp90), aes(y = num_of_snps, x = reorder(phenotypeID, num_of_snps), fill = 'All PGS catalog SNPs'))+
  geom_bar(stat='identity')+
  geom_bar(aes(x = reorder(phenotypeID, num_of_snps), y= num_of_pos_snsps, fill = 'Positively associated with phenotype'), alpha = 0.8, stat="identity")+
  geom_bar(aes(x = reorder(phenotypeID, num_of_snps), y= cis_eQTLs_total, fill = 'GTEx cis-eQTLs'), alpha = 0.8, stat="identity")+
  geom_bar(aes(x = reorder(phenotypeID, num_of_snps), y= cis_eQTLs_positive_association, fill = 'GTEx cis-eQTLs positively associated with phenotype'), alpha = 0.8, stat="identity")+
  theme_bw()+
  labs(x = "Phenotype",
       y = "number of SNPs",
       color = "Legend") +
  scale_fill_manual(values = colors, name ='')+
  theme(legend.position = 'top')
dev.off()


