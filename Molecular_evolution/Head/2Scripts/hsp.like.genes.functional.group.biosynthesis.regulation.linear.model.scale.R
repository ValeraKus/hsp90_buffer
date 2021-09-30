rm(list=ls(all=TRUE))
df <- read.csv('../../Body/2_Derived/enrichmentByNetworks.csv', header = T)
df <- as.data.frame(apply(df, 2, sub, pattern = ' ', replacement = ''))



group1 <- df$GO.0009889.regulation.of.biosynthetic.process


ids <- read.csv('../../Body/2_Derived/ensids.csv', header = T)

group1_ids <- as.vector(ids[grepl(paste(group1, collapse="|"), ids$Gene.Symbol), 'Ensembl.Gene.ID'])



knks <- read.table('../../Body/2_Derived/kn.ks.genes.like.hsp.txt')
knks[knks == 'n/a'] <- NA
knks[, -1] <- sapply(knks[, -1], as.vector)
knks[, -1] <- sapply(knks[, -1], as.numeric)

knks_gr1 <- knks[,grepl(paste(group1_ids, collapse = '|'), colnames(knks))]
knks_gr1$Species <- knks$Species
knks_gr1$Generation_Length <- knks$Generation_Length

hsp90 <- na.omit(knks_gr1[,c('dN.dS_ENSG00000096384', 'Species', 'Generation_Length')])
df_hsp90_sp <- knks_gr1[knks_gr1$Species %in% hsp90$Species,] #to make linear model for other genes using more or less same species as hsp90




cols <- colnames(df_hsp90_sp[,-c(1, 302)])
for (c in cols) {
  if (length(df_hsp90_sp[is.na(df_hsp90_sp[,c]),c]) >= 15){
    df_hsp90_sp[,c] <- NULL
  }
  
  
}



for (i in 1:(ncol(df_hsp90_sp)-2))
{# i = 71
  Temp=na.omit(df_hsp90_sp[,c(i,92)])
  Nrow=nrow(Temp)
  if (Nrow > 10)
    
  {
    Gene=as.character(names(Temp[1]))
    
    # Lm1
    Lm1<-lm(scale(Temp[,1])~scale(Temp$Generation_Length))
    res = data.frame(summary(Lm1)[4])
    Lm1Intercept = res[1,1]; Lm1InterceptP = res[1,4]
    Lm1Coeff = res[2,1]; Lm1CoeffP = res[2,4]
    
    # Lm2
    Lm2<-lm(scale(Temp[,1])~ 0 + scale(Temp$Generation_Length))
    res = data.frame(summary(Lm2)[4])
    Lm2Coeff = res[1,1]; Lm2CoeffP = res[1,4]
    
    # var
    CoeffOfVar = sd(Temp[,1])/mean(Temp[,1])
    
    # RankCor
    RankCor <-cor.test(Temp[,1],Temp$Generation_Length, method='spearman')
    RankCorP = RankCor$p.value
    RankCorRho = RankCor$estimate
    
    ResLine=data.frame(Gene,Nrow,CoeffOfVar,Lm1Intercept,Lm1InterceptP,Lm1Coeff,Lm1CoeffP,Lm2Coeff,Lm2CoeffP,RankCorP,RankCorRho)
    if (i==2) {FinalRes = ResLine}
    if (i>2) {FinalRes = rbind(FinalRes,ResLine)}
  }
}

pdf('../../Body/4_Figures/Regulation.of.biosynthesis.process.genes.lin.reg.scale.pdf')
#plotting LM1 results
plot(FinalRes$Lm1Coeff,-log10(FinalRes$Lm1CoeffP), xlab = 'Lm1Coeff', ylab = '-log10(Lm1CoeffP)'); 
points(FinalRes[FinalRes$Gene == "dN.dS_ENSG00000096384",]$Lm1Coeff,-log10(FinalRes[FinalRes$Gene == "dN.dS_ENSG00000096384",]$Lm1CoeffP), pch = 16, col = 'red', cex = 2);
abline(h=-log10(0.05),col = 'red')
legend('topleft', legend = 'HSP90AB1', col = 'red', pch = 19)
title('lm(scale(dN.dS) ~ scale(Generation_length)), n = 89')


plot(FinalRes$Lm1Intercept,FinalRes$Lm1InterceptP, xlab = 'Lm1Intercept', ylab = 'Lm1InterceptP'); 
points(FinalRes[FinalRes$Gene == "dN.dS_ENSG00000096384",]$Lm1Intercept,FinalRes[FinalRes$Gene == "dN.dS_ENSG00000096384",]$Lm1InterceptP, col = 'red', pch = 16, cex = 2);
legend('topleft', legend = 'HSP90AB1', col = 'red', pch = 19)
title('lm(scale(dN.dS) ~ scale(Generation_length)), n = 89')


plot(FinalRes$Lm1Intercept,FinalRes$Lm1Coeff, xlab = 'Lm1Intercept', ylab = 'Lm1Coeff'); par(new=TRUE)
points(FinalRes[FinalRes$Gene == "dN.dS_ENSG00000096384",]$Lm1Intercept,FinalRes[FinalRes$Gene == "dN.dS_ENSG00000096384",]$Lm1Coeff, col = 'red', pch = 16,  cex = 2);
legend('topleft', legend = 'HSP90AB1', col = 'red', pch = 19)
title('lm(scale(dN.dS) ~ scale(Generation_length)), n = 89')


## Lm2 is our favourite model till now (through the origin with scaled X)
plot(FinalRes$Lm2Coeff,-log10(FinalRes$Lm2CoeffP), xlab = 'Lm2Coeff', ylab = '-log10(Lm2CoeffP)'); 
points(FinalRes[FinalRes$Gene == "dN.dS_ENSG00000096384",]$Lm2Coeff,-log10(FinalRes[FinalRes$Gene == "dN.dS_ENSG00000096384",]$Lm1CoeffP), col = 'red', pch = 16, cex = 2);
abline(h=-log10(0.05),col = 'red')
legend('topleft', legend = 'HSP90AB1', col = 'red', pch = 19)
title('lm(scale(dN.dS) ~ 0 + scale(Generation_length)), n = 89')


plot(FinalRes$RankCorRho, -log10(FinalRes$RankCorP))
points(FinalRes[FinalRes$Gene == "dN.dS_ENSG00000096384",]$RankCorRho,-log10(FinalRes[FinalRes$Gene == "dN.dS_ENSG00000096384",]$RankCorP), col = 'red', pch = 16, cex = 2)
abline(h=-log10(0.05),col = 'red')
legend('topleft', legend = 'HSP90AB1', col = 'red', pch = 19)
title('Spearman rank correlation, n = 89')
dev.off()

