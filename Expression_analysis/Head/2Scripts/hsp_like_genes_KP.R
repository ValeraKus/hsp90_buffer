rm(list=ls(all=TRUE))

pdf('../../Body/4_Figures/hsp_like_genes_KP.R.pdf')

df <- read.table('../../Body/2_Derived/kn.ks.genes.like.hsp.txt')
df[df == 'n/a'] <- NA
names(df)
str(df)
df[, -1] <- sapply(df[, -1], as.vector)
df[, -1] <- sapply(df[, -1], as.numeric)
str(df)

############## 1:  RUN NAIVE ANALYSES WITH ALL AVAILABLE DATA

for (i in 2:(ncol(df)-1))
{# i = 71
  Temp=df[,c(i,302)]
  Temp=Temp[!is.na(Temp[,1]),]
  Nrow=nrow(Temp)
  if (Nrow > 10)
  
  {
  Gene=as.character(names(Temp[1]))
  
  # Lm1
  Lm1<-lm(Temp[,1]~Temp$Generation_Length)
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

plot(FinalRes$Lm1Coeff,-log10(FinalRes$Lm1CoeffP), xlab = 'Lm1Coeff', ylab = 'Lm1CoeffP'); 
points(FinalRes[FinalRes$Gene == "dN.dS_ENSG00000096384",]$Lm1Coeff,-log10(FinalRes[FinalRes$Gene == "dN.dS_ENSG00000096384",]$Lm1CoeffP), pch = 16, col = 'red', cex = 2);

plot(FinalRes$Lm1Intercept,-log10(FinalRes$Lm1InterceptP), xlab = 'Lm1Intercept', ylab = 'Lm1InterceptP'); 
points(FinalRes[FinalRes$Gene == "dN.dS_ENSG00000096384",]$Lm1Intercept,-log10(FinalRes[FinalRes$Gene == "dN.dS_ENSG00000096384",]$Lm1InterceptP), col = 'red', pch = 16, cex = 2);

plot(FinalRes$Lm1Intercept,FinalRes$Lm1Coeff, xlab = 'Lm1Intercept', ylab = 'Lm1Coeff'); par(new=TRUE)
points(FinalRes[FinalRes$Gene == "dN.dS_ENSG00000096384",]$Lm1Intercept,FinalRes[FinalRes$Gene == "dN.dS_ENSG00000096384",]$Lm1Coeff, col = 'red', pch = 16,  cex = 2);

plot(FinalRes[FinalRes$Lm1InterceptP < 0.01 & FinalRes$Lm1CoeffP < 0.01,]$Lm1Intercept,FinalRes[FinalRes$Lm1InterceptP < 0.01 & FinalRes$Lm1CoeffP < 0.01,]$Lm1Coeff, xlab = 'Lm1Intercept', ylab = 'Lm1Coeff'); par(new=TRUE)
points(FinalRes[FinalRes$Gene == "dN.dS_ENSG00000096384",]$Lm1Intercept,FinalRes[FinalRes$Gene == "dN.dS_ENSG00000096384",]$Lm1Coeff, col = 'red', pch = 16,  cex = 2);

## Lm2 is our favourite model till now (through the origin with scaled X)
plot(FinalRes$Lm2Coeff,-log10(FinalRes$Lm2CoeffP), xlab = 'Lm2Coeff', ylab = '-log10(Lm2CoeffP)'); 
points(FinalRes[FinalRes$Gene == "dN.dS_ENSG00000096384",]$Lm2Coeff,-log10(FinalRes[FinalRes$Gene == "dN.dS_ENSG00000096384",]$Lm1CoeffP), col = 'red', pch = 16, cex = 2);
abline(h=-log10(0.05),col = 'red')

### what about Nrow, CoeffOfVar and Pvalues from Lm2
cor.test(FinalRes$Nrow,FinalRes$Lm2Coeff, method = 'spearman') # nothing
cor.test(FinalRes$Nrow,FinalRes$Lm2CoeffP, method = 'spearman') # negative significant
boxplot(FinalRes[FinalRes$Lm2CoeffP>=0.05,]$Nrow,FinalRes[FinalRes$Lm2CoeffP<0.05,]$Nrow,ylab='Nrow',names = c('Lm2CoeffP>=0.05','Lm2CoeffP<0.05'), notch = TRUE)

plot(FinalRes$Nrow,-log10(FinalRes$Lm2CoeffP))
points(FinalRes[FinalRes$Gene == "dN.dS_ENSG00000096384",]$Nrow,-log10(FinalRes[FinalRes$Gene == "dN.dS_ENSG00000096384",]$Lm2CoeffP), col = 'red', pch = 16, cex = 2)
cor.test(FinalRes$Nrow,FinalRes$CoeffOfVar, method = 'spearman') # very negative and very significant
plot(FinalRes$Nrow,FinalRes$CoeffOfVar)
points(FinalRes[FinalRes$Gene == "dN.dS_ENSG00000096384",]$Nrow,FinalRes[FinalRes$Gene == "dN.dS_ENSG00000096384",]$CoeffOfVar, col = 'red', pch = 16, cex = 2)


plot(FinalRes$RankCorRho,-log10(FinalRes$RankCorP), xlab = 'RankCorRho', ylab = 'RankCorP'); 
points(FinalRes[FinalRes$Gene == "dN.dS_ENSG00000096384",]$RankCorRho,-log10(FinalRes[FinalRes$Gene == "dN.dS_ENSG00000096384",]$RankCorP), col = 'red', pch = 16, cex = 2);

plot(FinalRes$Lm2CoeffP,FinalRes$CoeffOfVar, xlab = 'Lm2CoeffP', ylab = 'CoeffOfVar'); 

############## 2:  CONTROL SUBSETS OF SPECIES

SpeciesHsp =  as.character(df[!is.na(df$dN.dS_ENSG00000096384),]$Species)
length(SpeciesHsp) # 47

DfShort = df[df$Species %in% SpeciesHsp,]

for (i in 2:(ncol(DfShort)-1))
{
  Temp=df[,c(i,302)]
  Temp=Temp[!is.na(Temp[,1]),]
  Nrow=nrow(Temp)
  if (Nrow == 47)
  {
    Gene=as.character(names(Temp[1]))
    
    # Lm1
    Lm1<-lm(Temp[,1]~Temp$Generation_Length)
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

plot(FinalRes$Lm2Coeff,-log10(FinalRes$Lm2CoeffP), xlab = 'Lm2Coeff', ylab = '-log10(Lm2CoeffP)'); 
points(FinalRes[FinalRes$Gene == "dN.dS_ENSG00000096384",]$Lm2Coeff,-log10(FinalRes[FinalRes$Gene == "dN.dS_ENSG00000096384",]$Lm1CoeffP), col = 'red', pch = 16, cex = 2);
abline(h=-log10(0.05),col = 'red')

plot(FinalRes$RankCorRho,-log10(FinalRes$RankCorP), xlab = 'RankCorRho', ylab = 'RankCorP'); 
points(FinalRes[FinalRes$Gene == "dN.dS_ENSG00000096384",]$RankCorRho,-log10(FinalRes[FinalRes$Gene == "dN.dS_ENSG00000096384",]$RankCorP), col = 'red', pch = 16, cex = 2);

# Now we delete species by species from 47 (start from species which are rarely presented in our dataset) and check at each step the percentile - where hsp is located? 

DfShort$na_count <- apply(DfShort, 1, function(x) sum(!is.na(x)))
summary(DfShort$na_count)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 56.0   222.5   233.0   230.3   251.0   280.0 

DfShort = DfShort[order(DfShort$na_count),] # dataset is sorted from 56 on the top till 280 at the end, we can cut line by line from the top and see results till we reach 10 for example.



dev.off()
