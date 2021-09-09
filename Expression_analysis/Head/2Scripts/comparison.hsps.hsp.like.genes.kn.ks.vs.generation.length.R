rm(list=ls(all=TRUE))


hsps <- read.table('../../Body/2_Derived/kn.ks.all.hsp.txt')
hshs[hshs == 'n/a'] <- NA
hsps[, -1] <- sapply(hsps[, -1], as.vector)
hsps[, -1] <- sapply(hsps[, -1], as.numeric)

df <- read.table('../../Body/2_Derived/kn.ks.genes.like.hsp.txt')
df[df == 'n/a'] <- NA
df[, -1] <- sapply(df[, -1], as.vector)
df[, -1] <- sapply(df[, -1], as.numeric)


#delete columns with too much NA
todrop <- c()
for (col in colnames(df[,-c(1,96)])) {
  if (sum(is.na(df[,col])) > 35) {
    todrop <- c(todrop, col)
  }
}
df <- df[, !(names(df) %in% todrop)]

todrop <- c()
for (col in colnames(hsps[,-c(1,96)])) {
  if (sum(is.na(hsps[,col])) > 35) {
    todrop <- c(todrop, col)
  }
}
hsps <- hsps[, !(names(hsps) %in% todrop)]

df <- df[,!(names(df[,-c(1,258)]) %in% names(hsps[,-c(1,75)]))]

df_mean_kn_ks <- apply(df[,-c(1,254)], 2, mean, na.rm = T)
hsps_mean_kn_ks <- apply(hsps[,-c(1,75)], 2, mean, na.rm = T)

boxplot(df_mean_kn_ks, hsps_mean_kn_ks)
hsps[,'dN.dS_ENSG00000096384']

wilcox.test(df_mean_kn_ks, hsps_mean_kn_ks)
#data:  df_mean_kn_ks and hsps_mean_kn_ks
#W = 5802, p-value = 1.562e-06
#alternative hypothesis: true location shift is not equal to 0



########## model for hsp-like genes
for (i in 2:(ncol(df)-1))
{# i = 71
  Temp=df[,c(i,254)]
  Temp=Temp[!is.na(Temp[,1]),]
  Nrow=nrow(Temp)
  if (Nrow > 10)
    
  {
    Gene=as.character(names(Temp[1]))
    
    # Lm1
    Lm1<-lm(Temp[,1]~log10(Temp$Generation_Length))
    res = data.frame(summary(Lm1)[4])
    Lm1Intercept = res[1,1]; Lm1InterceptP = res[1,4]
    Lm1Coeff = res[2,1]; Lm1CoeffP = res[2,4]
    
    # Lm2
    Lm2<-lm(scale(Temp[,1])~ 0 + log10(scale(Temp$Generation_Length)))
    res = data.frame(summary(Lm2)[4])
    Lm2Coeff = res[1,1]; Lm2CoeffP = res[1,4]
    
    # var
    CoeffOfVar = sd(Temp[,1])/mean(Temp[,1])
    
    # RankCor
    RankCor <-cor.test(Temp[,1],Temp$Generation_Length, method='spearman')
    RankCorP = RankCor$p.value
    RankCorRho = RankCor$estimate
    
    ResLine=data.frame(Gene,Nrow,CoeffOfVar,Lm1Intercept,Lm1InterceptP,Lm1Coeff,Lm1CoeffP,Lm2Coeff,Lm2CoeffP,RankCorP,RankCorRho)
    if (i==2) {FinalRes_hsp_like = ResLine}
    if (i>2) {FinalRes_hsp_like = rbind(FinalRes_hsp_like,ResLine)}
  }
}



######### models for hsps
for (i in 2:(ncol(hsps)-1))
{# i = 71
  Temp=hsps[,c(i,75)]
  Temp=Temp[!is.na(Temp[,1]),]
  Nrow=nrow(Temp)
  if (Nrow > 10)
    
  {
    Gene=as.character(names(Temp[1]))
    
    # Lm1
    Lm1<-lm(Temp[,1]~log10(Temp$Generation_Length))
    res = data.frame(summary(Lm1)[4])
    Lm1Intercept = res[1,1]; Lm1InterceptP = res[1,4]
    Lm1Coeff = res[2,1]; Lm1CoeffP = res[2,4]
    
    # Lm2
    Lm2<-lm(scale(Temp[,1])~ 0 + log10(scale(Temp$Generation_Length)))
    res = data.frame(summary(Lm2)[4])
    Lm2Coeff = res[1,1]; Lm2CoeffP = res[1,4]
    
    # var
    CoeffOfVar = sd(Temp[,1])/mean(Temp[,1])
    
    # RankCor
    RankCor <-cor.test(Temp[,1],Temp$Generation_Length, method='spearman')
    RankCorP = RankCor$p.value
    RankCorRho = RankCor$estimate
    
    ResLine=data.frame(Gene,Nrow,CoeffOfVar,Lm1Intercept,Lm1InterceptP,Lm1Coeff,Lm1CoeffP,Lm2Coeff,Lm2CoeffP,RankCorP,RankCorRho)
    if (i==2) {FinalRes_hsps = ResLine}
    if (i>2) {FinalRes_hsps = rbind(FinalRes_hsps,ResLine)}
  }
}

FinalRes_hsp_like$group <- rep('hsp-like_genes', dim(FinalRes_hsp_like)[1])
FinalRes_hsps$group <- rep('hsp', dim(FinalRes_hsps)[1])
FinalRes <- rbind(FinalRes_hsp_like, FinalRes_hsps)


pdf('comparison.hsps.hsp.like.genes.kn.ks.vs.generation.length.pdf')
library(ggplot2)
p <- ggplot(data = FinalRes, aes(group, Lm1Coeff))
p + geom_violin(aes(fill = group))+
  geom_boxplot(width=0.05)+
  ggtitle('Lm1: lm(Kn/Ks ~ log10(Generation_Length))')

p <- ggplot(data = FinalRes, aes(group, Lm1Intercept))
p + geom_violin(aes(fill = group))+
  geom_boxplot(width=0.05)+
  ggtitle('Lm1: lm(Kn/Ks ~ log10(Generation_Length))')


p <- ggplot(FinalRes, aes(Lm1Coeff, -log10(Lm1CoeffP), shape = group))
p + geom_point(aes(color = group), size = 3, alpha = 0.6)+
  ggtitle('Lm1: lm(Kn/Ks ~ log10(Generation_Length))')

p <- ggplot(FinalRes, aes(Lm1Intercept, -log10(Lm1InterceptP), shape = group))
p + geom_point(aes(color = group), size = 3, alpha = 0.6)+
  ggtitle('Lm1: lm(Kn/Ks ~ log10(Generation_Length))')

p <- ggplot(FinalRes, aes(Lm1Intercept, Lm1Coeff, shape = group))
p + geom_point(aes(color = group), size = 3, alpha = 0.6)+
  ggtitle('Lm1: lm(Kn/Ks ~ log10(Generation_Length))')

##Lm2
p <- ggplot(data = FinalRes, aes(group, Lm2Coeff))
p + geom_violin(aes(fill = group))+
  geom_boxplot(width=0.05)+
  ggtitle('Lm2: lm(Kn/Ks ~ 0 + log10(Generation_Length))')

p <- ggplot(FinalRes, aes(Lm2Coeff, -log10(Lm2CoeffP), shape = group))
p + geom_point(aes(color = group), size = 3, alpha = 0.6)+
  ggtitle('Lm2: lm(Kn/Ks ~ 0 + log10(Generation_Length))')

#Rankcor
p <- ggplot(data = FinalRes, aes(group, RankCorRho))
p + geom_violin(aes(fill = group))+
  geom_boxplot(width=0.05)+
  ggtitle('Spearman\'s Rank Correlation')

p <- ggplot(FinalRes, aes(RankCorRho, -log10(RankCorP), shape = group))
p + geom_point(aes(color = group), size = 3, alpha = 0.6)+
  ggtitle('Spearman\'s Rank Correlation')
dev.off()














