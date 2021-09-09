rm(list=ls(all=TRUE))
df <- read.table('../../Body/2_Derived/kn.ks.genes.like.hsp.txt')

df[df == 'n/a'] <- NA
df[, -1] <- sapply(df[, -1], as.vector)
df[, -1] <- sapply(df[, -1], as.numeric)

hsp90 <- na.omit(df[,c('dN.dS_ENSG00000096384', 'Species', 'Generation_Length')])
df_hsp90_sp <- df[df$Species %in% hsp90$Species,] #to make model for other genes using more or less same speciaes as hsp90




cols <- colnames(df_hsp90_sp[,-c(1, 302)])
for (c in cols) {
  if (length(df_hsp90_sp[is.na(df_hsp90_sp[,c]),c]) >= 15){
    df_hsp90_sp[,c] <- NULL
  }
  
  
}
rows <- rownames(df_hsp90_sp)
for (r in rows) {
  if (length(df_hsp90_sp[r, is.na(df_hsp90_sp[r,])]) >= 15){
    df_hsp90_sp <- df_hsp90_sp[-as.numeric(r),]
  }
}

rho <- c()
S <- c()
p_val <- c()
number_of_species <- c()
genes <- colnames(df_hsp90_sp)[-c(1,217)]



for (gene in genes){
  a <- na.omit(df_hsp90_sp[,c("Species", gene, 'Generation_Length')])
  
  species <- a[,'Species']
  number_of_species <- c(number_of_species, length(species))
  
  rcor <- cor.test(~ a[,gene] + a[,'Generation_Length'], method = 'spearman', continuity = F, conf.level = 0.95)
  rho <- c(rho, rcor$estimate)
  S <- c(S, rcor$statistic)
  p_val <- c(p_val, rcor$p.value)
}


results <- data.frame(genes, S, rho, p_val, number_of_species)
write.table(results, '../../Body/3_Results/hsp.like.genes.sperman.rank.cor.txt')


pdf('../../Body/4_Figures/hsp.like.genes.spearman.rank.cor.less.15.sp.pdf')
plot(~results$rho + results$p_val)
points(results[results$genes=='dN.dS_ENSG00000096384', 'rho'], results[results$genes == 'dN.dS_ENSG00000096384', 'p_val'], col = 'red', pch = 19)
legend('topright', legend = 'HSP90AB1', col = 'red', pch = 19)
title('rank Spearman correlation for genes like HSP90AB1, n = 215')


boxplot(results$rho, ylab = 'Spearman\'s rho')
points(results[results$genes=='dN.dS_ENSG00000096384', 'rho'], col = 'red', pch = 19)
legend('topright', legend = 'HSP90AB1', col = 'red', pch = 19)
title('rank Spearman correlation for genes like HSP90AB1, n = 215')
dev.off()
