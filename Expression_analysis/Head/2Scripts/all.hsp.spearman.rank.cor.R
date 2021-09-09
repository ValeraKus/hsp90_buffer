rm(list=ls(all=TRUE))


df <- read.table('../../Body/2_Derived/kn.ks.all.hsp.txt')


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

rho <- c()
S <- c()
p_val <- c()
number_of_species <- c()
genes <- colnames(df)[-c(1,75)]



for (gene in genes){
  a <- na.omit(df[,c("Species", gene, 'Generation_Length')])
  
  species <- a[,'Species']
  number_of_species <- c(number_of_species, length(species))
  
  rcor <- cor.test(~ a[,gene] + a[,'Generation_Length'], method = 'spearman', continuity = F, conf.level = 0.95)
  rho <- c(rho, rcor$estimate)
  S <- c(S, rcor$statistic)
  p_val <- c(p_val, rcor$p.value)
}


results <- data.frame(genes, S, rho, p_val, number_of_species)
write.table(results, '../../Body/3_Results/all.hsp.spearman.rank.cor.txt')


rm(list=ls(all=TRUE))
results <- read.table('../../Body/3_Results/all.hsp.spearman.rank.cor.txt')

#hsp_groups <- read.table('../../Body/2_Derived/human.hsp.ensID.group.txt', sep = '\t', header = T)
#results$genes <- sub('dN.dS_', '', results$genes)
#results$group <- hsp_groups$Group.name[match(results$genes, hsp_groups$Ensembl.gene.ID)]

pdf('../../Body/4_Figures/all.hsp.spearman.rank.cor.pdf')
plot(~results$rho + results$p_val, col = results$group, pch = 19)
legend('topright',legend = unique(results$group),col=1:length(results$group),pch=19)
title('rank Spearman correlation for all hsps')


par(yat)
boxplot(rho~group, data=results, ylab ='Spearman rho')
stripchart(rho~group, data=results,vertical=TRUE, add=TRUE, method="stack", pch =1)
title('rank Spearman correlation for all hsps')
dev.off()
