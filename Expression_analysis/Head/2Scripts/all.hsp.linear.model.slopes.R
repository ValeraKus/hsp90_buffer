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

#saving linear regression parameters
slopes <- c()
p_val_slope <- c()
intercept <- c()
p_val_intercept <- c()
R_sq <- c()
R_sq_adj <- c()
residual_std_err <- c()
number_of_species <- c()
genes <- colnames(df)[-c(1,75)]


for (gene in genes) {
  #a <- na.omit(df[,c(gene, 'Generation_Length')])
  if ((length(df[,gene]) - sum(is.na(df[,gene]))) < 10){
    print(gene)
    genes <- genes[genes != gene]
    next
  }
  a <- na.omit(df[,c("Species", gene, 'Generation_Length')])
  
  species <- a[,'Species']
  number_of_species <- c(number_of_species, length(species))
  
  lg <- lm(a[,gene]~a[,'Generation_Length'], data = a)
  sum <- summary(lg)
  b <- sum$coefficients
  
  slopes <- c(slopes, b[2,1])
  p_val_slope <- c(p_val_slope, b[2,4])
  intercept <- c(intercept, b[1,1])
  p_val_intercept <- c(p_val_intercept, b[1,4])
  R_sq <- c(R_sq, sum$r.squared)
  R_sq_adj <- c(R_sq_adj, sum$adj.r.squared)
  residual_std_err <- c(residual_std_err, sum$sigma)
  
}




#table with linear regression parameters
results <- data.frame(genes, slopes, intercept, p_val_slope, p_val_intercept, number_of_species, R_sq, R_sq_adj, residual_std_err)

write.table(results, '../../Body/3_Results/linear.regression.kn.ks.vs.generation.length.all.hsp.mammals.txt'            )




#### drawing plots
results <- read.table('../../Body/3_Results/linear.regression.kn.ks.vs.generation.length.all.hsp.mammals.txt')

#should we filter by p_value or not?

hsps  <- read.table('../../Body/2_Derived/human.hsp.ensID.group.txt', header = T, sep = '\t')
hsps$Ensembl.gene.ID <- paste('dN.dS_', hsps$Ensembl.gene.ID, sep = '')
results <- merge(x = results, y = hsps, by.x = 'genes', by.y = 'Ensembl.gene.ID')

library(ggplot2)

theme_set(theme_bw())
gg <- ggplot(results, aes(intercept, slopes))+
  geom_point(aes(col = Group.name))+
  xlab('intercept')+
  ylab('slope')+
  labs(title ='Slope vs intercept all hsp genes', subtitle = 'n=73')
print(gg)

pdf('../../Body/4_Figures/all.hsp.linear.model.slopes.pdf')
boxplot(results[,'slopes'], ylab = 'slope')
points(results[results$genes == 'dN.dS_ENSG00000096384', 'slopes'], col = 'red', pch = 19)
title('Slopes of linear regression \"Kn/Ks\" ~ \"generation length\"\n for all hsp genes\n n=73')
legend(x=1.3, y = 0.00015,legend = 'hsp90', col = 'red', pch = 19)


boxplot(results[,'intercept'], ylab = 'intercept')
points(results[results$genes == 'dN.dS_ENSG00000096384', 'intercept'], col = 'red', pch = 19)
title('Inercepts of linear regression \"Kn/Ks\" ~ \"generation length\"\n for all hsp genes\n n=73')
legend(x = 1.3, y = 0.37, legend = 'hsp90', col = 'red', pch = 19)





plot(results[,'intercept'], results[,'slopes'],
     xlab = 'intecept', ylab = 'slope')
points(results[results$genes == 'dN.dS_ENSG00000096384','intercept'], results[results$genes == 'dN.dS_ENSG00000096384','slopes'], col = 'red', pch = 19)
title('Slope vs intercept all hsp genes\n n=73')     
legend(x = 0.35, y = 0.00015, legend = 'hsp90', col = 'red', pch = 19)


print(gg)
#######


boxplot(results[results$p_val_slope < 0.05,'slopes'], ylab = 'slope')
points(results[results$genes == 'dN.dS_ENSG00000096384', 'slopes'], col = 'red', pch = 19)
title('Slopes of linear regression \"Kn/Ks\" ~ \"generation length\"\n for all hsp genes\n n=59, p_value_slope < 0.05')
legend(x = 1.3, y = 0.00015,legend = 'hsp90', col = 'red', pch = 19)


boxplot(results[results$p_val_slope < 0.05,'intercept'], ylab = 'intercept')
points(results[results$genes == 'dN.dS_ENSG00000096384', 'intercept'], col = 'red', pch = 19)
title('Inercepts of linear regression \"Kn/Ks\" ~ \"generation length\"\n for all hsp genes\n n=59, p_value_slope < 0.05')
legend(x = 1.3, y = 0.37, legend = 'hsp90', col = 'red', pch = 19)



plot(results[results$p_val_slope < 0.05,'intercept'], results[results$p_val_slope < 0.05,'slopes'],
     xlab = 'intecept', ylab = 'slope')
points(results[results$genes == 'dN.dS_ENSG00000096384','intercept'], results[results$genes == 'dN.dS_ENSG00000096384','slopes'], col = 'red', pch = 19)
title('Slope vs intercept for all hsp genes\n n=59, p_value_slope < 0.05')     
legend(x = 0.35, y = 0.00015, legend = 'hsp90', col = 'red', pch = 19)
dev.off()

