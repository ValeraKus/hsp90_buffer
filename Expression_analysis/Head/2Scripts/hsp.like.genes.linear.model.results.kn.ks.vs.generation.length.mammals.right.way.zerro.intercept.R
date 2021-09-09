rm(list=ls(all=TRUE))


df <- read.table('../../Body/2_Derived/kn.ks.genes.like.hsp.txt')
df[df == 'n/a'] <- NA
df[, -1] <- sapply(df[, -1], as.vector)
df[, -1] <- sapply(df[, -1], as.numeric)

hsp90 <- na.omit(df[,c('dN.dS_ENSG00000096384', 'Species', 'Generation_Length')])
df_hsp90_sp <- df[df$Species %in% hsp90$Species,] #to make linear model for other genes using more or less same speciaes as hsp90




cols <- colnames(df_hsp90_sp[,-c(1, 302)])
for (c in cols) {
  if (length(df_hsp90_sp[is.na(df_hsp90_sp[,c]),c]) >= 10){
    df_hsp90_sp[,c] <- NULL
  }
  
  
}
rows <- rownames(df_hsp90_sp)
for (r in rows) {
  if (length(df_hsp90_sp[r, is.na(df_hsp90_sp[r,])]) >= 10){
    df_hsp90_sp <- df_hsp90_sp[-as.numeric(r),]
  }
}

#perform linear model
slopes <- c()
p_val_slope <- c()

R_sq <- c()
R_sq_adj <- c()
number_of_species <- c()
genes <- colnames(df_hsp90_sp)[-c(1,179)]
residual_std_err <- c()
df <- c()


for (gene in genes) {
  
  a <- na.omit(df_hsp90_sp[,c("Species", gene, 'Generation_Length')])
  
  species <- a[,'Species']
  number_of_species <- c(number_of_species, length(species))
  
  lg <- lm(a[,gene]~0+a[,'Generation_Length'], data = a)
  sum <- summary(lg)
  b <- sum$coefficients
  
  slopes <- c(slopes, b[1])
  p_val_slope <- c(p_val_slope, b[4])

  R_sq <- c(R_sq, sum$r.squared)
  R_sq_adj <- c(R_sq_adj, sum$adj.r.squared)
  residual_std_err <- c(residual_std_err, sum$sigma)
  
}


#table with linear regression parameters
results <- data.frame(genes, slopes, p_val_slope, number_of_species, R_sq, R_sq_adj, residual_std_err)

write.table(results, '../../Body/3_Results/hsp.like.genes.linear.model.results.kn.ks.vs.generation.length.mammals.right.way.zerro.intercept.txt')



###drawing plots
results <- read.table('../../Body/3_Results/hsp.like.genes.linear.model.results.kn.ks.vs.generation.length.mammals.right.way.zerro.intercept.txt')


pdf('../../Body/4_Figures/hsp.like.genes.linear.model.slopes.right.way.zerro.intercept.pdf')

boxplot(results[,'slopes'], ylab = 'slope')
points(results[results$genes == 'dN.dS_ENSG00000096384', 'slopes'], col = 'red', pch = 19)
title('Slopes of linear regression \"Kn/Ks\" ~ \"generation length\"\n for genes closed to hsp90 (red)\n n = 177, intercept = 0')
legend(x = 0, y = 8e-05,legend = 'hsp90', col = 'red', pch = 19)

boxplot(results[results$p_val_slope < 0.05,'slopes'], ylab = 'slope')
points(results[results$genes == 'dN.dS_ENSG00000096384', 'slopes'], col = 'red', pch = 19)
title('Slopes of linear regression \"Kn/Ks\" ~ \"generation length\"\n for genes closed to hsp90 (red)\n n = 177, , p_value_slope < 0.05, intercept = 0')
legend(x = 50, y = 0.005,legend = 'hsp90', col = 'red', pch = 19)


##########





dev.off()
