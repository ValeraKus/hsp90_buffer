rm(list=ls(all=TRUE))


df <- read.table('../../Body/2_Derived/kn.ks.genes.like.hsp.txt')


#gene_ids <- gsub('dN.dS_', '', colnames(df)[-c(1,302)])

df[df == 'n/a'] <- NA
df[, -1] <- sapply(df[, -1], as.vector)
df[, -1] <- sapply(df[, -1], as.numeric)

slopes <- c()
p_val_slope <- c()
intercept <- c()
p_val_intercept <- c()

number_of_species <- c()
genes <- colnames(df)[-c(1,302)]


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

}





results <- data.frame(genes, slopes, intercept, p_val_slope, p_val_intercept, number_of_species)

write.table(results, '../../Body/3_Results/linear.regression.kn.ks.vs.generation.length.hsp.like.genes.mammals.txt'            )


results <- read.table('../../Body/3_Results/linear.regression.kn.ks.vs.generation.length.hsp.like.genes.mammals.txt')

pdf('../../Body/4_Figures/hsp.like.genes.linear.model.slopes.pdf')
boxplot(results[results$p_val_slope < 0.05,'slopes'], ylab = 'slope')
points(results[results$genes == 'dN.dS_ENSG00000096384', 'slopes'], col = 'red', pch = 19)
title('Slopes of linear regression \"Kn/Ks\" ~ \"generation length\"\n for genes closed to hsp90 (red)\n n=155, p_value_slope < 0.05')
legend(x = 50, y = 0.005,legend = 'hsp90', col = 'red', pch = 19)


boxplot(results[results$p_val_slope < 0.05,'intercept'], ylab = 'intercept')
points(results[results$genes == 'dN.dS_ENSG00000096384', 'intercept'], col = 'red', pch = 19)
title('Inercepts of linear regression \"Kn/Ks\" ~ \"generation length\"\n for genes closed to hsp90 (red)\n n=155, p_value_slope < 0.05')
legend(x = 10, y = 35, legend = 'hsp90', col = 'red', pch = 19)



plot(results[results$p_val_slope < 0.05,'intercept'], results[results$p_val_slope < 0.05,'slopes'],
     xlab = 'intecept', ylab = 'slope')
points(results[results$genes == 'dN.dS_ENSG00000096384','intercept'], results[results$genes == 'dN.dS_ENSG00000096384','slopes'], col = 'red', pch = 19)
title('Slope vs intercept genes like hsp90, n = 155, p_value_slope < 0.05')     
legend(x = 32, y = 0.0056, legend = 'hsp90', col = 'red', pch = 19)


#######


boxplot(results[results$p_val_slope < 0.01,'slopes'], ylab = 'slope')
points(results[results$genes == 'dN.dS_ENSG00000096384', 'slopes'], col = 'red', pch = 19)
title('Slopes of linear regression \"Kn/Ks\" ~ \"generation length\"\n for genes closed to hsp90 (red)\n n=109, p_value_slope < 0.01')
legend(x = 50, y = 0.005,legend = 'hsp90', col = 'red', pch = 19)


boxplot(results[results$p_val_slope < 0.01,'intercept'], ylab = 'intercept')
points(results[results$genes == 'dN.dS_ENSG00000096384', 'intercept'], col = 'red', pch = 19)
title('Inercepts of linear regression \"Kn/Ks\" ~ \"generation length\"\n for genes closed to hsp90 (red)\n n=109, p_value_slope < 0.01')
legend(x = 10, y = 35, legend = 'hsp90', col = 'red', pch = 19)



plot(results[results$p_val_slope < 0.01,'intercept'], results[results$p_val_slope < 0.01,'slopes'],
     xlab = 'intecept', ylab = 'slope', xlim = c(0, 30))
points(results[results$genes == 'dN.dS_ENSG00000096384','intercept'], results[results$genes == 'dN.dS_ENSG00000096384','slopes'], col = 'red', pch = 19)
title('Slope vs intercept genes like hsp90, n = 109, p_value_slope < 0.01')     
legend(x = 21, y = 0.006, legend = 'hsp90', col = 'red', pch = 19)
dev.off()


