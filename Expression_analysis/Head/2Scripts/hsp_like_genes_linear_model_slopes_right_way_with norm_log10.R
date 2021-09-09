rm(list=ls(all=TRUE))


df <- read.table('../../Body/2_Derived/kn.ks.genes.like.hsp.txt')
df[df == 'n/a'] <- NA
df[, -1] <- sapply(df[, -1], as.vector)
df[, -1] <- sapply(df[, -1], as.numeric)

hsp90 <- na.omit(df[,c('dN.dS_ENSG00000096384', 'Species', 'Generation_Length')])
hsp90_lm <- lm(scale(hsp90$dN.dS_ENSG00000096384) ~  scale(hsp90$Generation_Length), data = hsp90)
summary(hsp90_lm)
#Call:
#  lm(formula = scale(hsp90$dN.dS_ENSG00000096384) ~ scale(hsp90$Generation_Length), data = hsp90)


R version 3.4.4 (2018-03-15) -- "Someone to Lean On"
Copyright (C) 2018 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

Natural language support but running in an English locale

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

[Workspace loaded from ~/Desktop/HSP/Head/2Scripts/.RData]

> rm(list=ls(all=TRUE))
> df <- read.table('../../Body/2_Derived/kn.ks.genes.like.hsp.txt')
> df[df == 'n/a'] <- NA
> df[, -1] <- sapply(df[, -1], as.vector)
> df[, -1] <- sapply(df[, -1], as.numeric)
> hsp90 <- na.omit(df[,c('dN.dS_ENSG00000096384', 'Species', 'Generation_Length')])
> hsp90_lm <- lm(scale(hsp90$dN.dS_ENSG00000096384) ~  scale(hsp90$Generation_Length), data = hsp90)
> summary(hsp90_lm)

#Call:
#  lm(formula = scale(hsp90$dN.dS_ENSG00000096384) ~ scale(hsp90$Generation_Length), 
#     data = hsp90)

#Residuals:
#  Min      1Q  Median      3Q     Max 
#-1.2453 -0.4566 -0.3023 -0.1628  3.8155 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)  
#(Intercept)                    6.478e-17  1.413e-01   0.000    1.000  
#scale(hsp90$Generation_Length) 2.864e-01  1.428e-01   2.005    0.051 .
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.9687 on 45 degrees of freedom
#Multiple R-squared:  0.08201,	Adjusted R-squared:  0.06161 
#F-statistic:  4.02 on 1 and 45 DF,  p-value: 0.051

hsp90_lm1 <- lm(scale(hsp90$dN.dS_ENSG00000096384) ~ 0 + scale(hsp90$Generation_Length), data = hsp90)
summary(hsp90_lm1)
#Call:
#  lm(formula = scale(hsp90$dN.dS_ENSG00000096384) ~ 0 + scale(hsp90$Generation_Length), 
#     data = hsp90)

#Residuals:
#  Min      1Q  Median      3Q     Max 
#-1.2453 -0.4566 -0.3023 -0.1628  3.8155 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)  
#scale(hsp90$Generation_Length)   0.2864     0.1413   2.027   0.0485 *
  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.9581 on 46 degrees of freedom
#Multiple R-squared:  0.08201,	Adjusted R-squared:  0.06205 
#F-statistic: 4.109 on 1 and 46 DF,  p-value: 0.04846

plot(x = scale(hsp90$Generation_Length), y = scale(hsp90$dN.dS_ENSG00000096384), pch = 19, xlab = 'Generation Length, days', ylab = 'Kn/Ks')  
abline(hsp90_lm$coefficients[1], hsp90_lm$coefficients[2])  
abline(0, hsp90_lm1$coefficients[1], col = 'red') 
legend('topright', col = c('black', 'red'), pch = c('-', '-'), legend = c('Kn/Ks ~ Generation Length, R^2 = 0.082', 'Kn/Ks ~ 0 + Generation Length, R^2 = 0.082'))
title('Linear regression of Kn/Ks vs Generation Length for HSP90AB1')

library('ggplot2')
ggplot(hsp90, aes(Generation_Length, dN.dS_ENSG00000096384))+
  geom_point(size = 2)+
  geom_smooth(method = 'lm')



df_hsp90_sp <- df[df$Species %in% hsp90$Species,] #to make linear model for other genes using more or less same species as hsp90




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

#perform linear model
slopes <- c()
p_val_slope <- c()
#intercept <- c()
#p_val_intercept <- c()
R_sq <- c()
R_sq_adj <- c()
number_of_species <- c()
genes <- colnames(df_hsp90_sp)[-c(1,217)]
residual_std_err <- c()



for (gene in genes) {

  a <- na.omit(df_hsp90_sp[,c("Species", gene, 'Generation_Length')])
  
  species <- a[,'Species']
  number_of_species <- c(number_of_species, length(species))
  lg <- lm(scale(a[,'Generation_Length']) ~ 0 + scale(a[,gene]))
  sum <- summary(lg)
  b <- sum$coefficients
  
  slopes <- c(slopes, b[1])
  p_val_slope <- c(p_val_slope, b[4])
  #intercept <- c(intercept, b[1,1])
  #p_val_intercept <- c(p_val_intercept, b[1,4])
  R_sq <- c(R_sq, sum$r.squared)
  R_sq_adj <- c(R_sq_adj, sum$adj.r.squared)
  residual_std_err <- c(residual_std_err, sum$sigma)
  
}


#table with linear regression parameters
results <- data.frame(genes, slopes, p_val_slope,number_of_species, R_sq, R_sq_adj, residual_std_err)

write.table(results, '../../Body/3_Results/hsp.like.genes.linear.model.results.kn.ks.vs.generation.length.mammals.right.way.with.norm.log10.txt')



###drawing plots
results <- read.table('../../Body/3_Results/hsp.like.genes.linear.model.results.kn.ks.vs.generation.length.mammals.right.way.with.norm.log10.txt')


pdf('../../Body/4_Figures/hsp.like.genes.linear.model.slopes.right.way.with.scale.less.15.species.from.hsp.zerro.intercept.pdf')

boxplot(results[,'slopes'], ylab = 'slope')
points(results[results$genes == 'dN.dS_ENSG00000096384', 'slopes'], col = 'red', pch = 19)
title('Slopes of linear regression \"Kn/Ks\" ~ \"generation length\"\n for genes closed to hsp90 (red)\n n = 215')
legend('topright',legend = 'HSP90AB1', col = 'red', pch = 19)


boxplot(results[,'intercept'], ylab = 'intercept')
points(results[results$genes == 'dN.dS_ENSG00000096384', 'intercept'], col = 'red', pch = 19)
title('Inercepts of linear regression \"Kn/Ks\" ~ \"generation length\"\n for genes closed to hsp90 (red)\n n = 215')
legend('topright', legend = 'HSP90AB1', col = 'red', pch = 19)

plot(results[,'slopes'], results[,'p_val_slope'],
     xlab = 'slope', ylab = 'p_val_slope')
points(results[results$genes == 'dN.dS_ENSG00000096384','slopes'], results[results$genes == 'dN.dS_ENSG00000096384','p_val_slope'], col = 'red', pch = 19)
title('Slope vs intercept genes like hsp90, n = 215')     
legend('topright', legend = 'HSP90AB1', col = 'red', pch = 19)


plot(results[,'intercept'], results[,'slopes'],
     xlab = 'intecept', ylab = 'slope')
points(results[results$genes == 'dN.dS_ENSG00000096384','intercept'], results[results$genes == 'dN.dS_ENSG00000096384','slopes'], col = 'red', pch = 19)
title('Slope vs intercept genes like hsp90, n = 215')     
legend('topright', legend = 'HSP90AB1', col = 'red', pch = 19)

##########


plot(results$p_val_slope, results$p_val_intercept, xlab = 'p_val_slope', ylab= 'p_val_intercept')
points(results[results$genes == 'dN.dS_ENSG00000096384', 'p_val_slope'], results[results$genes == 'dN.dS_ENSG00000096384', 'p_val_intercept'], col = 'red', pch = 19)
legend('topright', legend = 'HSP90AB1', col = 'red', pch = 19)

dev.off()



slint <- lm(results$slopes~results$intercept)
cor.test(results[,'intercept'], results[,'slopes']) #-0.8930672 
