rm(list=ls(all=TRUE))


df <- read.table('../../Body/2_Derived/kn.ks.genes.like.hsp.txt')
df[df == 'n/a'] <- NA
df[, -1] <- sapply(df[, -1], as.vector)
df[, -1] <- sapply(df[, -1], as.numeric)

hsp90 <- na.omit(df[,c('dN.dS_ENSG00000096384', 'Species', 'Generation_Length')])
hsp90_lm <- lm(scale(hsp90$dN.dS_ENSG00000096384) ~  scale(hsp90$Generation_Length), data = hsp90)
summary(hsp90_lm)

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




lm_hsp90_log <- lm((hsp90$dN.dS_ENSG00000096384 + 0.00001) ~  log10(hsp90$Generation_Length + 0.00001), data = hsp90)
summary(lm_hsp90_log)

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                          -10.1334     1.4822  -6.837 1.78e-08 ***
#log(hsp90$Generation_Length + 1e-05)   0.6906     0.1877   3.679 0.000624 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.9313 on 45 degrees of freedom
#Multiple R-squared:  0.2312,	Adjusted R-squared:  0.2141 
#F-statistic: 13.53 on 1 and 45 DF,  p-value: 0.0006238


cor.test(hsp90$dN.dS_ENSG00000096384, hsp90$Generation_Length, method = 'spearman') #S = 8280.8, p-value = 0.0001722, rho = 0.5212286 
cor.test(log(hsp90$dN.dS_ENSG00000096384 + 0.00001), log(hsp90$Generation_Length + 0.00001), method = 'spearman') # S = 8280.8, p-value = 0.0001722, rho = 0.5212286 

pdf('../../Body/4_Figures/hsp.like.genes.linear.model.results.kn.ks.vs.generation.length.mammals.right.way.with.log.pdf')

library('ggplot2')
ggplot(hsp90, aes(log10(Generation_Length + 0.00001), log10(dN.dS_ENSG00000096384 + 0.00001)))+
  geom_point(size = 2)+
  geom_smooth(method = 'lm', se = F, color = 'red')+
  xlab('log(Generation length, days)')+
  ylab('relaxation ~ log(Kn/Ks)')+
  theme_bw()



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
intercept <- c()
p_val_intercept <- c()
R_sq <- c()
R_sq_adj <- c()
number_of_species <- c()
genes <- colnames(df_hsp90_sp)[-c(1,217)]
residual_std_err <- c()



for (gene in genes) {

  a <- na.omit(df_hsp90_sp[,c("Species", gene, 'Generation_Length')])
  
  species <- a[,'Species']
  number_of_species <- c(number_of_species, length(species))
  lg <- lm(log(a[,'Generation_Length'] + 0.00001) ~ log(a[,gene] + 0.00001))
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
results <- data.frame(genes, slopes, p_val_slope,number_of_species, R_sq, R_sq_adj, residual_std_err, intercept, p_val_intercept)

write.table(results, '../../Body/3_Results/hsp.like.genes.linear.model.results.kn.ks.vs.generation.length.mammals.right.way.with.log.txt')



###drawing plots
results <- read.table('../../Body/3_Results/hsp.like.genes.linear.model.results.kn.ks.vs.generation.length.mammals.right.way.with.log.txt')




boxplot(results[,'slopes'], ylab = 'slope')
points(results[results$genes == 'dN.dS_ENSG00000096384', 'slopes'], col = 'red', pch = 19)
title('Slopes of linear regression \"Kn/Ks\" ~ \"generation length\"\n for genes closed to hsp90 (red)\n n = 215')
legend('topright',legend = 'HSP90AB1', col = 'red', pch = 19)


boxplot(results[results$p_val_slope < 0.05,'slopes'], ylab = 'slope')
points(results[(results$genes == 'dN.dS_ENSG00000096384') & (results$p_val_slope < 0.05), 'slopes'], col = 'red', pch = 19)
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
abline(h = 0.01, col = 'red', lwd = 2)

boxplot(results[results$p_val_slope < 0.01,'slopes'], ylab = 'relaxarion ~ slope')
points(results[(results$genes == 'dN.dS_ENSG00000096384') & (results$p_val_slope < 0.01), 'slopes'], col = 'red', pch = 19)
title('Slopes of linear regression \"Kn/Ks\" ~ \"generation length\"\n for genes closed to hsp90 (red)\n n = 73')
legend('topright',legend = 'HSP90AB1', col = 'red', pch = 19)



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
cor.test(results[,'intercept'], results[,'slopes']) #0.9249095






#Spearman's rank corr
spear_cor <- cor.test(hsp90$dN.dS_ENSG00000096384, hsp90$Generation_Length, method = 'spearman')

statistic <- c()
p_val <- c()
rho <- c()
number_of_species <- c()

for (gene in genes) {
  
  a <- na.omit(df_hsp90_sp[,c("Species", gene, 'Generation_Length')])
  
  species <- a[,'Species']
  number_of_species <- c(number_of_species, length(species))
  rank_cor <- cor.test(a[,gene], a[,'Generation_Length'], method = 'spearman' )
  statistic <- c(statistic, rank_cor$statistic)
  p_val <- c(p_val, rank_cor$p.value)
  rho <- c(rho, rank_cor$estimate)


}

result_spearman <- data.frame(genes, statistic, p_val, rho, number_of_species)


boxplot(result_spearman$rho, ylab = 'Spearman\'s rho')
points(result_spearman[result_spearman$genes == 'dN.dS_ENSG00000096384', 'rho'], pch = 19, col = 'red')
legend('topright', legend = 'HSP90AB1', col = 'red', pch = 19)


plot(result_spearman$rho, result_spearman$p_val)
points(result_spearman[result_spearman$genes == 'dN.dS_ENSG00000096384', 'rho'], result_spearman[result_spearman$genes == 'dN.dS_ENSG00000096384', 'p_val'], pch = 19, col = 'red')
legend('topright', legend = 'HSP90AB1', col = 'red', pch = 19)


length(results[results$p_val_slope < 0.01,'slopes'])

hsp <- results[(results$genes == 'dN.dS_ENSG00000096384') & (results$p_val_slope < 0.01), 'slopes']
length(results[(results$p_val_slope < 0.01) & (results$slopes <= hsp), 'slopes'])/length(results[results$p_val_slope < 0.01,'slopes']) #0.2191781



###z-statistic


qqnorm(results$slopes)
qqline(results$slopes)
shapiro.test(results$slopes) #p-value = 1.413e-08

qqnorm(results[results$p_val_slope < 0.01, 'slopes'])
qqline(results[results$p_val_slope < 0.01, 'slopes'])
shapiro.test(results[results$p_val_slope < 0.01, 'slopes']) #p-value = 1.032e-06
shapiro.test(log(results[results$p_val_slope < 0.01, 'slopes'])) #p-value = 0.0297

library(TeachingDemos)

z.test(results[results$p_val_slope < 0.01, 'slopes'], mu = hsp, alternative = 'two.sided', stdev = sd(results[results$p_val_slope < 0.01, 'slopes']))
#One Sample z-test

#data:  results[results$p_val_slope < 0.01, "slopes"]
#z = 5.0595, n = 73.000000, Std. Dev. = 0.469030, Std. Dev. of the sample mean = 0.054896, p-value = 4.203e-07
#alternative hypothesis: true mean is not equal to 0.3804685
#95 percent confidence interval:
#  0.5506204 0.7658071

z.test(results[results$p_val_slope < 0.01, 'slopes'], mu = hsp, alternative = 'greater', stdev = sd(results[results$p_val_slope < 0.01, 'slopes']))
#One Sample z-test

#data:  results[results$p_val_slope < 0.01, "slopes"]
#z = 5.0595, n = 73.000000, Std. Dev. = 0.469030, Std. Dev. of the sample mean = 0.054896, p-value = 2.102e-07
#alternative hypothesis: true mean is greater than 0.3804685
#95 percent confidence interval:
#  0.5679186       Inf
#sample estimates:
#  mean of results[results$p_val_slope < 0.01, "slopes"] 
#0.6582137

z.test(na.omit(log(results[results$p_val_slope < 0.01, 'slopes'])), mu = log(hsp), alternative = 'two.sided', stdev = sd(na.omit(log(results[results$p_val_slope < 0.01, 'slopes']))))
#One Sample z-test

#data:  na.omit(log(results[results$p_val_slope < 0.01, "slopes"]))
#z = 8.2762, n = 71.000000, Std. Dev. = 0.492190, Std. Dev. of the sample mean = 0.058412, p-value < 2.2e-16
#alternative hypothesis: true mean is not equal to -0.9663519
#95 percent confidence interval:
#  -0.5974082 -0.3684369

z.test(na.omit(log(results[results$p_val_slope < 0.01, 'slopes'])), mu = log(hsp), alternative = 'greater', stdev = sd(na.omit(log(results[results$p_val_slope < 0.01, 'slopes']))))

#One Sample z-test

#data:  na.omit(log(results[results$p_val_slope < 0.01, "slopes"]))
#z = 8.2762, n = 71.000000, Std. Dev. = 0.492190, Std. Dev. of the sample mean = 0.058412, p-value < 2.2e-16
#alternative hypothesis: true mean is greater than -0.9663519
#95 percent confidence interval:
#  -0.579002       Inf



wilcox.test(results[results$p_val_slope < 0.01, 'slopes'], mu = hsp, alternative = "two.sided")
#Wilcoxon signed rank test with continuity correction

#data:  results[results$p_val_slope < 0.01, "slopes"]
#V = 2318, p-value = 1.788e-08
#alternative hypothesis: true location is not equal to 0.3804685


wilcox.test(results[results$p_val_slope < 0.01, 'slopes'], mu = hsp, alternative = "greater")

#data:  results[results$p_val_slope < 0.01, "slopes"]
#V = 2318, p-value = 8.94e-09
#alternative hypothesis: true location is greater than 0.3804685

pdf('../../Body/4_Figures/hsp_like_genes_log10_linear_model_slopes_with_log_nominally_sign.pdf')
boxplot(results[results$p_val_slope < 0.01,'slopes'], ylab = 'relaxarion ~ slope', col = 'white', cex.lab=2.2, cex.axis = 1.7)
points(results[(results$genes == 'dN.dS_ENSG00000096384') & (results$p_val_slope < 0.01), 'slopes'], col = 'red', pch = 19, cex = 1.3)
#title('Slopes of linear regression \"Kn/Ks\" ~ \"generation length\"\n for genes closed to hsp90 (red)\n n = 73')
legend('topright',legend = 'HSP90', col = 'red', pch = 19, cex = 2)
dev.off()
