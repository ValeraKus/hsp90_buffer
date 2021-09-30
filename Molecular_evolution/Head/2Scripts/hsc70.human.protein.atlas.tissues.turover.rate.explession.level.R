rm(list=ls(all=TRUE))

hsp_hpa_rna <- read.csv('../../Body/2_Derived/human_protein_atlas_hsp_expression_level.csv', index = 1)
hsc70_hpa_rna <- read.csv('../../Body/2_Derived/human_protein_atlas_hsc70_expression_level.csv')
colnames(hsc70_hpa_rna) <- c('N', 'tissue', 'expression_level', 'turnover_rate_days')
hsc70_hpa_prot <- read.csv('../../Body/2_Derived/human_protein_atlas_hsc70_protein_level.csv')



pdf('../../Body/4_Figures/hsc70.human.protein.atlas.tissues.turover.rate.explession.level.pdf')

library(ggplot2)

ggplot(data = hsc70_hpa_rna, aes(log10(turnover_rate_days), expression_level, label=tissue))+
  geom_point()+
  geom_text(check_overlap = TRUE, vjust = -0.2)+
  geom_smooth(method = lm)+
  ggtitle('Human Protein Atlas RNA for HSP90AB1')+
  theme_classic()

ggplot(data = na.omit(hsc70_hpa_rna), aes(x = (log10(turnover_rate_days) > 3) , y = expression_level))+
  geom_boxplot()+
  theme_classic()+
  ggtitle('Human Protein Atlas RNA for HSP90AB1')


ggplot(data = hsc70_hpa_prot, aes(x = Level, y = log10(Turnover_days)))+
  geom_boxplot()+
  theme_classic()+
  ggtitle('Human Protein Atlas protein for HSP90AB1')

dev.off()

cor.test(hsc70_hpa_rna$expression_level, hsc70_hpa_rna$turnover_rate_days, method = 'spearman', paired = T)
#Spearman's rank correlation rho

#data:  hsc70_hpa_rna$expression_level and hsc70_hpa_rna$turnover_rate_days
#S = 6831.1, p-value = 0.00863
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#0.4049555 

rna_lm <- lm(hsc70_hpa_rna$expression_level ~ log10(hsc70_hpa_rna$turnover_rate_days))
summary(rna_lm)
#Call:
#lm(formula = hsc70_hpa_rna$expression_level ~ log10(hsc70_hpa_rna$turnover_rate_days))

#Residuals:
#  Min      1Q  Median      3Q     Max 
#-58.195 -22.751  -7.986   0.873 120.366 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                               75.543     11.990    6.30 1.98e-07 ***
#  log10(hsc70_hpa_rna$turnover_rate_days)    7.669      3.995    1.92   0.0622 .  
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 38.92 on 39 degrees of freedom
#(20 observations deleted due to missingness)
#Multiple R-squared:  0.08634,	Adjusted R-squared:  0.06291 
#F-statistic: 3.685 on 1 and 39 DF,  p-value: 0.06223


ggplot(data = na.omit(hsc70_hpa_rna), aes(x = (log10(turnover_rate_days) > 3) , y = expression_level))+
  geom_boxplot()+
  theme_classic()

wilcox.test(hsc70_hpa_rna$expression_level ~ (log10(hsc70_hpa_rna$turnover_rate_days) > 3))
#WWilcoxon rank sum exact test

#data:  hsc70_hpa_rna$expression_level by log10(hsc70_hpa_rna$turnover_rate_days) > 3
#W = 113, p-value = 0.01542
#alternative hypothesis: true location shift is not equal to 0




