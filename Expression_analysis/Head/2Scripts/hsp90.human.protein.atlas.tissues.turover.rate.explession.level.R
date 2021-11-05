rm(list=ls(all=TRUE))

hsp90_hpa_rna <- read.csv('../../Body/2_Derived/human_protein_atlas_hsp_expression_level.csv')
hsp_hpa_prot <- read.csv('../../Body/2_Derived/human_protein_atlas_hsp_protein_level.csv')

pdf('../../Body/4_Figures/hsp90.human.protein.atlas.tissues.turover.rate.explession.level.pdf')

library(ggplot2)

p1 <- ggplot(data = hsp90_hpa_rna, aes(log10(turnover_rate_days), expression_level, label=tissue))+
  geom_point()+
  geom_text(check_overlap = TRUE, vjust = -0.2,size=6)+
  geom_smooth(method = lm)+
  ylab('Expression level of HSP90, NX')+
  xlab('Turnover rate, days')+
  xlim(0,5.2)+
  #ggtitle('Human Protein Atlas RNA for HSP90AB1')+
    theme_classic()+
  theme(axis.title = element_text(size = 27), axis.text = element_text(size = 25))

print(p1)
ggsave('../../Body/4_Figures/hsp90.human.protein.atlas.tissues.turover.rate.explession.level.linear.reg.pdf', p1)


ggplot(data = na.omit(hsp90_hpa_rna), aes(x = (log10(turnover_rate_days) > 3) , y = expression_level))+
  geom_boxplot()+
  theme_classic()+
  ggtitle('Human Protein Atlas RNA for HSP90AB1')


ggplot(data = hsp_hpa_prot, aes(x = Level, y = log10(Turnover_days)))+
  geom_boxplot()+
  theme_classic()+
  ggtitle('Human Protein Atlas protein for HSP90AB1')

dev.off()

rna_lm <- lm(hsp90_hpa_rna$expression_level ~ hsp90_hpa_rna$turnover_rate_days)
summary(rna_lm)
#Call:
#  lm(formula = hsp_hpa_rna$expression_level ~ hsp_hpa_rna$turnover_rate_days)

#Residuals:
#  Min      1Q  Median      3Q     Max 
#-52.494 -22.047 -11.866   6.097 136.663 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                    86.732505   8.461983  10.250 2.33e-12 ***
#  hsp_hpa_rna$turnover_rate_days  0.001115   0.000561   1.987   0.0544 .  
#---
#  Signif. codes:  0 â€˜***â€™ 0.001 â€˜**â€™ 0.01 â€˜*â€™ 0.05 â€˜.â€™ 0.1 â€˜ â€™ 1

#Residual standard error: 41.47 on 37 degrees of freedom
#(22 observations deleted due to missingness)
#Multiple R-squared:  0.09642,	Adjusted R-squared:  0.072 
#F-statistic: 3.948 on 1 and 37 DF,  p-value: 0.05436



ggplot(data = na.omit(hsp90_hpa_rna), aes(x = (log10(turnover_rate_days) > 3) , y = expression_level))+
  geom_boxplot()+
  theme_classic()

wilcox.test(hsp90_hpa_rna$expression_level ~ (log10(hsp90_hpa_rna$turnover_rate_days) > 3))
#Wilcoxon rank sum test

#data:  hsp_hpa_rna$expression_level by log10(hsp_hpa_rna$turnover_rate_days) > 3
#W = 79, p-value = 0.001728
#alternative hypothesis: true location shift is not equal to 0

fast <- na.omit(hsp90_hpa_rna[hsp90_hpa_rna$turnover_rate_days <= 31,'tissue'])
intermediate <- na.omit(hsp90_hpa_rna[hsp90_hpa_rna$turnover_rate_days > 31 & hsp90_hpa_rna$turnover_rate_days <= 3650, 'tissue'])
slow <- na.omit(hsp90_hpa_rna[hsp90_hpa_rna$turnover_rate_days > 3650, 'tissue'])

pdf('../../Body/4_Figures/hsp90.human.protein.atlas.tissues.turover.rate.explession.level.3.boxplots.pdf', width = 8, height = 6)
par(mar=c(5,6,3,1))
boxplot(hsp90_hpa_rna[hsp90_hpa_rna$tissue %in% fast, 'expression_level'], hsp90_hpa_rna[hsp90_hpa_rna$tissue %in% intermediate, 'expression_level'] ,
        hsp90_hpa_rna[hsp90_hpa_rna$tissue %in% slow, 'expression_level'], col = 'cyan3', names = c('fast', 'intermediate', 'slow'),
        ylab = 'Expression level of HSP90, NX', xlab = 'turnover rate', outline = F, ylim = c(40, 145), cex.lab = 2.1, cex.axis = 1.7)

dev.off()




hsc70_hsp_rna <- read.table('../../Body/2_Derived/human_protein_atlas_hsc70_expression_level.csv', header = T, sep = ',')
hsc70_hsp_rna <- merge(hsc70_hsp_rna,hsp90_hpa_rna[,c(2,4)], by = 'tissue')

ggplot(data = hsc70_hsp_rna, aes(log10(turnover_rate_days), expression_level, label=tissue))+
  geom_point()+
  geom_text(check_overlap = TRUE, vjust = -0.2)+
  geom_smooth(method = lm)+
  ggtitle('Human Protein Atlas RNA for HSPA8')+
  theme_classic()


rna_lm <- lm(hsc70_hsp_rna$expression_level ~ hsc70_hsp_rna$turnover_rate_days)
summary(rna_lm)
#Call:
#  lm(formula = hsc70_hsp_rna$expression_level ~ hsc70_hsp_rna$turnover_rate_days)

#Residuals:
#  Min      1Q  Median      3Q     Max 
#-60.539 -24.682  -9.966   3.642 130.217 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                      9.171e+01  8.405e+00  10.911 4.04e-13 ***
#  hsc70_hsp_rna$turnover_rate_days 4.853e-04  5.572e-04   0.871    0.389    
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 41.19 on 37 degrees of freedom
#(22 observations deleted due to missingness)
#Multiple R-squared:  0.0201,	Adjusted R-squared:  -0.006387 
#F-statistic: 0.7588 on 1 and 37 DF,  p-value: 0.3893


boxplot(hsc70_hsp_rna[hsc70_hsp_rna$tissue %in% fast, 'expression_level'], hsc70_hsp_rna[hsc70_hsp_rna$tissue %in% intermediate, 'expression_level'] ,
        hsc70_hsp_rna[hsc70_hsp_rna$tissue %in% slow, 'expression_level'], col = 'cyan3', names = c('fast', 'intermediate', 'slow'),
        ylab = 'Expression level of HSC70, NX', xlab = 'turnover rate', outline = F, ylim = c(30, 145), cex.lab = 1.5, cex.axis = 1.2)









