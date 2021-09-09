rm(list=ls(all=TRUE))

df <- read.csv('../../Body/1_Raw/hsp90_human_tissues_turnover_rate_expression_level.csv')

pdf('../../Body/4_Figures/hsp90.human.tissues.turnover.rate.expression.level.pdf')
plot(log10(df$Turnover..days.), log10(df$Expression.level.ofHSP90AB1.in.human..TPM..median.), col = 'blue', pch = 19)
text(log10(df$Turnover..days.)+0.1, log10(df$Expression.level.ofHSP90AB1.in.human..TPM..median.)+0.02, labels=df$Tissue, cex= 0.7)


plot(df$Turnover..days., df$Expression.level.ofHSP90AB1.in.human..TPM..median., col = 'blue', pch = 19)
text(df$Turnover..days.+500, df$Expression.level.ofHSP90AB1.in.human..TPM..median.+32, labels=df$Tissue, cex= 0.7)
dev.off()

cor.test(df$Turnover..days., df$Expression.level.ofHSP90AB1.in.human..TPM..median.)
#Pearson's product-moment correlation

#data:  df$Turnover..days. and df$Expression.level.ofHSP90AB1.in.human..TPM..median.
#t = 2.0158, df = 26, p-value = 0.05427
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# -0.006299571  0.651376489
#sample estimates:
#      cor 
#0.3676414 

cor.test(log10(df$Turnover..days.), log10(df$Expression.level.ofHSP90AB1.in.human..TPM..median.))
#	Pearson's product-moment correlation

#data:  log10(df$Turnover..days.) and log10(df$Expression.level.ofHSP90AB1.in.human..TPM..median.)
#t = 0.91171, df = 26, p-value = 0.3703
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  -0.2109170  0.5152526
#sample estimates:
#  cor 
#0.1760098 