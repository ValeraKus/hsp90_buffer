rm(list=ls(all=TRUE))

df <- read.table('../../Body/3_Results/hsp.like.genes.linear.model.results.kn.ks.vs.generation.length.mammals.right.way.txt')


pdf('../../Body/4_Figures/.hsp.like.genes.linear.model.parameters.dist.pdf')
par(mfrow = c(1,2))
hist(df$slopes)
hist(df$p_val_slope)
hist(df$intercept)
hist(df$p_val_intercept)
par(mfrow = c(1,1))
hist(df$number_of_species)
par(mfrow = c(1,2))
hist(df$R_sq)
hist(df$R_sq_adj)
par(mfrow = c(1,1))
hist(df$residual_std_err)
dev.off()
