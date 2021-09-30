rm(list=ls(all=TRUE))

data <- read.table('../../Body/2_Derived/kn.ks.genes.like.hsp.txt')
data[data == 'n/a'] <- NA 

df <- c()

genes <- colnames(data)[-c(1,302)]

for (i in genes){
  gene <-  gsub('dN.dS_', '',i)
  a <- na.omit(data[,c('Species', i, 'Generation_Length')])
  colnames(a) <- c('Species', 'kn.ks', 'Generation_Length')
  Gene_id <- rep(gene, dim(a)[1])
  a <- cbind(a, Gene_id)
  df <- rbind(df,a)
}

df[, -c(1,4)] <- sapply(df[, -1], as.numeric)
library('lsmeans')
model <- lm(kn.ks ~ Generation_Length*Gene_id, data = df)
anova(model)
lst <- lstrends(model, "Gene_id", var="Generation_Length")
pairs.model <- pairs(lst)
fff <- as.data.frame(pairs.model)

write.table(pairs.model, '../../Body/2_Derived/kn.ks.hsp.like.genes.pair.comparison.txt')






write.table(df, '../../Body/2_Derived/hsp.like.genes.Species.kn.ks.Generation.Length.Gene.id.txt')


library(ggplot2)
pdf('../../Body/4_Figures/somegraff.pdf')
ggplot(df, aes(x=Generation_Length, y=kn.ks, col=Gene_id, alpha = 0.1), show.legend = FALSE, alpha = 0.1) + 
  geom_smooth(method="lm", se=FALSE , show.legend = FALSE, alpha = 0.1)
dev.off()

