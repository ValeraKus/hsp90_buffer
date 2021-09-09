rm(list=ls(all=TRUE))
library(ggfortify)

pdf("../../Body/4_Figures/GTEx.All.sign.pairs.all.tissues.Haploins.Analyses.pdf")

qtl = read.table("../../Body/2_Derived/Gtex.All.sign.pairs.all.tissues.Haploins.txt", header = TRUE)
genes <- qtl$EnsemblId

qtl = qtl[,c(2:8)]; str(qtl) # keep only numeric
names(qtl)[7] = c('Haploins')
qtl[is.na(qtl)] <- 0
round(cor(qtl),2) # correlation matrix
plot(qtl)         # all pairwise plots


QtlPca <- prcomp(qtl, scale = TRUE)
print(QtlPca) 
summary(QtlPca)
a1 <- QtlPca$rotation[, 1]
a1

PC1 <- predict(QtlPca)[,1]
PC2 <- predict(QtlPca)[,2]
qtl_pred <- data.frame(genes, PC1, PC2)
plot(QtlPca)
biplot(QtlPca, col = c("gray", "black"))

qtl_pred$name <- ''
qtl_pred[qtl_pred$genes == "ENSG00000096384", 'name'] = "HSP90AB1"


autoplot(QtlPca, data = qtl, colour = 'gray', loadings = TRUE, loadings.label = TRUE, loadings.label.size = 3, scale = 0, 
         loadings.colour = 'black', loadings.label.colour = 'black')+
  geom_point(aes(qtl_pred[qtl_pred$genes == "ENSG00000096384", "PC1"], qtl_pred[qtl_pc$genes == "ENSG00000096384", "PC2"]), colour = 'red', size = 3)+
  geom_text(data=qtl_pred[qtl_pred$genes == 'ENSG00000096384',], aes(PC1,PC2+1, label = name))

dev.off()



boxplot(qtl$NumberGoeEqtl/qtl$NumberEqtl, ylab = 'Доля GOE вариантов')
points(qtl[qtl$EnsemblId == 'ENSG00000096384', 'NumberGoeEqtl']/qtl[qtl$EnsemblId == 'ENSG00000096384', 'NumberEqtl'], col = 'red', pch = 19, cex = 1.5)
legend('topright', legend = 'HSP90', col = 'red', pch = 19, cex = 1.5)




boxplot(qtl$NumberLoeEqtl/qtl$NumberEqtl)
points(qtl[qtl$EnsemblId == 'ENSG00000096384', 'NumberLoeEqtl']/qtl[qtl$EnsemblId == 'ENSG00000096384', 'NumberEqtl'], col = 'red', pch = 19, cex = 1.5)
legend('topright', legend = 'HSP90', col = 'red', pch = 19, cex = 1.5)

boxplot(qtl$NumberLoeEqtl/qtl$NumberGoeEqtl)
points(qtl[qtl$EnsemblId == 'ENSG00000096384', 'NumberLoeEqtl']/qtl[qtl$EnsemblId == 'ENSG00000096384', 'NumberGoeEqtl'], col = 'red', pch = 19, cex = 1.5)
legend('topright', legend = 'HSP90', col = 'red', pch = 19, cex = 1.5)

boxplot(qtl$NumberEqtl, ylim = c(0, 5000))
points(qtl[qtl$EnsemblId == 'ENSG00000096384', 'NumberEqtl'], col = 'red', pch = 19, cex = 1.5)
legend('topright', legend = 'HSP90', col = 'red', pch = 19, cex = 1.5)


hsp_goe_prop <- qtl[qtl$EnsemblId == 'ENSG00000096384', 'NumberGoeEqtl']/qtl[qtl$EnsemblId == 'ENSG00000096384', 'NumberEqtl']
qtl_goe_prop <- qtl$NumberGoeEqtl/qtl$NumberEqtl
sum(qtl_goe_prop >= hsp_goe_prop, na.rm = T)/length(qtl_goe_prop)


hsp_like_genes = read.table('../../Body/2_Derived/hsp.like.genes.pca.genes.ranged.by.distance.to.hsp.txt')
hsp_like_genes <- hsp_like_genes[1:300,]
qtl_hsp_like <- qtl[(qtl$EnsemblId %in% hsp_like_genes$EnsemblId),]

boxplot(qtl_hsp_like$NumberGoeEqtl/qtl_hsp_like$NumberEqtl)
points(qtl_hsp_like[qtl_hsp_like$EnsemblId == 'ENSG00000096384', 'NumberGoeEqtl']/qtl_hsp_like[qtl_hsp_like$EnsemblId == 'ENSG00000096384', 'NumberEqtl'], col = 'red', pch = 19, cex = 1.5)
legend('topright', legend = 'HSP90', col = 'red', pch = 19, cex = 1.5)




