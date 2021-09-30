rm(list=ls(all=TRUE))
library(ggfortify)
library(ggrepel)


hsps <- read.csv('../../Body/2_Derived/human.hsp.ensemblid.list.txt', header = F)
hsps <- as.vector(hsps$V1)


genes <- read.table('../../Body/1_Raw/gencode.v25.annotation.gtf.Genes.Shet.pLI.FIS.RVIS.GHIS.KnKs.GC.BrainSpecificRanking.Branch', header = 1)

genes_hsp <- genes[genes$EnsemblId %in% hsps, ]

numer <- genes_hsp[,c(1,2,6,8,9,10,11,12,13,14,15,16,17,19)]

num_withoun_na <- na.omit(numer)

hsp_pca <- prcomp(num_withoun_na[,-c(1,2)], scale. = T)

summary(hsp_pca)
print(hsp_pca)
a1 <- hsp_pca$rotation[, 1]
a1

PC1 <- predict(hsp_pca)[,1]
PC2 <- predict(hsp_pca)[,2]
qtl_pred <- data.frame(num_withoun_na$EnsemblId, PC1, PC2)
plot(hsp_pca)
biplot(hsp_pca, col = c("gray", "black"))

groups <- read.csv('../../Body/2_Derived/human.hsp.ensID.group.txt', sep = '\t')
qtl_pred <- merge(x = qtl_pred, y = groups, by.x = 'num_withoun_na.EnsemblId', by.y = 'Ensembl.gene.ID')


autoplot(hsp_pca, data = num_withoun_na[,-c(1,2)], colour = 'gray', loadings = TRUE, loadings.label = TRUE, loadings.label.size = 3, scale = 0, 
         loadings.colour = 'black', loadings.label.colour = 'black')+
  geom_point(aes(qtl_pred$PC1, qtl_pred$PC2, color = qtl_pred$Group.name))+
  geom_label_repel(data=qtl_pred[qtl_pred$num_withoun_na.EnsemblId == 'ENSG00000096384',], aes(PC1,PC2, label = name),
                   box.padding   = 0.6, 
                   point.padding = 0.8,
                   segment.color = 'grey50')+
  theme_bw()

qtl_pc