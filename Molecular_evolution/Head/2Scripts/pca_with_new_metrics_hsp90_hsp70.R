rm(list = ls(all = TRUE))
library(ggfortify)
library(dabestr)
set.seed(42)

metrics_table <- read.table('../../Body/2_Derived/gencode.v25.annotation.gtf.Genes.Shet.pLI.FIS.RVIS.GHIS.KnKs.GC.BrainSpecificRanking.Branch.gnomad.oe_lof_upper_bin.p', header = T)

genes_num_cols <- metrics_table[metrics_table$GeneType == 'protein_coding', c('EnsemblId', 'GeneLength', 'NumberOfTranscripts', 'AverageExonLength', 'AverageNumberOfExonsPerTranscript', 'SelCoeffHet',
                                                                              'ProbOfBeingLofIntolerant', 'FunctionalIndispensabilityScore', 'ResidualVariationIntoleranceScore', 'GenomeWideHaploinsufficiencyScore',
                                                                              'KnKsMouse', 'GcContent', 'Branch', 'oe_lof_upper_bin', 'p')]

summary(genes_num_cols)[7,]
genes_without_na <- na.omit(genes_num_cols) #11503


gene_ids <- genes_without_na$EnsemblId
genes_without_na <- genes_without_na[,-1]



genes_pca <- prcomp(genes_without_na, scale. = T)
plot(genes_pca)

PC1 <- predict(genes_pca)[,1]
PC2 <- predict(genes_pca)[,2]
PC1_var <- summary(genes_pca)$importance[2]
PC2_var <- summary(genes_pca)$importance[5]


genes_PC <- data.frame(gene_ids, PC1, PC2)

hsp90_int <- read.csv('../../Body/2_Derived/hsp90_clients.Elisa_more_0.2.WT_interaction_score.csv', header = T)







write.table(genes_PC, '../../Body/2_Derived/pca_with_new_metrics_hsp90_hsc70_clients.txt')

clients <- read.csv('../../Body/2_Derived/Chaperone-Protein_Interaction_hsp90_hsc70_clients.csv')

clients_PC <- merge(x = genes_PC, y = clients, by.x = 'gene_ids', by.y = 'clients')
clients_PC$chaperone <- ifelse((clients_PC$HSP90AB1 == 'True') & (clients_PC$HSPA8 == 'True'), 'both', NaN)
clients_PC[(clients_PC$HSP90AB1 == 'True') & (clients_PC$HSPA8 == 'False'), 'chaperone'] = 'HSP90AB1'
clients_PC[(clients_PC$HSP90AB1 == 'False') & (clients_PC$HSPA8 == 'True'), 'chaperone'] = 'HSPA8'






pdf('../../Body/4_Figures/pca_with_new_metrics_hsp90_hsc70_clients.pdf')
autoplot(genes_pca, data = genes_without_na, colour = 'gray', loadings = TRUE, loadings.label = TRUE, loadings.label.size = 3, scale = 0, 
         loadings.colour = 'black', loadings.label.colour = 'black')+
  geom_point(data = clients_PC, aes(PC1, PC2, color = chaperone, alpha = 0.5),size = 2)+
  ggtitle('Principal Component analisys with clients of hsp90 and hsc70')+
  geom_point(data = genes_PC[genes_PC$gene_ids == 'ENSG00000096384',], aes(PC1, PC2, size = 1, alpha = 0.5, color = 'hsp90'))+
  geom_point(data = genes_PC[genes_PC$gene_ids == 'ENSG00000109971',], aes(PC1, PC2, size = 1, alpha = 0.5, color = 'hsp70'))+         
  theme_bw()+
  scale_color_manual(name  ="Genes",
                     breaks=c("both", "hsp70", 'hsp90', 'HSP90AB1', 'HSPA8'),
                     labels=c("HSP90 and HSC70 clients", "HSPA8 (HSC70)", 'HSP90AB1 (HSP90)', 'HSP90 clients', 'HSC70 clients'),
                     values = c('indianred2', 'goldenrod1', 'red3', 'deeppink4', 'purple'))+
  
  guides(alpha = FALSE, size = FALSE)




#plotting density plots for pc1
genes_PC$client <- ifelse(!(genes_PC$gene_ids %in% clients_PC$gene_ids), 'No', 'Yes')
genes_PC[genes_PC$client == 'Yes',]$client <- ifelse((genes_PC[genes_PC$client == 'Yes',]$gene_ids %in% clients_PC[(clients_PC$HSP90AB1 == 'True') & (clients_PC$HSPA8 == 'False'),]$gene_ids), 'hsp90', 'Yes')
genes_PC[genes_PC$client == 'Yes',]$client <- ifelse((genes_PC[genes_PC$client == 'Yes',]$gene_ids %in% clients_PC[(clients_PC$HSP90AB1 == 'False') & (clients_PC$HSPA8 == 'True'),]$gene_ids), 'hsc70', 'Yes')
genes_PC[genes_PC$client == 'Yes',]$client <- ifelse((genes_PC[genes_PC$client == 'Yes',]$gene_ids %in% clients_PC[(clients_PC$HSP90AB1 == 'True') & (clients_PC$HSPA8 == 'True'),]$gene_ids), 'both', 'Yes')


par(mfcol = c(2,1))
g1 <- ggplot(genes_PC, aes(x=PC1, color=client)) +
  geom_density()+
  scale_color_manual(values=c('indianred2', 'deeppink4', 'purple', 'darkgray'))+
  theme_bw()
ggsave(g1, file = '../../Body/4_Figures/density_pc1_clients.png', height = 2, width = 13)



g2 <- ggplot(genes_PC, aes(x=PC2, color=client)) +
  geom_density()+
  scale_color_manual(values=c('indianred2', 'deeppink4', 'purple', 'darkgray'))+
  theme_bw()
ggsave(g2, file = '../../Body/4_Figures/density_pc2_clients.png', height = 2, width = 13)

dev.off()



  

##################### Mann_Whitney U test#############

#for PC1
not_clients <- genes_PC[genes_PC$client == 'No', 'PC1']
hsp90 <- genes_PC[(genes_PC$client == 'hsp90') | (genes_PC$client == 'both'), 'PC1']
hsc70 <- genes_PC[(genes_PC$client == 'hsc70') | (genes_PC$client == 'both'), 'PC1']
hsp90_only <- genes_PC[(genes_PC$client == 'hsp90'), 'PC1']
hsc70_only <- genes_PC[(genes_PC$client == 'hsc70'), 'PC1']
both <- genes_PC[(genes_PC$client == 'both'), 'PC1']

#hsp90 clients
wilcox.test(not_clients, hsp90)
#W = 1406300, p-value = 0.01111


#hsc70
wilcox.test(not_clients, hsc70)
#W = 1560600, p-value = 0.000487


#both
wilcox.test(not_clients, both)
#W = 994380, p-value = 0.01567


#only hsp90
wilcox.test(not_clients, hsp90_only)

#only hsc70
wilcox.test(not_clients, hsc70_only)


#for PC2
not_clients <- genes_PC[genes_PC$client == 'No', 'PC2']
hsp90 <- genes_PC[(genes_PC$client == 'hsp90') | (genes_PC$client == 'both'), 'PC2']
hsc70 <- genes_PC[(genes_PC$client == 'hsc70') | (genes_PC$client == 'both'), 'PC2']
hsp90_only <- genes_PC[(genes_PC$client == 'hsp90'), 'PC2']
hsc70_only <- genes_PC[(genes_PC$client == 'hsc70'), 'PC2']
both <- genes_PC[(genes_PC$client == 'both'), 'PC2']

#hsp90 clients
wilcox.test(not_clients, hsp90)
#W = 1336000, p-value = 0.2646

#hsc70
wilcox.test(not_clients, hsc70)
#W = 1428400, p-value = 0.3612

#both
wilcox.test(not_clients, both)
#W = 942960, p-value = 0.2434

#only hsp90
wilcox.test(not_clients, hsp90_only)

#only hsc70
wilcox.test(not_clients, hsc70_only)



shapiro.test(clients_PC$PC1)#p-value = 2.009e-12
wilcox.test(clients_PC$PC1, mu =0)
#V = 18136, p-value = 2.147e-05

shapiro.test(clients_PC$PC2)#p-value = 4.37e-08
wilcox.test(clients_PC$PC2, mu =0)
#V = 22341, p-value = 0.09651




library(reshape2)
data_long <- melt(clients_PC[,c(1,2,3)], id.vars = c('gene_ids'))

ggplot(data = data_long, aes(x = variable, y = value, fill = variable))+
  geom_violin()+
  theme_bw()+
  ggtitle('hsp90 and hsc70 clients')+
  geom_boxplot(width=0.1)+
  geom_hline(yintercept=00, color = 'gray')



shared.control <- 
  genes_PC %>%
  dabest(client, PC1, 
         idx = c("No", "hsp90", "hsc70"),
         paired = FALSE
  )
shared.control.mean_diff <- shared.control %>% mean_diff()
plot(unpaired_mean_diff)

two.group.unpaired <- 
  genes_PC %>%
  dabest(client, PC1,
         idx = c("No", "hsp90"), 
         paired = FALSE)
two.group.unpaired 
two.group.unpaired.meandiff <- mean_diff(two.group.unpaired, reps = 7000)



both <- sample(clients_PC[clients_PC$chaperone == 'both',]$gene_ids, 15)
hsp90 <- sample(clients_PC[clients_PC$chaperone == 'HSP90AB1',]$gene_ids, 15)
hsc70 <- sample(clients_PC[clients_PC$chaperone == 'HSPA8',]$gene_ids, 15)

dist <- read.table('distance.between.clients.and.nonclients.in.pca.txt')
dist_both <- cbind(dist$nonclients_id, dist[,(colnames(dist[,-1]) %in% both)])
dist_hsp90 <- cbind(dist$nonclients_id,dist[,(colnames(dist) %in% hsp90)])
dist_hsc70 <- cbind(dist$nonclients_id, dist[,(colnames(dist) %in% hsc70)])

closest_nonclient <- function(id) {
  ord <- dist_hsc70[order(dist_hsc70[,id]), id]
  return(dist_hsc70[which(dist_hsc70[,id] %in% ord[1:2]),1])
}
control_both <- unlist(lapply(colnames(dist_both[,-1]), FUN = closest_nonclient))
control_both <- control_both[!duplicated(control_both)] 

control_hsp90 <- unlist(lapply(colnames(dist_hsp90[,-1]), FUN = closest_nonclient))
control_hsp90 <- control_hsp90[!duplicated(control_hsp90)]

control_hsc70 <- unlist(lapply(colnames(dist_hsc70[,-1]), FUN = closest_nonclient))
control_hsc70 <- control_hsc70[!duplicated(control_hsc70)]

write(both, '../../Body/2_Derived/clients_hsp90_hsc70_15.txt', sep = '\n') 
write(hsp90, '../../Body/2_Derived/clients_hsp90_only_15.txt', sep = '\n') 
write(hsc70, '../../Body/2_Derived/clients_hsc70_only_15.txt', sep = '\n') 

write(control_both, '../../Body/2_Derived/nonclients_contor_hsp90_hsc70_15.txt', sep = '\n') 
write(control_hsp90, '../../Body/2_Derived/nonclients_contor_hsp90_only_15.txt', sep = '\n')
write(control_hsc70, '../../Body/2_Derived/nonclients_contor_hsc70_only_15.txt', sep = '\n')



#######################################
hsp90_int <- read.csv('../../Body/2_Derived/hsp90_clients.Elisa_more_0.2.WT_interaction_score.csv', header = T)
hsp90_int <- merge(hsp90_int, genes_PC, by.x = 'EnsemlID', by.y = 'gene_ids', all.x = T)
cor.test(hsp90_int$WT_interaction_score, hsp90_int$PC1, method = 'spearman', exact=FALSE) # rho=0.05791891, p-value = 0.1439
cor.test(hsp90_int$WT_interaction_score, hsp90_int$PC2, method = 'spearman', exact=FALSE) # rho=-0.1376778, p-value = 0.0004878

lm_hsp90 <- lm(hsp90_int$WT_interaction_score ~ hsp90_int$PC1 + hsp90_int$PC2)
summary(lm_hsp90)
#Residuals:
#  Min      1Q  Median      3Q     Max 
#-5.4600 -1.8195 -0.8485  0.6194 19.1097 
#
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)    1.45925    0.12604  11.578  < 2e-16 ***
#  hsp90_int$PC1  0.26870    0.07619   3.526 0.000451 ***
#  hsp90_int$PC2 -0.33556    0.10444  -3.213 0.001380 ** 
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 3.088 on 635 degrees of freedom
#(202 observations deleted due to missingness)
#Multiple R-squared:  0.03159,	Adjusted R-squared:  0.02854 
#F-statistic: 10.36 on 2 and 635 DF,  p-value: 3.75e-05



hsc70_int <- read.csv('../../Body/2_Derived/hsc70_clients.Elisa_more_0.2.WT_interaction_score.csv', header = T)
hsc70_int <- merge(hsc70_int, genes_PC, by.x = 'EnsemlID', by.y = 'gene_ids', all.x = T)
cor.test(hsc70_int$WT_interaction_score, hsc70_int$PC1, method = 'spearman', exact=FALSE) # rho=0.04525303, p-value = 0.2609
cor.test(hsc70_int$WT_interaction_score, hsc70_int$PC2, method = 'spearman', exact=FALSE) # rho=-0.1114708, p-value = 0.005496

lm_hsc70 <- lm(hsc70_int$WT_interaction_score ~ hsc70_int$PC1 + hsc70_int$PC2)
summary(lm_hsc70)
#Residuals:
#  Min      1Q  Median      3Q     Max 
#-5.4263 -1.2024 -0.4265  0.9647  8.6665 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)    1.03998    0.07915  13.140   <2e-16 ***
#  hsc70_int$PC1  0.07082    0.04738   1.495   0.1355    
#hsc70_int$PC2 -0.13602    0.06685  -2.034   0.0423 *  
#  ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 1.915 on 616 degrees of freedom
#(198 observations deleted due to missingness)
#Multiple R-squared:  0.009355,	Adjusted R-squared:  0.006138 
#F-statistic: 2.908 on 2 and 616 DF,  p-value: 0.05531


both <- merge(hsp90_int, hsc70_int, by = 'EnsemlID')
cor.test(both$WT_interaction_score.x, both$WT_interaction_score.y, method = 'spearman')
###rho=0.6108015, p-value < 2.2e-16

both <- both[,c(1,2,5,6,7)]
colnames(both) <- c('id', 'WT_interaction_score_hsp90', 'WT_interaction_score_hsc70', 'PC1', 'PC2')



autoplot(genes_pca, data = genes_without_na, colour = 'gray', loadings = TRUE, loadings.label = TRUE, loadings.label.size = 3, scale = 0, 
         loadings.colour = 'black', loadings.label.colour = 'black')+
  #geom_point(data = both, aes(PC1, PC2, color = WT_interaction_score_hsp90, alpha = 0.3),size = 2)+
  #geom_point(data = both, aes(PC1, PC2, color = WT_interaction_score_hsc70, alpha = 0.1),size = 1)+
  ggtitle('Principal Component analisys with clients of hsp90 and hsc70')+
  geom_point(data = genes_PC[genes_PC$ID == 'ENSG00000096384',], aes(PC1, PC2, size = 1, alpha = 0.5, color = 'hsp90'))+
  geom_point(data = genes_PC[genes_PC$ID == 'ENSG00000109971',], aes(PC1, PC2, size = 1, alpha = 0.5, color = 'hsc70'))+
  geom_point(data = genes_PC[genes_PC$ID == 'ENSG00000144381',], aes(PC1, PC2, size = 1, alpha = 0.5, color = 'hsp60'))+
  geom_point(data = genes_PC[genes_PC$ID == 'ENSG00000080824',], aes(PC1, PC2, size = 1, alpha = 0.5, color = 'hsp90aa1'))+
  theme_bw()
  #scale_color_gradient(low = "blue", high = "red")

  
  #guides(alpha = FALSE, size = FALSE)





