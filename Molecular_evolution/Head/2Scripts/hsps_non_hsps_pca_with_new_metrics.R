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
hsc70_int <- read.csv('../../Body/2_Derived/hsc70_clients.Elisa_more_0.2.WT_interaction_score.csv', header = T)

genes_PC <- merge(genes_PC, hsp90_int, by.x = 'gene_ids', by.y = 'EnsemlID', all.x = T)
genes_PC <- merge(genes_PC, hsc70_int, by.x = 'gene_ids', by.y = 'EnsemlID', all.x = T)
colnames(genes_PC) <- c('ID', 'PC1', 'PC2', 'WT_interaction_score_hsp90', 'WT_interaction_score_hsc70')

hsps <- read.csv('../../Body/2_Derived/human.hsp.ensemblid.list.txt', header = F)
hsps <- hsps$V1
groups <- read.csv('../../Body/2_Derived/human.hsp.ensID.group.txt', sep = '\t')
genes_PC <- merge(genes_PC, groups, by.x = 'ID', by.y = 'Ensembl.gene.ID', all.x = T)


length(na.omit(genes_PC$WT_interaction_score_hsc70))

dim(genes_PC[!is.na(genes_PC$WT_interaction_score_hsp90) & (genes_PC$WT_interaction_score_hsp90 >2.5),]) #133
dim(genes_PC[!is.na(genes_PC$WT_interaction_score_hsc70) & (genes_PC$WT_interaction_score_hsc70 >2.5),]) #118
autoplot(genes_pca, data = genes_without_na, colour = 'gray', loadings = TRUE, loadings.label = TRUE, loadings.label.size = 3, scale = 0, 
         loadings.colour = 'black', loadings.label.colour = 'black')+
  #geom_point(data = both, aes(PC1, PC2, color = WT_interaction_score_hsp90, alpha = 0.3),size = 2)+
  #geom_point(data = both, aes(PC1, PC2, color = WT_interaction_score_hsc70, alpha = 0.1),size = 1)+
  ggtitle('Principal Component analisys with clients of hsp90 and hsc70')+
  geom_point(data = genes_PC[genes_PC$ID %in% hsps,], aes(PC1, PC2, size = 1, alpha = 0.5, color = Group.name))+
  theme_bw()
  #theme(legend.position = 'Bottom')



for_paml <- genes_PC[(genes_PC$ID %in% hsps) & (genes_PC$PC1 > quantile(genes_PC$PC1, probs = 0.85)),]$ID
files <- list.files('../../../Kn.Ks.with.GODON/Body/1_Raw/omm_NT_fasta.v10b_116tax_CDS_final/')
files_hsp <- files[gsub('_[A-Z]*.*', '', files) %in% for_paml]
for_paml <- gsub('_[A-Z]*.*', '', files_hsp)

dist <- function(a_1,a_2, b_1, b_2){
  return(sqrt((a_1 - a_2)^2 + (b_1 - b_2)^2))
} 

genes_PC <- genes_PC[genes_PC$ID %in% gsub('_[A-Z]*.*', '', files),]

df_dist <- setNames(data.frame(matrix(ncol = 16, nrow = 10183)), c('non_hsp', for_paml))
df_dist$non_hsp <- genes_PC[!(genes_PC$ID %in% for_paml), 1]

f <- function(x, y){
pc1x <- genes_PC[genes_PC$ID == x, 'PC1']
pc2x <- genes_PC[genes_PC$ID == x, 'PC2']

pc1y <- genes_PC[genes_PC$ID == y, 'PC1']
pc2y <- genes_PC[genes_PC$ID == y, 'PC2']

d <-dist(pc1x, pc1y, pc2x, pc2y)
return(d)
}



for (hsp in colnames(df_dist)[-1]){
  d <- lapply(df_dist$non_hsp, f, y=hsp)
  for (i in d[lengths(d) > 1]) {d[lengths(d) > 1] <- i[1]}
  df_dist[,hsp] <- unlist(d, use.names=FALSE)
}

control <- c()
for (h in colnames(df_dist[,-1])){
  ord <- df_dist[order(df_dist[,h]), h]
  min_1 <- df_dist[which(df_dist[,h] == ord[1]),1]
  if (min_1 %in% control){
    min_2 <- df_dist[which(df_dist[,h] == ord[2]),1]
    control <- c(control, min_2)
    print(paste(h, min_2, sep = '  '))
  }
  else{control <- c(control, min_1)
  print(paste(h, min_1, sep = '  '))}
  
}

control_files <- files[gsub('_[A-Z]*.*', '', files) %in% control]

hsps_with_control <- data.frame(hsp <- colnames(df_dist[,-1]), control_gene <- control)
colnames(hsps_with_control) <- c('hsp', 'control_gene')

write.table(hsps_with_control, '../../../Kn.Ks.with.GODON/Body/2_Derived/hsps.with.control.nonhsps.txt')

write(files_hsp, '../../../Kn.Ks.with.GODON/Body/2_Derived/hsp_files_list.txt')

write(control_files, '../../../Kn.Ks.with.GODON/Body/2_Derived/nonhsp_files_list.txt')

