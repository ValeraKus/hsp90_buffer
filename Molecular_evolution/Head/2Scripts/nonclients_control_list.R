rm(list = ls(all = TRUE))

genes_PC <- read.table('../../Body/2_Derived/pca_with_new_metrics_hsp90_hsc70_clients.txt')



hsp90_int <- read.csv('../../Body/2_Derived/hsp90_clients.Elisa_more_0.2.WT_interaction_score.csv', header = T)
hsc70_int <- read.csv('../../Body/2_Derived/hsc70_clients.Elisa_more_0.2.WT_interaction_score.csv', header = T)
colnames(hsp90_int) <- c('gene_ids', 'WT_int_score_hsp90')
colnames(hsc70_int) <- c('gene_ids', 'WT_int_score_hsc70')

genes_PC <- merge(genes_PC, hsp90_int, by = 'gene_ids', all.x = T)
genes_PC <- merge(genes_PC, hsc70_int, by = 'gene_ids', all.x = T)

hsp90_clients <- genes_PC[!is.na(genes_PC$WT_int_score_hsp90) & (genes_PC$WT_int_score_hsp90 > 3), 'gene_ids'] #120
hsc70_clients <- genes_PC[!is.na(genes_PC$WT_int_score_hsc70) & (genes_PC$WT_int_score_hsc70 > 3), 'gene_ids'] #90

clients <- c(hsp90_clients, hsc70_clients)
clients <- clients[!duplicated(clients)] #141



files <- list.files('../../../Kn.Ks.with.GODON/Body/1_Raw/omm_NT_fasta.v10b_116tax_CDS_final/')
files_clients <- files[gsub('_[A-Z]*.*', '', files) %in% clients] #124
clients <- clients[clients %in% gsub('_[A-Z]*.*', '', files)]



dist <- function(a_1,a_2, b_1, b_2){
  return(sqrt((a_1 - a_2)^2 + (b_1 - b_2)^2))
}  
  
genes_PC <- genes_PC[genes_PC$gene_ids %in% gsub('_[A-Z]*.*', '', files),]
genes_PC <- genes_PC[!duplicated(genes_PC$gene_ids),]

df_dist <- setNames(data.frame(matrix(ncol = 125, nrow = 10030)), c('nonclients', clients))
df_dist$nonclients <- genes_PC[!(genes_PC$gene_ids %in% clients), 1]

f <- function(x, y){
  pc1x <- genes_PC[genes_PC$gene_ids == x, 'PC1']
  pc2x <- genes_PC[genes_PC$gene_ids == x, 'PC2']
  
  pc1y <- genes_PC[genes_PC$gene_ids == y, 'PC1']
  pc2y <- genes_PC[genes_PC$gene_ids == y, 'PC2']
  
  d <-dist(pc1x, pc1y, pc2x, pc2y)
  return(d)
}



for (client in colnames(df_dist)[-1]){
  d <- lapply(df_dist$nonclients, f, y=client)
  for (i in d[lengths(d) > 1]) {d[lengths(d) > 1] <- i[1]}
  df_dist[,client] <- unlist(d, use.names=FALSE)
}

control <- c()
for (h in colnames(df_dist[,-1])){
  ord <- df_dist[order(df_dist[,h]), h]
  min_1 <- df_dist[which(df_dist[,h] == ord[1]),1]
  if (min_1 %in% control){
    min_2 <- df_dist[which(df_dist[,h] == ord[2]),1]
    control <- c(control, min_2)
  }
  else{control <- c(control, min_1)}
  
}
control_files <- files[gsub('_[A-Z]*.*', '', files) %in% control]
clients_with_control <- data.frame(colnames(df_dist[,-1]), control)
colnames(clients_with_control) <- c('client', 'control_gene')
clients_with_control$hsp90_client <- ifelse(clients_with_control$client %in% hsp90_clients, T, F)
clients_with_control$hsc70_client <- ifelse(clients_with_control$client %in% hsc70_clients, T, F)
write.table(clients_with_control, '../../../Kn.Ks.with.GODON/Body/2_Derived/clients.with.control.nonclients.txt')


write(files_clients, '../../../Kn.Ks.with.GODON/Body/2_Derived/clients_files_list.txt')
write(control_files, '../../../Kn.Ks.with.GODON/Body/2_Derived/nonclients_files_list.txt')


























plot(genes_PC$PC1, genes_PC$PC2, col = 'grey')
points(clients_PC$PC1, clients_PC$PC2, col = 'blue', pch = 22)


length(genes_PC[genes_PC$client, 1]) #316
length(genes_PC[!genes_PC$client, 1]) #11187

df_dist <- setNames(data.frame(matrix(ncol = 317, nrow = 11187)), c('nonclients_id', genes_PC[genes_PC$client, 1]))
df_dist$nonclients_id <- genes_PC[!genes_PC$client, 1]


f <- function(x, y){
  pc1x <- genes_PC[genes_PC$gene_ids == x, 'PC1']
  pc2x <- genes_PC[genes_PC$gene_ids == x, 'PC2']
  
  pc1y <- genes_PC[genes_PC$gene_ids == y, 'PC1']
  pc2y <- genes_PC[genes_PC$gene_ids == y, 'PC2']
  
  d <-dist(pc1x, pc1y, pc2x, pc2y)
  return(d)
}



for (client in colnames(df_dist)[-1]){
  d <- lapply(df_dist$nonclients_id, f, y=client)
  for (i in d[lengths(d) > 1]) {d[lengths(d) > 1] <- i[1]}
  df_dist[,client] <- unlist(d, use.names=FALSE)
  }

write.table(df_dist, 'distance.between.clients.and.nonclients.in.pca.txt')  
  
closest_nonclient <- function(id) {
  #min_1 <- (df_dist[which.min(df_dist[,id]),1])
  ord <- df_dist[order(df_dist[,id]), id]
  min_1 <- df_dist[which(df_dist[,id] == ord[1]),1]
  min_2 <- df_dist[which(df_dist[,id] == ord[2]),1]
  return(c(min_1, min_2))
}
  














neighbours <- unlist(lapply(colnames(df_dist[,-1]), FUN = closest_nonclient))  
neighbours <- neighbours[!duplicated(neighbours)]  
write(neighbours, '../../Body/2_Derived/nonclients_control_list.txt', sep = '\n')  
  


plot(genes_PC$PC1, genes_PC$PC2, col ='grey')
points(clients_PC$PC1, clients_PC$PC2, col = 'blue', pch = 19)
points(genes_PC[genes_PC$gene_ids %in% neighbours, 'PC1'], genes_PC[genes_PC$gene_ids %in% neighbours, 'PC2'], col = rgb(red=0.5, green=1, blue=0, alpha=0.5), pch = 8)















  