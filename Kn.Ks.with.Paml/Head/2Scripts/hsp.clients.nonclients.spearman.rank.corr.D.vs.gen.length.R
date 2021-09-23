rm(list=ls(all=TRUE))


clients <-  read.table('../../Body/2_Derived/clients_godon_D.txt', header = T)
nonclients <- read.table('../../Body/2_Derived/nonclients_godon_D.txt', header = T)

cols <- colnames(clients[,-c(1, 335)])
for (c in cols) {
  if (length(clients[!is.na(clients[,c]),c]) <= 10){
    clients[,c] <- NULL
  }
}  


cols <- colnames(nonclients[,-c(1, 429)])
for (c in cols) {
  if (length(nonclients[!is.na(nonclients[,c]),c]) <= 10){
    nonclients[,c] <- NULL
  }
}


rho <- c()
S <- c()
p_val <- c()
number_of_species <- c()
genes <- colnames(clients)[-c(1,329)]



for (gene in genes){
  a <- na.omit(clients[,c("Species", gene, 'Generation_Length')])
  
  species <- a[,'Species']
  number_of_species <- c(number_of_species, length(species))
  
  rcor <- cor.test(~ a[,gene] + a[,'Generation_Length'], method = 'spearman', continuity = F, conf.level = 0.95)
  rho <- c(rho, rcor$estimate)
  S <- c(S, rcor$statistic)
  p_val <- c(p_val, rcor$p.value)
}


results_clients <- data.frame(genes, S, rho, p_val, number_of_species)


rho <- c()
S <- c()
p_val <- c()
number_of_species <- c()
genes <- colnames(nonclients)[-c(1,425)]



for (gene in genes){
  a <- na.omit(nonclients[,c("Species", gene, 'Generation_Length')])
  
  species <- a[,'Species']
  number_of_species <- c(number_of_species, length(species))
  
  rcor <- cor.test(~ a[,gene] + a[,'Generation_Length'], method = 'spearman', continuity = F, conf.level = 0.95)
  rho <- c(rho, rcor$estimate)
  S <- c(S, rcor$statistic)
  p_val <- c(p_val, rcor$p.value)
}


results_nonclients <- data.frame(genes, S, rho, p_val, number_of_species)


results_clients$client <- T
results_nonclients$client <- F
result <- rbind(results_clients, results_nonclients)




library(ggplot2)

pdf('../../Body/4_Figures/hsp.clients.nonclients.spearman.rank.cor.D.vs.gen.length.pdf')

ggplot(result, aes(x=rho, y = p_val, color = client))+
  geom_point()

ggplot(result, aes(y = rho, color = client))+
  geom_boxplot()

ggplot(na.omit(result[result$p_val <= 0.1,]), aes(y = rho, color = client))+
  geom_boxplot()+
  ggtitle('p_val < 0.1')

ggplot(na.omit(result[result$p_val <= 0.01,]), aes(y = rho, color = client))+
  geom_boxplot()+
  ggtitle('p_val < 0.01')

dev.off()

wilcox.test(result[(result$client == T) & (result$p_val <= 0.01), 'rho'], result[(result$client == F) & (result$p_val <= 0.01), 'rho'])


