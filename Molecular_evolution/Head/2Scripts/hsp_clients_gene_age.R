rm(list = ls(all = TRUE))
library('ggplot2')
library(tidyverse)

metrics_table <- read.table('../../Body/2_Derived/gencode.v25.annotation.gtf.Genes.Shet.pLI.FIS.RVIS.GHIS.KnKs.GC.BrainSpecificRanking.Branch.gnomad.oe_lof_upper_bin.p', header = T)
pdf('../../Body/4_Figures/Branch_clients_int_score.pdf')
########## hsp90 ###############
hsp90_int <- read.csv('../../Body/2_Derived/hsp90_clients.Elisa_more_0.2.WT_interaction_score.csv', header = T)

hsp90_int_age <- merge(hsp90_int, metrics_table[,c(1,19)], by.x = 'EnsemlID', by.y = 'EnsemblId', all = T)

ggplot(na.omit(hsp90_int_age[,]), aes(x = as.factor(Branch)))+
  #geom_bar()+
  geom_bar(color = 'blue', alpha = 0)+
  geom_boxplot(aes(y = WT_interaction_score/0.05), color = 'red', outlier.shape = NA, alpha = 0)+
  #ylim(-5.5, 600)+
  scale_y_continuous(name = 'Number of genes', sec.axis = sec_axis(trans =~.*0.05, name = 'Interaction score'))+
  theme_bw()+ggtitle('hsp90')

hsp90_int_age$client <- ifelse(hsp90_int_age$WT_interaction_score > 3, 1, 0)
hsp90_int_age[is.na(hsp90_int_age$client),'client'] <- 0

hsp90_int_age$Branch_zero <- ifelse(hsp90_int_age$Branch == 0, 1, 0)
mosaicplot(table(hsp90_int_age$client, hsp90_int_age$Branch_zero), xlab = 'hsp90 client', ylab = 'Branch is zero')
fisher.test(table(hsp90_int_age$client, hsp90_int_age$Branch_zero))
#p-value = 0.003979
#odds ratio 
#1.765529 

clients_branch <- data.frame(table(hsp90_int_age$Branch, hsp90_int_age$client))
table(hsp90_int_age$client)

clients_branch <- clients_branch %>%
  pivot_wider(names_from = Var2, values_from = Freq)
colnames(clients_branch) <- c('Branch', 'Nonclients', 'Clients')
clients_branch$Num_of_genes <- apply(clients_branch[,c(2,3)], 1, sum)
clients_branch$Fr_of_clients <- clients_branch$Clients/clients_branch$Num_of_genes


ggplot(clients_branch, aes(x = as.factor(Branch)))+
  #geom_bar()+
  geom_col(aes(y = Num_of_genes), fill = 'blue', alpha = 0.5)+
  geom_col(aes(y = Fr_of_clients/0.000001), fill = 'red', alpha = 0.5)+
  scale_y_continuous(name = 'Number of genes', sec.axis = sec_axis(trans =~.*0.000001, name = 'Fraction of clients'))+
  theme_bw()+ggtitle('hsp90')


ggplot(na.omit(hsp90_int_age), aes(x = as.factor(Branch_zero), y = WT_interaction_score))+
  geom_violin()+
  theme_bw()+ggtitle('hsp90')

######### hsc70 #########

hsc70_int <- read.csv('../../Body/2_Derived/hsc70_clients.Elisa_more_0.2.WT_interaction_score.csv', header = T)

hsc70_int_age <- merge(hsc70_int, metrics_table[,c(1,19)], by.x = 'EnsemlID', by.y = 'EnsemblId', all = T)

ggplot(na.omit(hsc70_int_age[,]), aes(x = as.factor(Branch)))+
  #geom_bar()+
  geom_bar(color = 'blue', alpha = 0)+
  geom_boxplot(aes(y = WT_interaction_score/0.01), color = 'red', outlier.shape = NA, alpha = 0)+
  #ylim(-5.5, 600)+
  scale_y_continuous(name = 'Number of genes', sec.axis = sec_axis(trans =~.*0.01, name = 'Interaction score'))+
  theme_bw()+ggtitle('hsc70')

hsc70_int_age$client <- ifelse(hsc70_int_age$WT_interaction_score > 3, 1, 0)
hsc70_int_age[is.na(hsc70_int_age$client),'client'] <- 0

hsc70_int_age$Branch_zero <- ifelse(hsc70_int_age$Branch == 0, 1, 0)
mosaicplot(table(hsc70_int_age$client, hsc70_int_age$Branch_zero), xlab = 'hsc70 client', ylab = 'Branch is zero')

fisher.test(table(hsc70_int_age$client, hsc70_int_age$Branch_zero))
#p-value = 0.06418
#odds ratio 
#1.669555 

hsc70_clients_branch <- data.frame(table(hsc70_int_age$Branch, hsc70_int_age$client))


hsc70_clients_branch <- hsc70_clients_branch %>%
  pivot_wider(names_from = Var2, values_from = Freq)
colnames(hsc70_clients_branch) <- c('Branch', 'Nonclients', 'Clients')
hsc70_clients_branch$Num_of_genes <- apply(hsc70_clients_branch[,c(2,3)], 1, sum)
hsc70_clients_branch$Fr_of_clients <- hsc70_clients_branch$Clients/hsc70_clients_branch$Num_of_genes


ggplot(hsc70_clients_branch, aes(x = as.factor(Branch)))+
  #geom_bar()+
  geom_col(aes(y = Num_of_genes), fill = 'blue', alpha = 0.5)+
  geom_col(aes(y = Fr_of_clients/0.000001), fill = 'red', alpha = 0.5)+
  scale_y_continuous(name = 'Number of genes', sec.axis = sec_axis(trans =~.*0.000001, name = 'Fraction of clients'))+
  theme_bw()+
  ggtitle('hsc70')


ggplot(na.omit(hsc70_int_age), aes(x = as.factor(Branch_zero), y = WT_interaction_score))+
  geom_violin()+
  theme_bw()+ggtitle('hsc70')

dev.off()

