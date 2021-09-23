rm(list=ls(all=TRUE))
library(tidyverse)
library(tidyr)
library(ggplot2)

clients <-  read.table('../../Body/2_Derived/clients_godon_D.txt', header = T)
nonclients <- read.table('../../Body/2_Derived/nonclients_godon_D.txt', header = T)



cols <- colnames(clients[,-c(1, 335)])
for (c in cols) {
  if (length(clients[!is.na(clients[,c]),c]) <= 70){
    clients[,c] <- NULL
  }
}  


cols <- colnames(nonclients[,-c(1, 429)])
for (c in cols) {
  if (length(nonclients[!is.na(nonclients[,c]),c]) <= 70){
    nonclients[,c] <- NULL
  }
}

clients[,-c(1,295)] <- log10(clients[,-c(1,295)]+0.0001)
nonclients[,-c(1,387)] <- log10(nonclients[,-c(1,387)] + 0.0001)
clients[,295] <- log10(clients[,295])
nonclients[,387] <- log10(nonclients[,387])


#################################
###PICs


library(ape)
library(geiger)
library(caper)
library(stringr)

trees_clietns <- list.files('../../Body/1_Raw/clients/trees/')

rho <- c()
p_val_rho <- c()

#slopes <- c()
#p_val_slope <- c()
#intercept <- c()
#p_val_intercept <- c()
#R_sq <- c()
number_of_species <- c()
genes <- colnames(clients)[-c(1,295)]

#species <- c()
#residuals_clients <- data.frame(clients$Species)



for (gene in genes) {
  
  tree_file <- paste('../../Body/1_Raw/clients/trees/', trees_clietns[grepl(gsub('D_', '', gene), trees_clietns)], sep = '')
  tree <- read.tree(tree_file)
  
  a <- na.omit(clients[,c("Species", gene, 'Generation_Length')])
  a <- a[!duplicated(a[,c(2,3)]),]
  
  row.names(a) = a$Species
  n <- dim(a)[1]
  number_of_species <- c(number_of_species, n)
  
  tree$tip.label <- word(tree$tip.label, 1, 2, sep = '_')
  tree_w = treedata(tree, a, sort=T, warnings=T)$phy
  
  data<-as.data.frame(treedata(tree_w, a, sort=T, warnings=T)$data)
  
  
  data$Species = as.character(data$Species)
  data[,2] = as.numeric(as.character(data[,2]))
  data$Generation_Length = as.numeric(as.character(data$Generation_Length))
  

  
  pic_cor <- cor.test(pic(data[,2], tree_w), pic(data$Generation_Length, tree_w), method = 'spearman')
  r <- pic_cor$estimate
  p_rho <- pic_cor$p.value 
  rho <- c(rho, r)
  p_val_rho <- c(p_val_rho, p_rho)
  
  #a_comp <- comparative.data(tree_w, data, Species, vcv=TRUE)
  
  #model = pgls(Generation_Length ~ a_comp$data[,1], a_comp, lambda="ML")
  #model_sum <- summary(model)
  
  #inter <- model_sum$coefficients[1]
  #intercept <- c(intercept, inter)
  
  #sl <- model_sum$coefficients[2]
  #slopes <- c(slopes, sl)
  
#  p_val_int <- model_sum$coefficients[7]
 # p_val_intercept <- c(p_val_intercept, p_val_int)
  
#  p_val_sl <- model_sum$coefficients[8]
 # p_val_slope <- c(p_val_slope, p_val_sl)
  
  #residuals <- data.frame(rownames(a_comp$data), model_sum$residuals)
  #colnames(residuals) <- c('Species', gsub('D_','', gene))
  #residuals_clients <- merge(residuals_clients, residuals, by.x = 'clients.Species', by.y = 'Species', all = T)
  
  #R_sq <- c(R_sq, model_sum$r.squared)
  
}

#lm_clients <- data.frame(genes, slopes, intercept, p_val_slope, p_val_intercept, number_of_species, R_sq, rho, p_val_rho)

#write.table(lm_clients, '../../Body/3_Results/hsp.clients.linear.model.D.vs.gen.length.txt')


pics_spearman <- data.frame(genes, rho, p_val_rho, number_of_species)
pics_spearman$client <- T



##########################
##############nonclients

trees_nonclietns <- list.files('../../Body/1_Raw/nonclients/trees/')

rho <- c()
p_val_rho <- c()


#slopes <- c()
#p_val_slope <- c()
#intercept <- c()
#p_val_intercept <- c()
#R_sq <- c()
number_of_species <- c()
genes <- colnames(nonclients)[-c(1,387)]

#species <- c()
#residuals_clients <- data.frame(clients$Species)



for (gene in genes) {
  
  tree_file <- paste('../../Body/1_Raw/nonclients/trees/', trees_nonclietns[grepl(gsub('D_', '', gene), trees_nonclietns)], sep = '')
  tree <- read.tree(tree_file)
  
  a <- na.omit(nonclients[,c("Species", gene, 'Generation_Length')])
  a <- a[!duplicated(a[,c(2,3)]),]
  
  row.names(a) = a$Species
  n <- dim(a)[1]
  number_of_species <- c(number_of_species, n)
  
  tree$tip.label <- word(tree$tip.label, 1, 2, sep = '_')
  tree_w = treedata(tree, a, sort=T, warnings=T)$phy
  
  data<-as.data.frame(treedata(tree_w, a, sort=T, warnings=T)$data)
  
  
  data$Species = as.character(data$Species)
  data[,2] = as.numeric(as.character(data[,2]))
  data$Generation_Length = as.numeric(as.character(data$Generation_Length))
  
  pic_cor <- cor.test(pic(data[,2], tree_w), pic(data$Generation_Length, tree_w), method = 'spearman')
  r <- pic_cor$estimate
  p_rho <- pic_cor$p.value 
  rho <- c(rho, r)
  p_val_rho <- c(p_val_rho, p_rho)
  
  #a_comp <- comparative.data(tree_w, data, Species, vcv=TRUE)
  
  #model = pgls(Generation_Length ~ a_comp$data[,1], a_comp, lambda="ML")
  #model_sum <- summary(model)
  
  #inter <- model_sum$coefficients[1]
  #intercept <- c(intercept, inter)
  
  #sl <- model_sum$coefficients[2]
  #slopes <- c(slopes, sl)
  
  #  p_val_int <- model_sum$coefficients[7]
  # p_val_intercept <- c(p_val_intercept, p_val_int)
  
  #  p_val_sl <- model_sum$coefficients[8]
  # p_val_slope <- c(p_val_slope, p_val_sl)
  
  #residuals <- data.frame(rownames(a_comp$data), model_sum$residuals)
  #colnames(residuals) <- c('Species', gsub('D_','', gene))
  #residuals_clients <- merge(residuals_clients, residuals, by.x = 'clients.Species', by.y = 'Species', all = T)
  
  #R_sq <- c(R_sq, model_sum$r.squared)
  
}

nonclients_pics <- data.frame(genes, rho, p_val_rho, number_of_species)
nonclients_pics$client <- F
pics_spearman <- rbind(pics_spearman, nonclients_pics)

ggplot(pics_spearman, aes(x = client, y = rho, fill = client))+
  geom_violin(trim = F, width = 0.6)+
  geom_boxplot(width=0.1)+
  theme_bw()+
  theme(legend.position="none", axis.text=element_text(size=15), axis.title=element_text(size=14,face="bold"))+
  xlab('')+scale_x_discrete(breaks=c(FALSE, TRUE),
                            labels=c('nonclients', 'clients'))


ggplot(pics_spearman[pics_spearman$p_val_rho <= 0.05 & pics_spearman$rho >= 0,], aes(x = client, y = rho, fill = client))+
  geom_violin(trim = F, width = 0.6)+
  geom_boxplot(width = 0.1)+
  theme_bw()+
  theme(legend.position="none", axis.text=element_text(size=15), axis.title=element_text(size=14,face="bold"))+
  xlab('')+scale_x_discrete(breaks=c(FALSE, TRUE),
                            labels=c('nonclients', 'clients'))+
  ggtitle('p_val <= 0.05')


wilcox.test(pics_spearman[pics_spearman$client & pics_spearman$p_val_rho <= 0.05, 'rho'], 
            pics_spearman[!pics_spearman$client & pics_spearman$p_val_rho <= 0.05, 'rho'])

#W = 109, p-value = 0.846








################################################
###All genes together => multiple linear model
###################################################


clients <- clients %>% 
  pivot_longer(colnames(clients[,-c(1,295)]), names_to = "gene", values_to = "D")

nonclients <- nonclients %>% 
  pivot_longer(colnames(nonclients[,-c(1,387)]), names_to = "gene", values_to = "D")

clients$client <- T
nonclients$client <- F
clients <- rbind(clients, nonclients)
write.table(clients, '../../Body/2_Derived/clients.nonclients.D.genlength.all.genes.txt')

ggplot(clients[clients$gene %in% c('D_ENSG00000001084', 'D_ENSG00000004455', 'D_ENSG00000007168', 'D_ENSG00000007312', 'D_ENSG00000009765', 'D_ENSG00000001561', 'D_ENSG00000002330', 'D_ENSG00000002919', 'D_ENSG00000003436', 'D_ENSG00000004534'),], aes(x = Generation_Length, y = D, color = client))+
  geom_point()+
  ylab('log10(D)') + xlab('log10(generation_length)')+
  geom_smooth(method = 'lm')



ggplot(clients, aes(x = Generation_Length, y = D, color = client))+
  geom_point()+
  theme_bw()+
  ylab('log10(D)') + xlab('log10(generation_length)')+
  geom_smooth(method = "glm", se = FALSE)



##################################
model <- glm(D ~ Generation_Length + client, data = clients)
summary(model)
#Deviance Residuals: 
#Min       1Q   Median       3Q      Max  
#-0.7482  -0.6317  -0.5918  -0.5421   5.7112  

#Coefficients:
#                   Estimate Std. Error t value Pr(>|t|)    
#(Intercept)       -2.82713    0.05622 -50.287   <2e-16 ***
#  Generation_Length -0.16745    0.01634 -10.249   <2e-16 ***
#  clientTRUE         0.02787    0.01095   2.546   0.0109 *     

###########



model <-  glm(D ~ Generation_Length + client + as.factor(gene), data = clients)
print(summary(model))

#Call:
#glm(formula = D ~ Generation_Length + client + as.factor(gene), 
#    data = clients)

#Deviance Residuals: 
#  Min       1Q   Median       3Q      Max  
#-1.9062  -0.7076  -0.4916  -0.2570   5.7512  

#Coefficients: (1 not defined because of singularities)
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                      -2.714112   0.163767 -16.573  < 2e-16 ***
#  Generation_Length                -0.167675   0.016105 -10.411  < 2e-16 ***
#  clientTRUE                       -0.278489   0.204811  -1.360 0.173918    
#as.factor(gene)D_ENSG00000001561 -0.113880   0.204399  -0.557 0.577430    
#as.factor(gene)D_ENSG00000002330 -0.109256   0.206529  -0.529 0.596802    
#as.factor(gene)D_ENSG00000002919 -0.367547   0.205656  -1.787 0.073911 .  
#as.factor(gene)D_ENSG00000003436  0.004943   0.208370   0.024 0.981074    
#as.factor(gene)D_ENSG00000004455 -0.020890   0.198994  -0.105 0.916395    
#as.factor(gene)D_ENSG00000004534  0.071052   0.205230   0.346 0.729190    
#as.factor(gene)D_ENSG00000004838  0.261979   0.206089   1.271 0.203665    
#as.factor(gene)D_ENSG00000004961 -0.420243   0.206529  -2.035 0.041876 *  
#  as.factor(gene)D_ENSG00000006042 -0.579872   0.204812  -2.831 0.004638 ** 
#  as.factor(gene)D_ENSG00000006530  0.076983   0.205656   0.374 0.708161    
#as.factor(gene)D_ENSG00000006712 -0.325816   0.204811  -1.591 0.111656    
#as.factor(gene)D_ENSG00000007062  0.376908   0.205656   1.833 0.066849 .  
#as.factor(gene)D_ENSG00000007168 -0.211700   0.191234  -1.107 0.268289    
#as.factor(gene)D_ENSG00000007312  0.105237   0.192156   0.548 0.583925    
#as.factor(gene)D_ENSG00000009709 -0.420598   0.208371  -2.019 0.043543 *  
#  as.factor(gene)D_ENSG00000009765  0.378817   0.189908   1.995 0.046075 *  
#  as.factor(gene)D_ENSG00000010278 -0.197929   0.205230  -0.964 0.334835    
#as.factor(gene)D_ENSG00000010671 -0.178327   0.201451  -0.885 0.376046    
#as.factor(gene)D_ENSG00000010704  0.447067   0.195640   2.285 0.022307 *  
#  as.factor(gene)D_ENSG00000011198  0.020164   0.190343   0.106 0.915634    
#as.factor(gene)D_ENSG00000011376  0.664996   0.190785   3.486 0.000491 ***
#  as.factor(gene)D_ENSG00000012061  0.137588   0.191234   0.719 0.471850    
#as.factor(gene)D_ENSG00000012124  0.866289   0.192628   4.497 6.90e-06 ***
#  as.factor(gene)D_ENSG00000012174 -0.428781   0.206089  -2.081 0.037478 *  
#  as.factor(gene)D_ENSG00000014123 -0.247658   0.205656  -1.204 0.228504    
#as.factor(gene)D_ENSG00000017427 -0.249612   0.204400  -1.221 0.222015    
#as.factor(gene)D_ENSG00000018236 -0.227308   0.204400  -1.112 0.266109    
#as.factor(gene)D_ENSG00000025293 -0.279958   0.204811  -1.367 0.171659    
#as.factor(gene)D_ENSG00000025434 -0.066120   0.189908  -0.348 0.727716    
#as.factor(gene)D_ENSG00000028137  0.762167   0.193109   3.947 7.93e-05 ***
#  as.factor(gene)D_ENSG00000028277 -0.202377   0.206529  -0.980 0.327142    
#as.factor(gene)D_ENSG00000029364 -0.474808   0.204400  -2.323 0.020185 *  
#  as.factor(gene)D_ENSG00000030110 -0.136308   0.208852  -0.653 0.513982    
#as.factor(gene)D_ENSG00000037749 -0.394826   0.205230  -1.924 0.054381 .  
#as.factor(gene)D_ENSG00000038945  0.534909   0.191234   2.797 0.005157 ** 
#  as.factor(gene)D_ENSG00000044115 -0.421878   0.206530  -2.043 0.041086 *  
#  as.factor(gene)D_ENSG00000047457  0.384978   0.206529   1.864 0.062321 .  
#as.factor(gene)D_ENSG00000047579  0.527647   0.192156   2.746 0.006035 ** 








