rm(list = ls(all = TRUE))
library(tidyverse)
library(readxl)
library(phytools)

clients <- read.csv('../../Body/2_Derived/clients/dN.dS.paml.output.clients.nonclients.csv', header = T)
controls <- read.table('../../Body/2_Derived/clients.with.control.nonclients.txt', header = T)
clients <- clients[clients$dN.dS <= 20,]

gene_len <- read_excel("../../../Molecular_evolution/Body/1_Raw/Generation Lenght for Mammals.xlsx")
gene_len$Scientific_name <- gsub(' ', '_', gene_len$Scientific_name)
gene_len <- gene_len[,c(2,3,5,6,8,10,14)]

clients <- merge(clients, gene_len, by.x = 'Species', by.y = 'Scientific_name', all.x = T)
clients$gene <- gsub('_[A-Z|0-9]*', '', clients$gene)

clients$role <- ifelse(clients$gene %in% controls$control_gene, 'nonclient', 'client')

clients$hsp_group <- NA
for (gene in unique(clients$gene)){
  if (gene %in% controls$client){
    if (controls[(controls$client == gene),'hsp90_client'] & controls[(controls$client == gene),'hsc70_client']){
      clients[clients$gene == gene,'hsp_group'] <- 'client_both'}
    if (controls[(controls$client == gene),'hsp90_client'] & !controls[(controls$client == gene),'hsc70_client']){
      clients[clients$gene == gene,'hsp_group'] <- 'client_hsp90'}
    if (!controls[(controls$client == gene),'hsp90_client'] & controls[(controls$client == gene),'hsc70_client']){
      clients[clients$gene == gene,'hsp_group'] <- 'client_hsc70'}
  }
  if (gene %in% controls$control_gene){
    if (controls[(controls$control_gene == gene),'hsp90_client'] & controls[(controls$control_gene == gene),'hsc70_client']){
      clients[clients$gene == gene,'hsp_group'] <- 'control_both'}
    if (controls[(controls$control_gene == gene),'hsp90_client'] & !controls[(controls$control_gene == gene),'hsc70_client']){
      clients[clients$gene == gene,'hsp_group'] <- 'control_hsp90'}
    if (!controls[(controls$control_gene == gene),'hsp90_client'] & controls[(controls$control_gene == gene),'hsc70_client']){
      clients[clients$gene == gene,'hsp_group'] <- 'control_hsc70'}
  }
}




slopes <- c()
p_val_slope <- c()
intercept <- c()
p_val_intercept <- c()
R_sq <- c()
R_sq_adj <- c()
number_of_species <- c()
genes <- unique(clients$gene)
residual_std_err <- c()
residuals_clients <- data.frame(Species = unique(clients$Species))


for (gene in genes) {
  
  a <- na.omit(clients[clients$gene == gene,c(1,6,12)])
  
  species <- a[,'Species']
  number_of_species <- c(number_of_species, length(species))
  
  lg <- lm(a[,'dN.dS']~log10(a[,'GenerationLength_d']), data = a)
  sum <- summary(lg)
  b <- sum$coefficients
  
  slopes <- c(slopes, b[2,1])
  p_val_slope <- c(p_val_slope, b[2,4])
  intercept <- c(intercept, b[1,1])
  p_val_intercept <- c(p_val_intercept, b[1,4])
  R_sq <- c(R_sq, sum$r.squared)
  R_sq_adj <- c(R_sq_adj, sum$adj.r.squared)
  residual_std_err <- c(residual_std_err, sum$sigma)
  res <- data.frame(species, resid(lg))
  colnames(res) <- c('Species', gene)
  residuals_clients <- merge(residuals_clients, res, by = 'Species', all = T)
}


#table with linear regression parameters
lm_clients <- data.frame(genes, slopes, intercept, p_val_slope, p_val_intercept, number_of_species, R_sq, R_sq_adj, residual_std_err)
lm_clients <- merge(lm_clients, clients[,c(3,13,14)][!duplicated(clients[,c(3,13)]),], by.x = 'genes', by.y = 'gene', all.x = T)


lm_clients <- lm_clients[!is.na(lm_clients$hsp_group),]

library(ggplot2)

ggplot(lm_clients[lm_clients$p_val_slope <= 0.05,], aes(x = intercept, y = slopes, color = role, size = p_val_slope), alpha = 0.6)+
  geom_point()+geom_smooth(method = 'lm', se = F)+
  theme_bw()

ggplot(lm_clients[,], aes(x = intercept, y = slopes, color = hsp_group, size = p_val_slope))+
  geom_point(alpha = 0.5)+geom_smooth(method = 'lm', se = F, alpha = 0.5)+
  theme_bw()

ggplot(lm_clients[lm_clients$p_val_slope <= 0.01,], aes(x = log10(p_val_slope), y = slopes, color = hsp_group), alpha = 0.6)+
  geom_point()+geom_smooth(method = 'lm', se = F)+
  theme_bw()


ggplot(lm_clients[lm_clients$hsp_group %in% c('client_both', 'client_hsp90', 'control_both', 'control_hsp90'),], aes(x = intercept, y = slopes, color = role, size = p_val_slope))+
  geom_point(alpha = 0.5)+geom_smooth(method = 'lm', se = F, alpha = 0.5)+
  theme_bw()


ggplot(lm_clients[lm_clients$hsp_group %in% c('client_both', 'client_hsp90', 'control_both', 'control_hsp90') & lm_clients$p_val_slope <= 0.05,], aes(x = intercept, y = slopes, color = role, size = p_val_slope))+
  geom_point(alpha = 0.5)+geom_smooth(method = 'lm', se = F, alpha = 0.5)+
  theme_bw()


ggplot(lm_clients[lm_clients$hsp_group %in% c('client_both', 'client_hsp90', 'control_both', 'control_hsp90') & 
                    (lm_clients$p_val_slope <= 0.1) & (lm_clients$slopes < 2),], 
       aes(x = intercept, y = slopes, color = role, size = p_val_slope))+
  geom_point(alpha = 0.5)+geom_smooth(method = 'lm', se = F, alpha = 0.5)+
  theme_bw()


ggplot(lm_clients[lm_clients$hsp_group %in% c('client_both', 'client_hsc70', 'control_both', 'control_hsc70') & 
                    (lm_clients$p_val_slope <= 1) & (lm_clients$slopes < 1),], 
       aes(x = intercept, y = slopes, color = role, size = p_val_slope))+
  geom_point(alpha = 0.5)+geom_smooth(method = 'lm', se = F, alpha = 0.5)+
  theme_bw()


wilcox.test(lm_clients[lm_clients$role == 'client',]$slopes, lm_clients[lm_clients$role == 'nonclient',]$slopes)


ggplot(lm_clients[,], aes(x = role, y = slopes, fill = role))+
         geom_boxplot()+
  theme_bw()+
  #scale_fill_discrete(name = '√руппа', labels = c('клиенты', 'неклиенты'))+
  #scale_x_discrete(name = '', labels = c('клиенты', 'неклиенты'))+
  #ylab('наклон регрессионной пр€мой')+
  theme(axis.text.x = element_text(size = 15), axis.title.y = element_text(size = 15), 
        legend.title = element_text(size = 16, face = 'bold'), legend.text = element_text(size = 14))







lm_clients$hsp_group <- factor(lm_clients$hsp_group , levels=c("client_hsp90", "control_hsp90", "client_hsc70", "control_hsc70",
                                                               'client_both', 'control_both'))


ggplot(lm_clients[lm_clients$hsp_group %in% c('client_both', 'client_hsp90', 'control_both', 'control_hsp90') &
                    (lm_clients$p_val_slope <= 0.1),], 
       aes(x = role,y = slopes, fill = role))+
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+ylim(c(-0.4,1))+
  #scale_x_discrete(name = "", labels = c('клиенты HSP90', 'неклиенты'))+
  #ylab('наклон регрессионной пр€мой')+
  theme(axis.text.x = element_text(size = 15), axis.title.y = element_text(size = 15), 
        legend.position = 'None', axis.text.y = element_text(size = 13))

wilcox.test(lm_clients[lm_clients$hsp_group %in% c('client_both', 'client_hsp90'), 'slopes'],
            lm_clients[lm_clients$hsp_group %in% c('control_both', 'control_hsp90'), 'slopes'])
#W = 2483, p-value = 0.458

ggplot(lm_clients[lm_clients$hsp_group %in% c('client_both', 'client_hsc70', 'control_both', 'control_hsc70') &
                    (lm_clients$p_val_slope <= 1),], 
       aes(x = role,y = slopes, fill = role))+
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+ylim(c(-0.4,1))+
  scale_x_discrete(name = "", labels = c('клиенты HSC70', 'неклиенты'))+
  ylab('наклон регрессионной пр€мой')+
  theme(axis.text.x = element_text(size = 15), axis.title.y = element_text(size = 15), 
        legend.position = 'None', axis.text.y = element_text(size = 13))


wilcox.test(lm_clients[lm_clients$hsp_group %in% c('client_both', 'client_hsc70'), 'slopes'],
            lm_clients[lm_clients$hsp_group %in% c('control_both', 'control_hsc70'), 'slopes'])
#W = 1213, p-value = 0.7754


ggplot(lm_clients, aes(x = hsp_group,y = slopes, fill = role))+
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+
  ylim(c(-0.2, 1))+
  scale_x_discrete(name = "", labels = c('HSP90', 'HSP90', 'HSC70', 'HSC70',
                                                  'HSP90 и HSC70', 'HSP90 и HSC70'))+
  scale_fill_discrete(name = '√руппа', labels = c('клиент', 'неклиент'))+
  ylab('наклон регрессионной пр€мой')+
  xlab('aaa')+
    theme(axis.text.x = element_text(size = 12, angle = 15, hjust = 0.5 ), axis.title.y = element_text(size = 15), 
          legend.title = element_text(size = 16, face = 'bold'), legend.text = element_text(size = 14))



ggplot(lm_clients, aes(x = hsp_group,y = intercept, fill = role))+
  geom_boxplot(outlier.shape = 1)+
  theme_bw()+
  #ylim(c(-0.2, 1))+
  scale_x_discrete(name = "", labels = c('HSP90', 'HSP90', 'HSC70', 'HSC70',
                                         'HSP90 и HSC70', 'HSP90 и HSC70'))+
  scale_fill_discrete(name = '√руппа', labels = c('клиент', 'неклиент'))+
  ylab('интерсепта')+
  theme(axis.text.x = element_text(size = 12, angle = 15, hjust = 0.5 ), axis.title.y = element_text(size = 15), 
        legend.title = element_text(size = 16, face = 'bold'), legend.text = element_text(size = 14))

table(lm_clients[lm_clients$hsp_group %in% c('client_hsp90', 'client_both', 'control_hsp90', 'control_both') & !(is.na(lm_clients$slopes)),'role'])
#client nonclient 
#75        85 

p1 <- ggplot(lm_clients[lm_clients$hsp_group %in% c('client_hsp90', 'client_both', 'control_hsp90', 'control_both'),],
       aes(x = role, y = slopes, fill = role))+
  geom_boxplot(outlier.shape  = NA)+
  ylim(0,0.75)+
  theme_bw()+
  theme(axis.title = element_text(size = 27),
        axis.text = element_text(size = 22), legend.position = "none")+
  scale_x_discrete(name = "", labels = c('HSP90 clients (N = 75)', 'nonclients (N = 85)'))+
  ylab('slope')
ggsave('../../Body/4_Figures/clients.nonclients.paml.dN.dS.vs.genlen.hsp90.pdf', p1)

wilcox.test(lm_clients[lm_clients$hsp_group %in% c('client_hsp90', 'client_both'), 'slopes'], 
            lm_clients[lm_clients$hsp_group %in% c('control_hsp90', 'control_both'), 'slopes'], alternative = 'greater')

#W = 3538, p-value = 0.1157



###################
table(lm_clients[lm_clients$hsp_group %in% c('client_hsc70', 'client_both', 'control_hsc70', 'control_both') & !(is.na(lm_clients$slopes)),'role'])
#client nonclient 
#52        64

p2 <- ggplot(lm_clients[lm_clients$hsp_group %in% c('client_hsc70', 'client_both', 'control_hsc70', 'control_both') & lm_clients$slopes >=0.05,],
       aes(x = role, y = slopes, fill = role))+
  geom_boxplot(outlier.shape  = NA)+
  ylim(0,0.75)+
  theme_bw()+
  theme(axis.title = element_text(size = 27),
        axis.text = element_text(size = 22), legend.position = "none")+
  scale_x_discrete(name = "", labels = c('HSC70 clients (N = 52)', 'nonclients (N = 64)'))+
  ylab('slope')
ggsave('../../Body/4_Figures/clients.nonclients.paml.dN.dS.vs.genlen.hsc70.pdf', p2)

wilcox.test(lm_clients[lm_clients$hsp_group %in% c('client_hsc70', 'client_both'), 'slopes'], 
            lm_clients[lm_clients$hsp_group %in% c('control_hsc70', 'control_both'), 'slopes'], alternative = 'greater')

#W = 1748, p-value = 0.3215



######

wilcox.test(lm_clients[lm_clients$hsp_group == 'client_hsp90',]$slopes, lm_clients[lm_clients$hsp_group == 'control_hsp90',]$slopes)
#W = 127, p-value = 0.7626

wilcox.test(lm_clients[lm_clients$hsp_group == 'client_hsc70',]$slopes, lm_clients[lm_clients$hsp_group == 'control_hsc70',]$slopes)
#W = 27, p-value = 0.5908

wilcox.test(lm_clients[lm_clients$hsp_group == 'client_both',]$slopes, lm_clients[lm_clients$hsp_group == 'control_both',]$slopes)
#W = 207, p-value = 0.1556

wilcox.test(lm_clients[lm_clients$hsp_group == 'client_hsp90' | lm_clients$hsp_group == 'client_both',]$slopes, lm_clients[lm_clients$hsp_group == 'control_hsp90' | lm_clients$hsp_group == 'control_both',]$slopes)

wilcox.test(lm_clients[lm_clients$hsp_group == 'client_hsc70' | lm_clients$hsp_group == 'client_both',]$slopes, lm_clients[lm_clients$hsp_group == 'control_hsc70' | lm_clients$hsp_group == 'control_both',]$slopes)


length(lm_clients[lm_clients$p_val_slope<0.001 & lm_clients$role == 'client',]$slopes)/length(lm_clients[lm_clients$role == 'client',]$slopes)

length(lm_clients[lm_clients$p_val_slope<0.001 & lm_clients$role == 'nonclient',]$slopes)/length(lm_clients[lm_clients$role == 'nonclient',]$slopes)


wilcox.test(lm_clients[lm_clients$role == 'client',]$slopes, lm_clients[lm_clients$role == 'nonclient',]$slopes)

lg <- lm(clients$dN.dS ~ log10(clients$GenerationLength_d) * clients$role)
summary(lg)
#Residuals:
#  Min      1Q  Median      3Q     Max 
#-0.4341 -0.2496 -0.1292  0.0038 18.4794 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                                             -0.59281    0.10511  -5.640 1.77e-08 ***
#  log10(clients$GenerationLength_d)                        0.25351    0.03076   8.240  < 2e-16 ***
#  clients$rolenonclient                                   -0.09887    0.14454  -0.684    0.494    
#log10(clients$GenerationLength_d):clients$rolenonclient  0.02994    0.04230   0.708    0.479    
#---
#  Signif. codes:  0 С***Т 0.001 С**Т 0.01 С*Т 0.05 С.Т 0.1 С Т 1

#Residual standard error: 0.6048 on 7122 degrees of freedom
#(1623 observations deleted due to missingness)
#Multiple R-squared:  0.02241,	Adjusted R-squared:  0.022 
#F-statistic: 54.41 on 3 and 7122 DF,  p-value: < 2.2e-16






lg2 <- lm(clients$dN.dS ~ log10(clients$GenerationLength_d) * clients$hsp_group)
summary(lg2)

clients$hsp_group2 <- ifelse(clients$role == 'nonclient', 'control', clients$hsp_group)
lg3 <- lm(clients$dN.dS ~ log10(clients$GenerationLength_d) * clients$hsp_group2)
summary(lg3)

clients$hsp90_client <- ifelse(clients$hsp_group == 'client_hsp90' | clients$hsp_group == 'client_both', 1, 0)
clients$hsc70_client <- ifelse(clients$hsp_group == 'client_hsc70' | clients$hsp_group == 'client_both', 1, 0)
lg4 <- lm(clients$dN.dS ~ log10(clients$GenerationLength_d) + clients$hsp90_client * clients$hsc70_client)
summary(lg4)

########### t-test ###################

lm_clients$log.slope <- log10(lm_clients$slopes + 0.5)

shapiro.test(lm_clients$log.slope)
plot(density(lm_clients$log.slope))

qqnorm(lm_clients$log.slope)
qqline(lm_clients$log.slope)


t.test(lm_clients[lm_clients$role == 'client',]$log.slope, lm_clients[lm_clients$role == 'nonclient',]$log.slope)
#t = -0.2335, df = 186.67, p-value = 0.8156

t.test(lm_clients[lm_clients$hsp_group == 'client_hsp90',]$log.slope, lm_clients[lm_clients$hsp_group == 'control_hsp90',]$log.slope)
#t = 0.97277, df = 80.007, p-value = 0.3336

t.test(lm_clients[lm_clients$hsp_group == 'client_hsc70',]$log.slope, lm_clients[lm_clients$hsp_group == 'control_hsc70',]$log.slope)
#t = -0.46223, df = 30.96, p-value = 0.6471

t.test(lm_clients[lm_clients$hsp_group == 'client_both',]$log.slope, lm_clients[lm_clients$hsp_group == 'control_both',]$log.slope)
#t = -0.80187, df = 67.361, p-value = 0.4254

t.test(lm_clients[lm_clients$hsp_group == 'client_hsp90' | lm_clients$hsp_group == 'client_both',]$log.slope, lm_clients[lm_clients$hsp_group == 'control_hsp90' | lm_clients$hsp_group == 'control_both',]$log.slope)
#t = -0.023671, df = 151.84, p-value = 0.9811

t.test(lm_clients[lm_clients$hsp_group == 'client_hsc70' | lm_clients$hsp_group == 'client_both',]$log.slope, lm_clients[lm_clients$hsp_group == 'control_hsc70' | lm_clients$hsp_group == 'control_both',]$log.slope)
#t = -0.9311, df = 101.19, p-value = 0.354






############### paired wilcox ###########################

clients_lm_paired <- merge(controls, lm_clients[,c(1,2)], by.x = 'client', by.y = 'genes', all.x = T)
clients_lm_paired <- merge(clients_lm_paired, lm_clients[,c(1,2)], by.x = 'control_gene', by.y = 'genes', all.x = T)
colnames(clients_lm_paired) <- c( "control_gene", "client"   ,    "hsp90_client", "hsc70_client", "slope_client",     "slope_nonclient")
clients_lm_paired <- na.omit(clients_lm_paired)

wilcox.test(clients_lm_paired$slope_client, clients_lm_paired$slope_nonclient, paired = T)
#V = 1638, p-value = 0.3663

wilcox.test(clients_lm_paired[clients_lm_paired$hsp90_client,]$slope_client, clients_lm_paired[clients_lm_paired$hsp90_client,]$slope_nonclient, paired = T)
#V = 1167, p-value = 0.2779

wilcox.test(clients_lm_paired[clients_lm_paired$hsc70_client,]$slope_client, clients_lm_paired[clients_lm_paired$hsc70_client,]$slope_nonclient, paired = T)
#V = 427, p-value = 0.6143

wilcox.test(clients_lm_paired[clients_lm_paired$hsp90_client & !clients_lm_paired$hsc70_client,]$slope_client, clients_lm_paired[clients_lm_paired$hsp90_client & !clients_lm_paired$hsc70_client,]$slope_nonclient, paired = T)
#V = 389, p-value = 0.5809


wilcox.test(clients_lm_paired[!clients_lm_paired$hsp90_client & clients_lm_paired$hsc70_client,]$slope_client, clients_lm_paired[!clients_lm_paired$hsp90_client & clients_lm_paired$hsc70_client,]$slope_nonclient, paired = T)
#V = 45, p-value = 1

wilcox.test(clients_lm_paired[clients_lm_paired$hsp90_client & clients_lm_paired$hsc70_client,]$slope_client, clients_lm_paired[clients_lm_paired$hsp90_client & clients_lm_paired$hsc70_client,]$slope_nonclient, paired = T)
#V = 205, p-value = 0.4678





######## PGLS ##########

library(ape)
library(geiger)
library(caper)


slopes <- c()
p_val_slope <- c()
intercept <- c()
p_val_intercept <- c()
R_sq <- c()
R_sq_adj <- c()
number_of_species <- c()
genes <- unique(clients$gene)
residual_std_err <- c()
residuals_clients <- data.frame(Species = unique(clients$Species))
trees_clients <- list.files('../../Body/2_Derived/clients/tree/')
trees_nonclients <- list.files('../../Body/2_Derived/nonclients/tree/')

genes_pgls <- c()

for (g in genes){

  if (sum(grepl(g, trees_clients)) > 0){
    
    tree <- read.tree(paste('../../Body/2_Derived/clients/tree/',trees_clients[grepl(g, trees_clients)], sep = ''))
  }
  
  else if (sum(grepl(g, trees_nonclients)) > 0) {
    tree <- read.tree(paste('../../Body/2_Derived/nonclients/tree/',trees_nonclients[grepl(g, trees_nonclients)], sep = ''))
  }
  
  else {next}
  
  genes_pgls <- c(genes_pgls, g)
  
  gene <- clients[clients$gene == g,]
  row.names(gene) = gene$Species
  gene <- na.omit(gene[,c(1,6,12)])
  
  tree_w = treedata(tree, gene, sort=T, warnings=T)$phy
  data<-as.data.frame(treedata(tree_w, gene, sort=T, warnings=T)$data)
  
  
  data$Species = as.character(data$Species)
  
  data$dN.dS = as.numeric(as.character(data$dN.dS))
  data$GenerationLength_d = as.numeric(as.character(data$GenerationLength_d))
  
  gene_comp = comparative.data(tree_w, data, Species, vcv=TRUE)

  model = pgls(dN.dS ~ log10(GenerationLength_d), gene_comp, lambda="ML")
  sum_pgls <- summary(model)
  
  slopes <- c(slopes, sum_pgls$coefficients[2])
  intercept <- c(intercept, sum_pgls$coefficients[1])
  p_val_slope <- c(p_val_slope, sum_pgls$coefficients[8])
  p_val_intercept <- c(p_val_intercept, sum_pgls$coefficients[7])
  R_sq <- c(R_sq, sum_pgls$r.squared)
  R_sq_adj <- c(R_sq_adj, sum_pgls$adj.r.squared[1])
  
  number_of_species <- c(number_of_species, length(data$Species))
  residual_std_err <- c(residual_std_err, sum_pgls$sigma)
  
  res <- data.frame(resid(model))
  res$Species = rownames(res)
  colnames(res) <- c(g, 'Species')
  residuals_clients <- merge(residuals_clients, res, by = 'Species', all = T)  
  
}



pgls_result <- data.frame(data.frame(genes_pgls, slopes, intercept, p_val_slope, p_val_intercept, number_of_species, R_sq, R_sq_adj, residual_std_err))
pgls_result <- merge(pgls_result, clients[,c(3,13,14)][!duplicated(clients[,c(3,13)]),], by.x = 'genes_pgls', by.y = 'gene', all.x = T)
pgls_result <- pgls_result[!is.na(pgls_result$hsp_group),]


ggplot(pgls_result, aes(x = intercept, y = slopes, size = p_val_slope, color = hsp_group))+
  #geom_point()+
  geom_smooth(method = 'lm')




wilcox.test(pgls_result[pgls_result$role == 'client',]$slopes, pgls_result[pgls_result$role == 'nonclient',]$slopes, alternative = 'greater')
#W = 5497, p-value = 0.101

ggplot(pgls_result, aes(x = role, y = slopes, fill = role))+
  geom_boxplot()+
  theme_bw()+
  scale_fill_discrete(name = '√руппа', labels = c('клиенты', 'неклиенты'))+
  scale_x_discrete(name = '', labels = c('клиенты', 'неклиенты'))+
  ylab('наклон регрессионной пр€мой')+
  theme(axis.text.x = element_text(size = 15), axis.title.y = element_text(size = 15), 
        legend.position = 'None', legend.text = element_text(size = 14))









pgls_result$hsp_group <- factor(pgls_result$hsp_group , levels=c("client_hsp90", "control_hsp90", "client_hsc70", "control_hsc70",
                                                               'client_both', 'control_both'))
ggplot(pgls_result, aes(x = hsp_group,y = slopes, fill = role))+
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+
  ylim(c(-0.2, 1))+
  scale_x_discrete(name = "", labels = c('HSP90', 'HSP90', 'HSC70', 'HSC70',
                                         'HSP90 и HSC70', 'HSP90 и HSC70'))+
  scale_fill_discrete(name = '√руппа', labels = c('клиенты', 'неклиенты'))+
  ylab('наклон линейной регрессии')+
  xlab('aaa')+
  theme(axis.text.x = element_text(size = 12, angle = 15, hjust = 0.5 ), axis.title.y = element_text(size = 15), 
        legend.title = element_text(size = 16, face = 'bold'), legend.text = element_text(size = 14))






wilcox.test(pgls_result[pgls_result$hsp_group == 'client_hsp90',]$slopes, pgls_result[pgls_result$hsp_group == 'control_hsp90',]$slopes, alternative = 'greater')
#W = 1026, p-value = 0.09932

wilcox.test(pgls_result[pgls_result$hsp_group == 'client_hsc70',]$slopes, pgls_result[pgls_result$hsp_group == 'control_hsc70',]$slopes, alternative = 'greater')
#W = 195, p-value = 0.5375

wilcox.test(pgls_result[pgls_result$hsp_group == 'client_both',]$slopes, pgls_result[pgls_result$hsp_group == 'control_both',]$slopes, alternative = 'greater')
#W = 767, p-value = 0.2927

wilcox.test(pgls_result[pgls_result$hsp_group == 'client_hsp90' | pgls_result$hsp_group == 'client_both',]$slopes, 
            pgls_result[pgls_result$hsp_group == 'control_hsp90' | pgls_result$hsp_group == 'control_both',]$slopes, alternative = 'greater')
#W = 3607, p-value = 0.07597


wilcox.test(pgls_result[pgls_result$hsp_group == 'client_hsc70' | pgls_result$hsp_group == 'client_both',]$slopes, 
            pgls_result[pgls_result$hsp_group == 'control_hsc70' | pgls_result$hsp_group == 'control_both',]$slopes, alternative = 'greater')
#W = 1753, p-value = 0.3116



length(pgls_result[pgls_result$p_val_slope<0.01 & pgls_result$role == 'client',]$slopes)/length(pgls_result[pgls_result$role == 'client',]$slopes)

length(pgls_result[pgls_result$p_val_slope<0.01 & pgls_result$role == 'nonclient',]$slopes)/length(pgls_result[pgls_result$role == 'nonclient',]$slopes)


ggplot(pgls_result, aes(x = hsp_group, y = slopes, fill = role))+
  geom_boxplot()+
  theme(axis.text.x = element_text(size = 15), axis.title.y = element_text(size = 15), 
        legend.position = 'None', legend.text = element_text(size = 14))


##paired

pgls_paired <- merge(controls, pgls_result[,c(1,2)], by.x = 'client', by.y = 'genes_pgls', all.x = T)
pgls_paired <- merge(pgls_paired, pgls_result[,c(1,2)], by.x = 'control_gene', by.y = 'genes_pgls', all.x = T)
colnames(pgls_paired) <- c( "control_gene", "client"   ,    "hsp90_client", "hsc70_client", "slope_client",     "slope_nonclient")
pgls_paired <- na.omit(pgls_paired)

wilcox.test(pgls_paired$slope_client, pgls_paired$slope_nonclient, paired = T, alternative = 'greater')
#V = 1940, p-value = 0.09449

wilcox.test(pgls_paired[pgls_paired$hsp90_client,]$slope_client, pgls_paired[pgls_paired$hsp90_client,]$slope_nonclient, paired = T, alternative = 'greater')
#V = 1382, p-value = 0.06491

wilcox.test(pgls_paired[pgls_paired$hsc70_client,]$slope_client, pgls_paired[pgls_paired$hsc70_client,]$slope_nonclient, paired = T, alternative = 'greater')
#V = 532, p-value = 0.2419

wilcox.test(pgls_paired[pgls_paired$hsp90_client & !pgls_paired$hsc70_client,]$slope_client, pgls_paired[pgls_paired$hsp90_client & !pgls_paired$hsc70_client,]$slope_nonclient, paired = T, alternative = 'greater')
#V = 404, p-value = 0.2186


wilcox.test(pgls_paired[!pgls_paired$hsp90_client & pgls_paired$hsc70_client,]$slope_client, pgls_paired[!pgls_paired$hsp90_client & pgls_paired$hsc70_client,]$slope_nonclient, paired = T, alternative = 'greater')
#V = 44, p-value = 0.5537

wilcox.test(pgls_paired[pgls_paired$hsp90_client & pgls_paired$hsc70_client,]$slope_client, pgls_paired[pgls_paired$hsp90_client & pgls_paired$hsc70_client,]$slope_nonclient, paired = T, alternative = 'greater')
#V = 222, p-value = 0.1236


ggplot(pgls_result[pgls_result$hsp_group %in% c('client_hsp90', 'client_both', 'control_hsp90', 'control_both'),],
       aes(x = role, y = slopes, fill = role))+
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+
  ylim(c(-0.2, 1))+
  theme(axis.text.x = element_text(size = 29), axis.title.y = element_text(size = 29), 
        legend.position = 'None', axis.text.y = element_text(size = 23))+
  xlab('')+
  ylab(paste('slope of lm(dN/dS ~', 'log10(Generation_length))', sep = '\n'))+
  scale_x_discrete(name = "", labels = c('HSP90 clients', 'nonclients'))

ggplot(pgls_result[pgls_result$hsp_group %in% c('client_hsc70', 'client_both', 'control_hsc70', 'control_both'),],
       aes(x = role, y = slopes, fill = role))+
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+
  ylim(c(-0.2, 1))

pgls_result$log.slope <- log10(pgls_result$slopes + 0.5)
t.test(pgls_result[pgls_result$hsp_group %in% c('client_hsp90', 'control_hsp90'), 'log.slope'] ~ 
       pgls_result[pgls_result$hsp_group %in% c('client_hsp90', 'control_hsp90'), 'role'], alternative = 'greater')




##################redbook and residuals #######################

residuals_clients <- residuals_clients[!duplicated(residuals_clients),]
residuals_nonclients <- residuals_clients[,colnames(residuals_clients) %in% controls$control_gene]
residuals_nonclients$Species <- residuals_clients$Species
residuals_clients <- residuals_clients[,colnames(residuals_clients) %in% controls$client]
residuals_clients$Species <- residuals_nonclients$Species


redbook <- read.csv('../../Body/1_Raw/IUCN.csv', sep = ';', header = T)
redbook$specieName <- gsub(' ', '_', redbook$specieName)
redbook <- redbook[redbook$class == 'MAMMALIA',]
redbook <- redbook[,c(3,14)]

residuals_clients <- merge(residuals_clients, redbook, by.x  = 'Species', by.y = 'specieName', all.x = T)
residuals_clients$category <- factor(residuals_clients$category , levels=c("LC", 'NT', 'VU', 'EN', 'CR', 'DD'))
residuals_clients <- residuals_clients[!is.na(residuals_clients$category),]

residuals_nonclients <- merge(residuals_nonclients, redbook, by.x  = 'Species', by.y = 'specieName', all.x = T)
residuals_nonclients$category <- factor(residuals_nonclients$category , levels=c("LC", 'NT', 'VU', 'EN', 'CR', 'DD'))
residuals_nonclients <- residuals_nonclients[!is.na(residuals_nonclients$category),]

residuals_clients <- residuals_clients %>%
  pivot_longer(colnames(residuals_clients)[-c(1,84)], names_to = 'gene', values_to = 'residual')

residuals_clients <- merge(residuals_clients, clients[,c(3,14)], by = 'gene', all.x = T)
residuals_clients <- residuals_clients[!duplicated(residuals_clients),]

residuals_nonclients <- residuals_nonclients %>%
  pivot_longer(colnames(residuals_nonclients)[-c(1,87)], names_to = 'gene', values_to = 'residual')

residuals_nonclients <- merge(residuals_nonclients, clients[,c(3,14)], by = 'gene', all.x = T)
residuals_nonclients <- residuals_nonclients[!duplicated(residuals_nonclients),]

ggplot(residuals_clients, aes(x = category, y = residual, fill = hsp_group))+
  geom_boxplot(outlier.shape = NA)+
  ylim(c(-0.7, 0.5))

ggplot(residuals_nonclients, aes(x = category, y = residual, fill = hsp_group))+
  geom_boxplot(outlier.shape = NA)+
  ylim(c(-0.7, 0.5))

residuals_all <- rbind(residuals_clients, residuals_nonclients)

residuals_all$hsp_group <- factor(residuals_all$hsp_group , levels=c("client_hsp90", 'control_hsp90', 'client_hsc70', 'control_hsc70',
                                                                     'client_both', 'control_both'))

ggplot(residuals_all, aes(x = category, y = residual, fill = hsp_group))+
  geom_boxplot(outlier.shape = NA)+
  ylim(c(-0.7, 0.5))

residuals_all <- residuals_all[residuals_all$category != 'DD',]
residuals_all$binar_cat <- ifelse(residuals_all$category == 'LC', 'LC', 'others')

ggplot(residuals_all, aes(x = hsp_group, y = residual, fill = binar_cat))+
  geom_boxplot(outlier.shape = NA)+
  ylim(c(-0.7, 0.5))

ggplot(residuals_all, aes(x = binar_cat, y = residual, fill = hsp_group))+
  geom_boxplot(outlier.shape = NA)+
  ylim(c(-0.8,0.5))



