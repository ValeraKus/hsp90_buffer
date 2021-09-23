rm(list = ls(all = TRUE))
library(tidyverse)
library(readxl)
library(phytools)

hsps <- read.csv('../../Body/2_Derived/hsps/dN.dS.paml.output.hsps.csv', header = T)
controls <- read.table('../../Body/2_Derived/hsps.with.control.nonhsps.txt')
hsps <- hsps[hsps$dN.dS <= 20,]
hsps$gene <- gsub('_[A-Z|0-9]*', '', hsps$gene)
hsps$role <- ifelse(hsps$gene %in% controls$hsp, 'hsp', 'nonhsp')
hsps$Genus <- gsub('_[a-z]*', '', hsps$Species)

gene_len <- read_excel("../../../Molecular_evolution/Body/1_Raw/Generation Lenght for Mammals.xlsx")
gene_len$Scientific_name <- gsub(' ', '_', gene_len$Scientific_name)

hsps <- merge(hsps, gene_len[,c(2,3,5,6,8,10,14)], by.x = 'Species', by.y = 'Scientific_name', all.x = T)






hsp90 <- hsps[hsps$gene == 'ENSG00000096384' | hsps$gene ==  controls[controls$hsp == 'ENSG00000096384',]$control_gene,]
hsp90_control <- hsps[hsps$gene ==  controls[controls$hsp == 'ENSG00000096384',]$control_gene,]


summary(lm(hsp90$dN.dS ~ log10(hsp90$GenerationLength_d) * hsp90$gene))
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)                      0.0341252  0.1140386   0.299    0.766
#log10(hsp90$GenerationLength_d) -0.0008238  0.0331299  -0.025    0.980


summary(lm(hsp90_control$dN.dS ~ log10(hsp90_control$GenerationLength_d)))
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)  
#(Intercept)                              -1.1354     0.5690  -1.996   0.0492 *
#  log10(hsp90_control$GenerationLength_d)   0.3948     0.1671   2.363   0.0204 *


hsc70 <- hsps[hsps$gene == 'ENSG00000109971' | hsps$gene == 'ENSG00000009335',]
summary(lm(hsc70[hsc70$gene == 'ENSG00000109971',]$dN.dS ~ log10(hsc70[hsc70$gene == 'ENSG00000109971',]$GenerationLength_d)))
#Coefficients:
#                                                                     Estimate Std. Error t value Pr(>|t|)
#(Intercept)                                                         -0.3395     0.9742  -0.348    0.729
#log10(hsc70[hsc70$gene == "ENSG00000109971", ]$GenerationLength_d)   0.1561     0.2818   0.554    0.581

#Residual standard error: 0.7342 on 70 degrees of freedom
#(17 observations deleted due to missingness)
#Multiple R-squared:  0.004365,	Adjusted R-squared:  -0.009858 

summary(lm(hsc70[hsc70$gene == 'ENSG00000009335',]$dN.dS ~ log10(hsc70[hsc70$gene == 'ENSG00000009335',]$GenerationLength_d)))
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)                                                        -0.006234   0.120991  -0.052    0.959
#log10(hsc70[hsc70$gene == "ENSG00000009335", ]$GenerationLength_d)  0.030153   0.035361   0.853    0.396

#Residual standard error: 0.1111 on 86 degrees of freedom
#(19 observations deleted due to missingness)
#Multiple R-squared:  0.008384,	Adjusted R-squared:  -0.003147 

summary(lm(hsc70$dN.dS ~ log10(hsc70$GenerationLength_d) * hsc70$gene))
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)                                               -0.006234   0.542973  -0.011    0.991
#log10(hsc70$GenerationLength_d)                            0.030153   0.158691   0.190    0.850
#hsc70$geneENSG00000109971                                 -0.333258   0.855983  -0.389    0.698
#log10(hsc70$GenerationLength_d):hsc70$geneENSG00000109971  0.125970   0.248645   0.507    0.613

#Residual standard error: 0.4987 on 156 degrees of freedom
#(36 observations deleted due to missingness)
#Multiple R-squared:  0.01482,	Adjusted R-squared:  -0.004129



library(ggplot2)

ggplot(hsps[hsps$gene == 'ENSG00000096384' | hsps$gene == controls[controls$hsp == 'ENSG00000096384',]$control_gene,], 
       aes(x = log10(GenerationLength_d), y = dN.dS, color = gene))+
  geom_smooth(method = lm, size = 1.9)+
  geom_point()+
  theme_bw()+
  ylab('dN/dS')+
  xlab('log10(Generation length, days)')+
  scale_color_discrete(name = "Gene", labels = c('HSP90', 'COPS3'))+
  theme(axis.title = element_text(size = 27),
        axis.text = element_text(size = 22),
        legend.title = element_text(size = 25, face = 'bold'), legend.text = element_text(size = 20))+
  labs(colour="Gene")
  
  
  

ggplot(hsps[hsps$gene == 'ENSG00000109971' | hsps$gene == 'ENSG00000009335',],aes(x = log10(GenerationLength_d), y = dN.dS, color = role))+
  geom_smooth(method = lm, size = 1.6)+
  geom_point()+
  theme_bw()+
  ylab('dN/dS')+
  xlab('log10(Generation length, days)')+
  scale_color_discrete(name = "Gene", labels = c('HSÑ70', 'UBE3C'))+
  theme(axis.title = element_text(size = 19),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 17, face = 'bold'), legend.text = element_text(size = 14))



for (gene in unique(hsps[hsps$role == 'hsp',]$gene)) {
  if (controls[controls$hsp == gene,]$control_gene %in% unique(hsps[hsps$role == 'nonhsp',]$gene)){
  
  gg = ggplot(hsps[hsps$gene == gene | hsps$gene == controls[controls$hsp == gene,]$control_gene,], 
         aes(x = log10(GenerationLength_d), y = dN.dS, color = role))+
    geom_smooth(method = lm)+
    geom_point()+
    theme_bw()
  ggsave(paste('../../Body/4_Figures/hsps.nonhsps.paml.dN.dS.vs.genlen.', gene, '.png', sep = ''), gg)
  }
}

summary(lm(hsps[hsps$gene == 'ENSG00000109971',]$dN.dS ~ log10(hsps[hsps$gene == 'ENSG00000109971',]$GenerationLength_d)))
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)
#(Intercept)                                                       -0.3395     0.9742  -0.348    0.729
#log10(hsps[hsps$gene == "ENSG00000109971", ]$GenerationLength_d)   0.1561     0.2818   0.554    0.581


summary(lm(hsps[hsps$gene == 'ENSG00000009335',]$dN.dS ~ log10(hsps[hsps$gene == 'ENSG00000009335',]$GenerationLength_d)))



slopes <- c()
p_val_slope <- c()
intercept <- c()
p_val_intercept <- c()
R_sq <- c()
R_sq_adj <- c()
number_of_species <- c()
genes <- unique(unique(hsps[hsps$role == 'hsp',]$gene))
residual_std_err <- c()
residuals_hsps <- data.frame(Species = unique(unique(hsps[hsps$role == 'hsp',]$Species)))



for (gene in genes){
  a <- na.omit(hsps[hsps$gene == gene,])
  
  species <- unique(a[,'Species'])
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
  residuals_hsps <- merge(residuals_hsps, res, by = 'Species', all = T)
}

lm_hsps <- data.frame(genes, slopes, intercept, p_val_slope, p_val_intercept, number_of_species, R_sq, R_sq_adj, residual_std_err)

hsp_groups <- read.table('../../../Molecular_evolution/Body/1_Raw/human.hsp.genes.txt', header = T, sep = '\t')
lm_hsps <- merge(lm_hsps, hsp_groups[,c(2,10,13)], by.x = 'genes', by.y = 'Ensembl.gene.ID', all.x = T)

ggplot(lm_hsps, aes(x = slopes, y = p_val_slope, col = Group.name, label = Approved.symbol))+
  geom_point()+geom_text()


length(lm_hsps[lm_hsps$p_val_slope < 0.05,'genes'])/length(lm_hsps[,'genes']) #0.1428571



slopes <- c()
p_val_slope <- c()
intercept <- c()
p_val_intercept <- c()
R_sq <- c()
R_sq_adj <- c()
number_of_species <- c()
genes <- unique(hsps[hsps$role == 'nonhsp',]$gene)
residual_std_err <- c()
residuals_nonhsps <- data.frame(Species = unique(unique(hsps[hsps$role == 'nonhsp',]$Species)))



for (gene in genes){
  a <- na.omit(hsps[hsps$gene == gene,])
  
  species <- unique(a[,'Species'])
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
  residuals_nonhsps <- merge(residuals_nonhsps, res, by = 'Species', all = T)
}

lm_nonhsps <- data.frame(genes, slopes, intercept, p_val_slope, p_val_intercept, number_of_species, R_sq, R_sq_adj, residual_std_err)

length(lm_nonhsps[lm_nonhsps$p_val_slope < 0.05,'genes'])/length(lm_nonhsps[,'genes']) #0.2

hsps <- hsps[hsps$gene %in% c('ENSG00000109971', 'ENSG00000009335', 'ENSG00000096384', 'ENSG00000141030'),]

hsps$hsp <- ifelse(hsps$role == 'hsp', 1, 0)
hsps$nonhsp <- ifelse(hsps$role == 'nonhsp', 1, 0)
summary(lm(hsps$dN.dS ~ log10(hsps$GenerationLength_d) + hsps$hsp+ hsps$nonhsp))




#################redbook##############################
redbook <- read.csv('../../Body/1_Raw/IUCN.csv', sep = ';', header = T)
redbook$specieName <- gsub(' ', '_', redbook$specieName)
redbook <- redbook[redbook$class == 'MAMMALIA',]
redbook <- redbook[,c(3,14)]

residuals_hsps <- merge(residuals_hsps, redbook, by.x  = 'Species', by.y = 'specieName', all.x = T)
residuals_hsps$category <- factor(residuals_hsps$category , levels=c("LC", 'NT', 'VU', 'EN', 'CR', 'DD'))
residuals_hsps <- residuals_hsps[!is.na(residuals_hsps$category),]


residuals_nonhsps <- merge(residuals_nonhsps, redbook, by.x  = 'Species', by.y = 'specieName', all.x = T)
residuals_nonhsps$category <- factor(residuals_nonhsps$category , levels=c("LC", 'NT', 'VU', 'EN', 'CR', 'DD'))
residuals_nonhsps <- residuals_nonhsps[!is.na(residuals_nonhsps$category),]

ggplot(residuals_hsps, aes(y = ENSG00000096384, x = category, fill = category))+
  geom_boxplot()+
  theme_bw()+
  ggtitle('hsp90')+
  ylab('residuals')

ggplot(residuals_nonhsps, aes(y = ENSG00000141030, x = category, fill = category))+
  geom_boxplot()+
  theme_bw()+
  ggtitle('hsp90_control (COSPS3)')+
  ylab('residuals')

ggplot(residuals_hsps, aes(y = ENSG00000109971, x = category, fill = category))+
  geom_boxplot()+
  theme_bw()+
  ggtitle('hsc70')+
  ylab('residuals')

ggplot(residuals_nonhsps, aes(y = ENSG00000009335, x = category, fill = category))+
  geom_boxplot()+
  theme_bw()+
  ggtitle('hsc70')+
  ylab('residuals')


ggplot(residuals_hsps, aes(y = ENSG00000096384, x = category, fill = category))+
  geom_boxplot(outlier.shape = NA)+
  ylim(c(-0.05,0.1))+
  theme_bw()+
  ggtitle('hsp90')+
  ylab('residuals')

residuals_all <- merge(residuals_hsps, residuals_nonhsps, by = 'Species', all = T)
residuals_all <- residuals_all %>%
  pivot_longer(colnames(residuals_all)[-c(1,16,27)], names_to = 'gene', values_to = 'residuals')

residuals_all$binar_cat <- ifelse(residuals_all$category.x == 'LC', 'LC', 'others')


ggplot(residuals_all[(residuals_all$gene == 'ENSG00000096384') | (residuals_all$gene == 'ENSG00000109971'),],
       aes(y = residuals, x = gene, fill = binar_cat))+
  geom_boxplot()+
  theme_bw()+
  ggtitle('hsp90')+
  ylab('residuals')
  #ylim(c(-0.5, 0.3))


ggplot(residuals_all[(residuals_all$gene == 'ENSG00000096384'),],
       aes(y = residuals, x = gene, fill = binar_cat))+
  geom_boxplot()+
  theme_bw()+
  ggtitle('hsp90')+
  ylab('residuals')


ggplot(residuals_all[(residuals_all$gene == 'ENSG00000096384') | (residuals_all$gene == 'ENSG00000141030') | (residuals_all$gene == 'ENSG00000109971'),],
       aes(y = residuals, x = binar_cat, fill = gene))+
  geom_boxplot()+
  theme_bw()+
  ylab('residuals')



############PGLS##############
library(ape)
library(geiger)
library(caper)

tree_hsp90 <- read.tree("../../Body/2_Derived/hsps/tree/ENSG00000096384_HSP90AB1_NT.fasta.nw")
hsp90 <- hsps[hsps$gene == 'ENSG00000096384',]
row.names(hsp90) = hsp90$Species
hsp90_2 <- na.omit(hsp90[,c(1,6,14)])

tree_w = treedata(tree_hsp90, hsp90_2, sort=T, warnings=T)$phy
data<-as.data.frame(treedata(tree_w, hsp90_2, sort=T, warnings=T)$data)


data$Species = as.character(data$Species)

data$dN.dS = as.numeric(as.character(data$dN.dS))
data$GenerationLength_d = as.numeric(as.character(data$GenerationLength_d))



hsp90_comp = comparative.data(tree_w, data, Species, vcv=TRUE)
phylosig(tree_w, log10(data$GenerationLength_d), method = "lambda", test = TRUE) #0.999934
phylosig(tree_w, data$dN.dS, method = "lambda", test = TRUE) #6.6107e-05


model = pgls(dN.dS ~ log10(GenerationLength_d), hsp90_comp, lambda="ML")
summary(model)
#Branch length transformations:
  
#  kappa  [Fix]  : 1.000
#lambda [ ML]  : 0.000
#lower bound : 0.000, p = 1    
#upper bound : 1.000, p = 0.0084484
#95.0% CI   : (NA, 0.851)
#delta  [Fix]  : 1.000

#Coefficients:
#  Estimate  Std. Error t value Pr(>|t|)
#(Intercept)                0.03412526  0.11403867  0.2992   0.7656
#log10(GenerationLength_d) -0.00082385  0.03312995 -0.0249   0.9802

#Residual standard error: 0.007133 on 70 degrees of freedom
#Multiple R-squared: 8.834e-06,	Adjusted R-squared: -0.01428 
#F-statistic: 0.0006184 on 1 and 70 DF,  p-value: 0.9802 

par(mfrow=c(2,2))
plot(model)

par(mfrow=c(1,1))
plot(data$GenerationLength_d, data$dN.dS)
abline(summ$coefficients[1:2])


tree_hsp90_cont <- read.tree("../../Body/2_Derived/nonhsp/tree/ENSG00000141030_COPS3_NT.fasta.nw")
row.names(hsp90_control) = hsp90_control$Species
hsp90_control_2 <- na.omit(hsp90_control[,c(1,6,14)])

tree_w = treedata(tree_hsp90_cont, hsp90_control_2, sort=T, warnings=T)$phy
data<-as.data.frame(treedata(tree_w, hsp90_control_2, sort=T, warnings=T)$data)


data$Species = as.character(data$Species)

data$dN.dS = as.numeric(as.character(data$dN.dS))
data$GenerationLength_d = as.numeric(as.character(data$GenerationLength_d))

MutComp = comparative.data(tree_w, data, Species, vcv=TRUE)

model = pgls(dN.dS ~ log10(GenerationLength_d), MutComp, lambda="ML")
summary(model)
#Branch length transformations:
  
#  kappa  [Fix]  : 1.000
#lambda [ ML]  : 0.277
#lower bound : 0.000, p = 0.016206
#upper bound : 1.000, p = < 2.22e-16
#95.0% CI   : (0.044, 0.547)
#delta  [Fix]  : 1.000

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)               -0.62687    0.68161 -0.9197   0.3604
#log10(GenerationLength_d)  0.26346    0.19690  1.3381   0.1845

#Residual standard error: 0.04125 on 84 degrees of freedom
#Multiple R-squared: 0.02087,	Adjusted R-squared: 0.009213 
#F-statistic:  1.79 on 1 and 84 DF,  p-value: 0.1845











