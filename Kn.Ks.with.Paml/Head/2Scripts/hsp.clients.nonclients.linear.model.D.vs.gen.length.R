rm(list=ls(all=TRUE))


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

clients[,-1] <- log10(clients[,-1]+0.0001)
nonclients[,-1] <- log10(nonclients[,-1]+0.0001)


#perform linear model for clients
slopes <- c()
p_val_slope <- c()
intercept <- c()
p_val_intercept <- c()
R_sq <- c()
R_sq_adj <- c()
number_of_species <- c()
genes <- colnames(clients)[-c(1,295)]
residual_std_err <- c()
residuals_clients <- data.frame(Species = clients$Species)


for (gene in genes) {
  
  a <- na.omit(clients[,c("Species", gene, 'Generation_Length')])
  
  species <- a[,'Species']
  number_of_species <- c(number_of_species, length(species))
  
  lg <- lm(a[,gene]~a[,'Generation_Length'], data = a)
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

write.table(lm_clients, '../../Body/3_Results/hsp.clients.linear.model.D.vs.gen.length.txt')







##########################################
############################################


#perform linear model for nonclients
slopes <- c()
p_val_slope <- c()
intercept <- c()
p_val_intercept <- c()
R_sq <- c()
R_sq_adj <- c()
number_of_species <- c()
genes <- colnames(nonclients)[-c(1,387)]
residual_std_err <- c()
residuals_nonclients <- data.frame(Species = nonclients$Species)


for (gene in genes) {
  
  a <- na.omit(nonclients[,c("Species", gene, 'Generation_Length')])
  
  species <- a[,'Species']
  number_of_species <- c(number_of_species, length(species))
  
  lg <- lm(a[,gene]~a[,'Generation_Length'], data = a)
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
  residuals_nonclients <- merge(residuals_nonclients, res, by = 'Species', all = T)
  
}


#table with linear regression parameters
lm_nonclients <- data.frame(genes, slopes, intercept, p_val_slope, p_val_intercept, number_of_species, R_sq, R_sq_adj, residual_std_err)

write.table(lm_nonclients, '../../Body/3_Results/hsp.nonclients.linear.model.D.vs.gen.length.txt')


lm_clients$client <- T
lm_nonclients$client <- F
result <- rbind(lm_clients, lm_nonclients)


pdf('../../Body/4_Figures/hsp.clients.nonclients.linear.model.D.vs.gen.length.pdf')
library(ggplot2)
ggplot(data = result, aes(x = slopes, y = p_val_slope, color = client))+
  geom_point()

ggplot(result, aes(y = slopes, color = client))+
  geom_boxplot()

ggplot(result[result$p_val_slope <= 0.05,], aes(y = slopes, x = client, fill = client))+
  geom_violin(trim = F, width = 0.6)+
  geom_boxplot(width = 0.1)+
  theme_bw()+
  ggtitle('p_val_slope < 0.05')+
  theme(legend.position="none", axis.text=element_text(size=15), axis.title=element_text(size=14,face="bold"))+
  xlab('')+scale_x_discrete(breaks=c(FALSE, TRUE),
                            labels=c('nonclients', 'clients'))

ggplot(result, aes(y = intercept, color = client))+
  geom_boxplot()

ggplot(data = result, aes(x = intercept, y = slopes, color = client))+
  geom_point()
dev.off()


##################################################################
#########Species-specific-residuals
################################################################

residuals_clients$mean <- apply(residuals_clients[,-1], 1, mean, na.rm=T)
residuals_nonclients$mean <- apply(residuals_nonclients[,-1], 1, mean, na.rm=T)

residuals_clients$count_positives <- apply(residuals_clients[,-c(1,295)], 1, function(x) {sum(x>=0, na.rm = T)})
residuals_nonclients$count_positives <- apply(residuals_nonclients[,-c(1,387)], 1, function(x) {sum(x>=0, na.rm = T)})

residuals_clients$count_negatives <- apply(residuals_clients[,-c(1,295, 296)], 1, function(x) {sum(x<0, na.rm = T)})
residuals_nonclients$count_negatives <- apply(residuals_nonclients[,-c(1,387,388)], 1, function(x) {sum(x<0, na.rm = T)})

residuals_clients[residuals_clients$count_positives > residuals_clients$count_negatives, 'Species'] #0
residuals_nonclients[residuals_nonclients$count_positives > residuals_nonclients$count_negatives, 'Species'] #0

residuals_clients[residuals_clients$mean > 0, 'Species'] #49
residuals_nonclients[residuals_nonclients$mean > 0, 'Species'] #46

residuals_sp <- merge(residuals_clients[,c(1,295,296,297)], residuals_nonclients[,c(1,387,388,389)], by = 'Species', all = T)
library(tidyverse)
library(tidyr)

residuals_sp_mean <- residuals_sp %>%
  pivot_longer(c('mean.x', 'mean.y'), names_to = "group", values_to = "mean")

residuals_sp_pos <- residuals_sp %>%
  pivot_longer(c('count_positives.x', 'count_positives.y'), names_to = "group", values_to = "count_positives")

residuals_sp_neg <- residuals_sp %>%
  pivot_longer(c('count_negatives.x', 'count_negatives.y'), names_to = "group", values_to = "count_negatives")

residuals_sp_mean$client <- ifelse(residuals_sp_mean$group == 'mean.x', T, F)
residuals_sp_pos$client <- ifelse(residuals_sp_pos$group == 'count_positives.x', T, F)
residuals_sp_neg$client <- ifelse(residuals_sp_neg$group == 'count_negatives.x', T, F)

residuals_sp <- merge(residuals_sp_mean[,c(1,7,8)], residuals_sp_pos[,c(1,7,8)], by = c('Species', 'client'), all = T)
residuals_sp <- merge(residuals_sp, residuals_sp_neg[,c(1,7,8)], by = c('Species', 'client'), all = T)


residuals_sp$positives_prop <- residuals_sp$count_positives / (residuals_sp$count_positives + residuals_sp$count_negatives)
residuals_sp$pos_to_neg <- residuals_sp$count_positives / residuals_sp$count_negatives

wilcox.test(residuals_sp[residuals_sp$client,]$mean, residuals_sp[!residuals_sp$client,]$mean, paired = T)
#V = 2870, p-value = 0.701

wilcox.test(residuals_sp[residuals_sp$client,]$positives_prop, residuals_sp[!residuals_sp$client,]$positives_prop, paired = T)
#V = 3927, p-value = 0.004976

wilcox.test(residuals_sp[residuals_sp$client,]$pos_to_neg, residuals_sp[!residuals_sp$client,]$pos_to_neg, paired = T)
#V = 3937, p-value = 0.004528

ggplot(residuals_sp, aes(x = Species, y = positives_prop, color = client))+
  geom_point()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90), legend.position = 'top')+
  ylab('proportion of positive residuals')+
  ggtitle('Species-specific residuals from linear model')


ggplot(residuals_sp, aes(x = Species, y = pos_to_neg, color = client))+
  geom_point()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90), legend.position = 'top')+
  ylab('positive to negative residuals ratio')+
  ggtitle('Species-specific residuals from linear model')


residuals_sp <- merge(residuals_sp, clients[,c(1,295)], by = 'Species', all.x = T)

#corr for all
cor.test(residuals_sp$positives_prop, residuals_sp$Generation_Length, method = 'spearman')
#rho = -0.2339365, S = 2130600, p-value = 0.0004965

#corr for clients
cor.test(residuals_sp[residuals_sp$client,]$positives_prop, residuals_sp[residuals_sp$client,]$Generation_Length, method = 'spearman')
#rho = -0.2641133, S = 272821, p-value = 0.005519

#corr for nonclients
cor.test(residuals_sp[!residuals_sp$client,]$positives_prop, residuals_sp[!residuals_sp$client,]$Generation_Length, method = 'spearman')
#rho = -0.2095157, S = 261038, p-value = 0.02878



glm_model <- glm(residuals_sp$positives_prop ~ residuals_sp$Generation_Length + residuals_sp$client)
summary(glm_model)

#Deviance Residuals: 
#  Min        1Q    Median        3Q       Max  
#-0.11758  -0.04680  -0.01327   0.03949   0.23401  

#Coefficients:
#                             Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                     0.31575    0.04564   6.919 5.14e-11 ***
#  residuals_sp$Generation_Length -0.04275    0.01326  -3.223  0.00146 ** 
#  residuals_sp$clientTRUE         0.01019    0.00878   1.161  0.24700    
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#(Dispersion parameter for gaussian family taken to be 0.004201343)

#Null deviance: 0.95260  on 217  degrees of freedom
#Residual deviance: 0.90329  on 215  degrees of freedom
#AIC: -569.34

#Number of Fisher Scoring iterations: 2












