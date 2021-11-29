rm(list = ls(all=T))
library(ggplot2)


setwd('~/Desktop/Work/hsp90_buffer/Human.genetic.variation/head/2scripts/')

hsp90_phenotypes_files <- list.files('../../body/2derived/hsp90.PGS.phenotypes/')

phenotypes_without_snps <- c()

phenotypeID <- c()
num_of_snps <- c()
num_of_pos_snsps <- c()
proportion_of_pos_snps <- c()
mean_beta <- c()
pos_to_neg_effect_weight <- c()

diff_cols <- c()

pdf('../../body/4figures/pgs.effect.size.distribution.for.each.phenotypes.pdf')
for (f in hsp90_phenotypes_files){
  
  phenotype = gsub('.txt.snps.txt', '', gsub('ENSG[0-9]*\\.', '', f))
  
  filename <-  paste('../../body/2derived/hsp90.PGS.phenotypes/', f, sep ='')
  
  if (file.info(filename)$size != 0){
    
  pgs.table <- read.table(filename, header = F, sep = '\t')
  if  ((ncol(pgs.table) >= 6) && (typeof(pgs.table[,6]) != "character")){
    
    #colnames(pgs.table) <- c('SnpID', 'chr', 'pos', 'effect_allele', 'reference_allele', 'effect_weight', 'weight_type')

    
    boxplot(pgs.table[,6], ylab = 'effect weight', main = phenotype)
    
    phenotypeID <- c(phenotypeID, phenotype)
    num_of_snps <- c(num_of_snps, nrow(pgs.table))
    num_of_pos_snsps <- c(num_of_pos_snsps, nrow(pgs.table[pgs.table[,6] >0,]))
    proportion_of_pos_snps <- c(proportion_of_pos_snps, nrow(pgs.table[pgs.table[,6] >0,])/nrow(pgs.table))
    mean_beta <- c(mean_beta, mean(pgs.table[,6]))
    pos_to_neg_effect_weight <- c(pos_to_neg_effect_weight, 
                                  sum(pgs.table[pgs.table[,6] > 0, 6])/abs(sum(pgs.table[pgs.table[,6] < 0, 6])))
  }
  else{diff_cols <- c(diff_cols, phenotype)}}
  
  else {
    phenotypes_without_snps <- c(phenotypes_without_snps, phenotype)
  }}
  
  
  
dev.off()  
  

pgs.summary.table <- data.frame(phenotypeID, num_of_snps, num_of_pos_snsps, proportion_of_pos_snps, mean_beta, pos_to_neg_effect_weight)
write.table(pgs.summary.table, '../../body/2derived/pgs.phenotypes.hsp90.association.summary.txt')
write(phenotypes_without_snps, '../../body/3results/pgs.phenotypes.with.no.association.hhp90.txt')

pgs.summary.table <- read.table('../../body/2derived/pgs.phenotypes.hsp90.association.summary.txt')

ggplot(pgs.summary.table, aes(mean_beta, log10(num_of_snps), color = proportion_of_pos_snps))+
  geom_point()

ggplot(pgs.summary.table, aes(pos_to_neg_effect_weight, num_of_snps, color = proportion_of_pos_snps))+
  geom_point()

ggplot(pgs.summary.table, aes(x = num_of_snps, color = mean_beta > 0))+
  geom_boxplot(outlier.shape = NA)+
  xlim(0,50)+
  theme_bw()

ggplot(pgs.summary.table, aes(x = num_of_snps, color = pos_to_neg_effect_weight > 1))+
  geom_boxplot(outlier.shape = NA)+
  xlim(0,50)+
  theme_bw()

ggplot(pgs.summary.table, aes(mean_beta, fill = num_of_snps))+
  geom_histogram(bins=60)

wilcox.test(pgs.summary.table$mean_beta, mu=0, alternative = 'greater')

##################
######Phenotype Ids decoding
##################


phenotypes_decode <- read.csv('../../body/1raw/pgs_scores_data.csv')
phenotypes_decode <- phenotypes_decode[,c(1,3)]
phenotypes_decode$Polygenic.Score.ID...Name <- gsub('\\(.*\\)', '', phenotypes_decode$Polygenic.Score.ID...Name)
colnames(phenotypes_decode) <- c('phenotypeID', 'Trait')

pgs.summary.table <- merge(pgs.summary.table, phenotypes_decode, by = 'phenotypeID', all.x = T)


####################
#####Look at top 19 (num of snps > 1000) traits which hsp90 is associated with (by number of snps)
#####################

top <- pgs.summary.table[order(pgs.summary.table$num_of_snps, decreasing = T) & (pgs.summary.table$num_of_snps>1000), c(1,2,7)]

#"Coronary artery disease"              "Ischaemic stroke"                    
#"Coronary heart disease"               "Type 2 diabetes"                     
# "Atrial fibrillation"                  "Breast cancer"                       
# "Prostate cancer"                      "Atrial fibrillation"                 
# "Chronic kidney disease"               "Type 2 diabetes"                     
# "Diabetic retinopathy"                 "Estimated glomerular filtration rate"
# "Major depressive disorder"            "Insomnia"                            
# "Headache"                             "Body mass index"                     
# "Coronary artery disease"              "Atrial fibrillation"                 
# "Type 2 diabetes" 

write.table(top, '../../body/3results/pgs.phenotypes.hsp90.association.top.phenotypes.num.of.snsps.more.1000.txt')




