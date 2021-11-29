rm(list = ls(all=T))
library(ggplot2)


setwd('~/Desktop/Work/hsp90_buffer/Human.genetic.variation/head/2scripts/')

egen <- read.table(gzfile('../../body/1raw/GTEx_Analysis_v8_eQTL/Adipose_Subcutaneous.v8.egenes.txt.gz'), sep = '\t', header = 1)
length(egen$gene_id) #24665
length(unique(egen$gene_id)) #24665
egen.hsp90 <- egen[grepl('ENSG00000096384' ,egen$gene_id),] #top variant of each egene
egen.hsp90$num_var #7563
summary(egen)

##File extensions: *.signif_variant_gene_pairs.txt.gz, *.sqtl_signifpairs.txt.gz
##Column headers:
  
##  variant_id:               variant ID in the format {chr}_{pos_first_ref_base}_{ref_seq}_{alt_seq}_b38
##gene_id:                  GENCODE/Ensembl gene ID
##tss_distance:             distance between variant and transcription start site. Positive when variant is downstream of the TSS, negative otherwise
##ma_samples:               number of samples carrying the minor allele
##ma_count:                 total number of minor alleles across individuals
##maf:                      minor allele frequency observed in the set of donors for a given tissue
#pval_nominal:             nominal p-value
##slope:                    regression slope
##slope_se:                 standard error of the regression slope
##pval_nominal_threshold:   nominal p-value threshold for calling a variant-gene pair significant for the gene
##min_pval_nominal:         smallest nominal p-value for the gene
##pval_beta:                beta-approximated permutation p-value for the gene
sign_pairs <- read.table(gzfile('../../body/1raw/GTEx_Analysis_v8_eQTL/Adipose_Subcutaneous.v8.signif_variant_gene_pairs.txt.gz'), sep = '\t', header = 1)
sign_pairs.hsp90 <- sign_pairs[grepl('ENSG00000096384', sign_pairs$gene_id),]


sign_pairs_files <- list.files('../../body/1raw/GTEx_Analysis_v8_eQTL/')
sign_pairs_files <- sign_pairs_files[grepl('.signif_variant_gene_pairs.txt.gz', sign_pairs_files)]

hsp90_sign_pairs <- data.frame()


for (tiss_file in sign_pairs_files){
  
  tissue = gsub('.v8.signif_variant_gene_pairs.txt.gz', '', tiss_file)
  df <- read.table(gzfile(paste('../../body/1raw/GTEx_Analysis_v8_eQTL/', tiss_file,sep = '')), sep = '\t', header = 1)
  df$tissue <- tissue
  df <- df[grepl('ENSG00000096384', df$gene_id),]
  hsp90_sign_pairs <- rbind(hsp90_sign_pairs, df)
  }

#209
write.table(hsp90_sign_pairs, '../../body/2derived/Gtex.signif.variant.gene.pairs.hsp90.txt')

hsp90_sign_pairs <- read.table('../../body/2derived/Gtex.signif.variant.gene.pairs.hsp90.txt')
length(unique(hsp90_sign_pairs$variant_id)) #126

pdf('../../body/4figures/Gtex.signif.variant.gene.pairs.hsp90.pdf')
hist(hsp90_sign_pairs$maf)

ggplot(hsp90_sign_pairs, aes(x = slope>0, y = maf, fill = slope>0))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = 'None')

ggplot(hsp90_sign_pairs, aes(x = tissue, fill = tissue))+
  geom_bar()+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size = 10), legend.position = 'None', axis.title.x = element_blank())


ggplot(hsp90_sign_pairs, aes(slope, maf, color = tissue))+
  geom_point()+
  theme_bw()

sum(hsp90_sign_pairs$slope > 0)/length(hsp90_sign_pairs$slope)
sum(hsp90_sign_pairs$slope < 0)

ggplot(hsp90_sign_pairs, aes(x = reorder(tissue, slope), y = slope, fill = tissue))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size = 10), legend.position = 'None', axis.title.x = element_blank())

ggplot(hsp90_sign_pairs, aes(x = reorder(tissue, maf), y = maf, fill = tissue))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size = 10), legend.position = 'None', axis.title.x = element_blank())


ggplot(hsp90_sign_pairs, aes(y = slope*maf, x = reorder(tissue, slope*maf), fill = tissue))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size = 10), legend.position = 'None', axis.title.x = element_blank())


tissue_longivity <- read.table('../../../Expression_analysis/Body/1_Raw/hsp90_human_tissues_turnover_rate_expression_level.csv', sep =',', header = 1)
hsp90_sign_pairs$Turnover <- NA
unique(hsp90_sign_pairs$tissue)

fast <- c('Colon_Sigmoid', 'Esophagus_Mucosa', 'Small_Intestine_Terminal_Ileum')
intermediate <- c('Liver', 'Skin_Not_Sun_Exposed_Suprapubic', 'Skin_Sun_Exposed_Lower_leg', 'Testis')
slow <- c('Adipose_Subcutaneous', 'Artery_Aorta', 'Artery_Tibial', 'Brain_Cortex', 'Brain_Nucleus_accumbens_basal_ganglia', 'Brain_Spinal_cord_cervical_c-1')

hsp90_sign_pairs[hsp90_sign_pairs$tissue %in% fast,]$Turnover <- 'fast'
hsp90_sign_pairs[hsp90_sign_pairs$tissue %in% intermediate,]$Turnover <- 'intermediate'
hsp90_sign_pairs[hsp90_sign_pairs$tissue %in% slow,]$Turnover <- 'slow'

ggplot(hsp90_sign_pairs[!is.na(hsp90_sign_pairs$Turnover),], aes(tissue, fill = Turnover))+
  geom_bar()+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size = 10), axis.title.x = element_blank())



ggplot(hsp90_sign_pairs[!is.na(hsp90_sign_pairs$Turnover),], aes(x = Turnover, y = slope, fill = Turnover))+
  geom_violin()+
  geom_boxplot(width=0.05)
  theme_bw()+
  theme(axis.text.x = element_text(size = 15), legend.position = 'None')

ggplot(hsp90_sign_pairs[!is.na(hsp90_sign_pairs$Turnover),], aes(x = Turnover, y = maf, fill = Turnover))+
  geom_boxplot()
  theme_bw()+
  theme(axis.text.x = element_text(size = 15), legend.position = 'None')

ggplot(hsp90_sign_pairs[!is.na(hsp90_sign_pairs$Turnover),], aes(x = Turnover, y = maf*slope, fill = Turnover))+
  geom_violin()+
  geom_boxplot(width=0.05)
  theme_bw()+
  theme(axis.text.x = element_text(size = 15), legend.position = 'None')





dev.off()

