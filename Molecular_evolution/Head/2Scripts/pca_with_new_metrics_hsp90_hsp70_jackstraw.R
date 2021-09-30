rm(list = ls(all = TRUE))
library('jackstraw')

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
genes_jackstraw_pca <- jackstraw_pca(scale(t(genes_without_na)), r = 1, s = 3, B = 200)
genes_jackstraw_pca$p.value
result <- data.frame(feature = colnames(genes_without_na), p.val <- genes_jackstraw_pca$p.value)
