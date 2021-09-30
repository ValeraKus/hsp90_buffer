rm(list = ls(all = TRUE))

gnomad <- read.table('../../Body/1_Raw/supplementary_dataset_11_full_constraint_metrics.tsv', sep = '\t', header = T)

metrics_table <- read.table('../../Body/1_Raw/gencode.v25.annotation.gtf.Genes.Shet.pLI.FIS.RVIS.GHIS.KnKs.GC.BrainSpecificRanking.Branch', sep = '\t', header = T)

gnomad <- gnomad[gnomad$canonical == 'true',]

gnomad <- gnomad[,c('gene_id', 'oe_lof_upper_bin', 'p')]

new_metrics_table <- merge(metrics_table, gnomad, by.x = 'EnsemblId', by.y = 'gene_id', all.x = T, all.y = T)


write.table(new_metrics_table, '../../Body/2_Derived/gencode.v25.annotation.gtf.Genes.Shet.pLI.FIS.RVIS.GHIS.KnKs.GC.BrainSpecificRanking.Branch.gnomad.oe_lof_upper_bin.p')
