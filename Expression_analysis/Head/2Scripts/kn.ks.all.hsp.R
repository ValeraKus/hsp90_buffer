rm(list=ls(all=TRUE))

tables <- list.files('../../Body/1_Raw/Ensemble_Compara_all_hsp/')
whole_table <- setNames(data.frame(matrix(ncol = 1, nrow = 0)), c('Species'))
for (file in tables) {
  table <- paste('../../Body/1_Raw/Ensemble_Compara_all_hsp/', file, sep = '')
  table <- read.table(table, sep = ',', header = T)
  table <- table[table$Type == '1-to-1View Gene Tree',c("Species",'dN.dS')]
  ens <- gsub('.csv', '', gsub('^.*orthologues-ComparaOrthologs-Homo_sapiens_Gene_Compara_Ortholog_', '',file))
  colnames(table) <- c('Species', paste('dN.dS_', ens, sep = ''))
  whole_table <- merge(whole_table, table, by = 'Species', all = T)
}


glenght_table <- read.table('../../Body/2_Derived/GenerationLenghtforMammals.csv', sep = ',', header = T)
#glenght_table <- glenght_table[, c('Scientific_name', 'Calculated_GL_d')]
Species <- gsub('\\)', '',gsub('^([a-z]|[A-Z]|\ |-|\'|[0-9])*\\(', '', whole_table$Species))
library(stringr)
whole_table$Species <- word(Species, 1, 2)

whole_table$Generation_Length <- NA

for (species in whole_table$Species){
  if (species %in% glenght_table$Scientific_name){
    whole_table[whole_table$Species == species,'Generation_Length'] <-  glenght_table[glenght_table$Scientific_name == species, 'GenerationLength_d']
  } 
}
whole_table <- whole_table[!is.na(whole_table$Generation_Length),]


write.table(whole_table, '../../Body/2_Derived/kn.ks.all.hsp.txt')
