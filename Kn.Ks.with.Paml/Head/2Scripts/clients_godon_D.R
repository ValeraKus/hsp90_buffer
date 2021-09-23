rm(list=ls(all=TRUE))

tables <- list.files('../../Body/2_Derived/Godon/clients/parsed/')
whole_table <- setNames(data.frame(matrix(ncol = 1, nrow = 0)), c('Species'))
for (file in tables) {
  table <- paste('../../Body/2_Derived/Godon/clients/parsed/', file, sep = '')
  table <- read.table(table, sep = '\t', header = T)
  #table <- table[table$Type == '1-to-1View Gene Tree',c("species",'D')]
  ens <- gsub('_[A-Z, 0-9]*_[A-Z]*\\.txt$','' ,gsub('parsed_', '', file))
  colnames(table) <- c('Species', paste('D_', ens, sep = ''))
  whole_table <- merge(whole_table, table, by = 'Species', all = T)
}


glenght_table <- read.table('../../Body/2_Derived/GenerationLenghtforMammals.csv', sep = ',', header = T)
#glenght_table <- glenght_table[, c('Scientific_name', 'GenerationLength_d')]
#Species <- gsub('\\)', '',gsub('^([a-z]|[A-Z]|\ |-|\'|[0-9])*\\(', '', whole_table$Species))
library(stringr)
whole_table$Species <-  word(whole_table$Species, 1, 2, sep = '_')
glenght_table$Scientific_name <- gsub(' ', '_', glenght_table$Scientific_name)
whole_table$Generation_Length <- NA

for (species in whole_table$Species){
  if (species %in% glenght_table$Scientific_name){
    whole_table[whole_table$Species == species,'Generation_Length'] <-  glenght_table[glenght_table$Scientific_name == species, 'GenerationLength_d']
  } 
  else {
    if (word(species, 1, sep = '_') %in% glenght_table$Genus){
      mean_glen = mean(glenght_table[glenght_table$Genus == word(species, 1, sep = '_'), 'GenerationLength_d'])
      whole_table[whole_table$Species == species,'Generation_Length'] <- mean_glen
    }
  }
}



whole_table <- whole_table[!is.na(whole_table$Generation_Length),]

write.table(whole_table, '../../Body/2_Derived/clients_godon_D.txt')
