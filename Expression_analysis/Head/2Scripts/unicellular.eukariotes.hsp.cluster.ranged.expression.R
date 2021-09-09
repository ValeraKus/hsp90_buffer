rm(list=ls(all=TRUE))

data <- read.csv('../../Body/1_Raw/unicelullar_expression.csv')
uniref_cluster <- 'UniRef50_P07900'

HSP_cluster <- t(data[data$UniRef50.ID == uniref_cluster,-c(1,2)])
colnames(HSP_cluster) <- uniref_cluster
sample_names <- read.csv('../../Body/1_Raw/unicellular_expression_species_samles.csv')

HSP_cluster_ordered <- as.data.frame(HSP_cluster[order(HSP_cluster),])
colnames(HSP_cluster_ordered)<-uniref_cluster

HSP_cluster_ordered$species <- sample_names[rownames(HSP_cluster_ordered), 'species']

species <- c()
for (x in rownames(HSP_cluster_ordered))
  { for (y in sample_names$sample) 
    {if (x == y) 
      {HSP_cluster_ordered[x,'species'] <- sample_names[sample_names$sample == x,2]
    }
  }
}
#Дальше вручную, так как в таблице расшифровка не совсе стандартная
rownames(HSP_cluster_ordered[is.na(HSP_cluster_ordered$species),])
#"MMETSP0924" "MMETSP0398" "MMETSP0451" "MMETSP0132" "MMETSP0922" "MMETSP0378" "MMETSP0196"
HSP_cluster_ordered$species <- as.character(HSP_cluster_ordered$species)
HSP_cluster_ordered['MMETSP0924',]$species <- 'Perkinsus chesapeaki'
HSP_cluster_ordered['MMETSP0398',]$species <- 'Amphidinium carterae'
HSP_cluster_ordered['MMETSP0451',]$species <- 'Oxyrrhis marina'
HSP_cluster_ordered['MMETSP0132',]$species <- 'Symbiodinium kawagutii'
HSP_cluster_ordered['MMETSP0922',]$species <- 'Perkinsus marinus'
HSP_cluster_ordered['MMETSP0378',]$species <- 'Alexandrium tamarense'
HSP_cluster_ordered['MMETSP0196',]$species <- 'Alexandrium fundyense'

write.table(HSP_cluster_ordered, '../../Body/2_Derived/unicellular.euckariotes.ranged.expression.HSP.cluster.with.species.names.txt')
