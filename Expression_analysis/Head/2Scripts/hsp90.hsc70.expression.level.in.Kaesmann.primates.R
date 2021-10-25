rm(list = ls(all = T))

ortho <- read.csv('../../Body/1_Raw/RPKM_1to1Orthologs/NormalizedRPKM_ConstitutiveAlignedExons_Primate1to1Orthologues.txt', sep = '\t')


hsp90 <- ortho[ortho$hsa == 'ENSG00000096384',]
hsc70 <- ortho[ortho$hsa == 'ENSG00000109971',]

library(tidyr)

hsp90 <- hsp90[,-c(1,2,3,4,5)] %>%
  pivot_longer(colnames(hsp90[,-c(1,2,3,4,5)]), names_to = 'species', values_to = 'RPKM')
hsp90 <- separate(hsp90, col = 1, sep = '\\.', into = c ('species', 'tissue', 'sex', 'replicate'))

hsp90[hsp90$species == 'hsa','species'] <- 'Homo sapiens'
hsp90[hsp90$species == 'ptr','species'] <- 'Pan troglodytes'
hsp90[hsp90$species == 'ppa','species'] <- 'Pan paniscus'
hsp90[hsp90$species == 'ggo','species'] <- 'Gorilla gorilla'
hsp90[hsp90$species == 'ppy','species'] <- 'Pongo pygmaeus'
hsp90[hsp90$species == 'mml','species'] <- 'Macaca mulatta'

hsp90[hsp90$tissue == 'br', 'tissue'] <- 'Brain'
hsp90[hsp90$tissue == 'cb', 'tissue'] <- 'Cerebellum'
hsp90[hsp90$tissue == 'ht', 'tissue'] <- 'Heart'
hsp90[hsp90$tissue == 'kd', 'tissue'] <- 'Kidney'
hsp90[hsp90$tissue == 'lv', 'tissue'] <- 'Liver'
hsp90[hsp90$tissue == 'ts', 'tissue'] <- 'Testis'


hsc70 <- hsc70[,-c(1,2,3,4,5)] %>%
  pivot_longer(colnames(hsc70[,-c(1,2,3,4,5)]), names_to = 'species', values_to = 'RPKM')
hsc70 <- separate(hsc70, col = 1, sep = '\\.', into = c ('species', 'tissue', 'sex', 'replicate'))

hsc70[hsc70$species == 'hsa','species'] <- 'Homo sapiens'
hsc70[hsc70$species == 'ptr','species'] <- 'Pan troglodytes'
hsc70[hsc70$species == 'ppa','species'] <- 'Pan paniscus'
hsc70[hsc70$species == 'ggo','species'] <- 'Gorilla gorilla'
hsc70[hsc70$species == 'ppy','species'] <- 'Pongo pygmaeus'
hsc70[hsc70$species == 'mml','species'] <- 'Macaca mulatta'

hsc70[hsc70$tissue == 'br', 'tissue'] <- 'Brain'
hsc70[hsc70$tissue == 'cb', 'tissue'] <- 'Cerebellum'
hsc70[hsc70$tissue == 'ht', 'tissue'] <- 'Heart'
hsc70[hsc70$tissue == 'kd', 'tissue'] <- 'Kidney'
hsc70[hsc70$tissue == 'lv', 'tissue'] <- 'Liver'
hsc70[hsc70$tissue == 'ts', 'tissue'] <- 'Testis'

library(ggplot2)

#####hsp90
p1 <- ggplot(hsp90, aes(y = RPKM, x = reorder(species, RPKM), fill = species))+
  geom_boxplot()+ 
  facet_wrap( ~ tissue, nrow = 2)+
  theme_bw()+
  xlab('')+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.8, hjust = 0.8, size = 12), legend.position = 'None',
        text = element_text(size = 13, face = 'bold'))
  
ggsave('../../Body/4_Figures/hsp90.expression.primates.Kaesmann.pdf', p1)  
  

ggplot(hsp90[(hsp90$species == 'Macaca mulatta' | hsp90$species == 'Gorilla gorilla') & hsp90$tissue == 'Heart',], aes(y = RPKM, x = reorder(species, RPKM), fill = species))+
  geom_boxplot()+ 
  facet_wrap( ~ tissue, nrow = 2)+
  theme_bw()+
  theme(axis.text.x = element_text(angle =335, vjust = 0.6, hjust=0.2, size = 27), text = element_text(size = 30, face = 'bold'),
         axis.title.y = element_text(size = 25), legend.position = 'None')+
  xlab('')+
  ylab(paste('Expression level of HSP90', 'RPKM', sep = '\n'))

ggplot(hsp90[hsp90$species == 'Macaca mulatta' | hsp90$species == 'Homo sapiens',], aes(y = RPKM, x = reorder(species, RPKM), fill = species))+
  geom_boxplot()+ 
  facet_wrap( ~ tissue, nrow = 2)+
  theme_bw()

ggplot(hsp90[hsp90$species == 'Macaca mulatta' | hsp90$species == 'Pan paniscus',], aes(y = RPKM, x = reorder(species, RPKM), fill = species))+
  geom_boxplot()+ 
  facet_wrap( ~ tissue, nrow = 2)+
  theme_bw()

ggplot(hsp90[hsp90$species == 'Macaca mulatta' | hsp90$species == 'Pan troglodytes',], aes(y = RPKM, x = reorder(species, RPKM), fill = species))+
  geom_boxplot()+ 
  facet_wrap( ~ tissue, nrow = 2)+
  theme_bw()


ggplot(hsp90[hsp90$species == 'Macaca mulatta' | hsp90$species == 'Pongo pygmaeus',], aes(y = RPKM, x = reorder(species, RPKM), fill = species))+
  geom_boxplot()+ 
  facet_wrap( ~ tissue, nrow = 2)+
  theme_bw()



#####hsc70
p1 <- ggplot(hsc70, aes(y = RPKM, x = reorder(species, -RPKM), fill = species))+
  geom_boxplot()+ 
  facet_wrap( ~ tissue, nrow = 2)+
  theme_bw()+
  xlab('')+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.8, hjust = 0.8, size = 12), legend.position = 'None',
        text = element_text(size = 13, face = 'bold'))

ggsave('../../Body/4_Figures/hsc70.expression.primates.Kaesmann.pdf', p1)  

ggplot(hsc70[hsc70$species == 'Macaca mulatta' | hsc70$species == 'Gorilla gorilla',], aes(y = RPKM, x = reorder(species, -RPKM), fill = species))+
  geom_boxplot()+ 
  facet_wrap( ~ tissue, nrow = 2)+
  theme_bw()+
  theme(axis.text.x = element_text(angle =335, vjust = 0.6, hjust=0.2, size = 13), text = element_text(size = 15, face = 'bold'),
        axis.title.y = element_text(size = 15), legend.position = 'None')+
  xlab('')+
  ylab('Expression level of HSC70, RPKM')

ggplot(hsc70[hsc70$species == 'mml' | hsc70$species == 'hsa',], aes(y = RPKM, x = reorder(species, RPKM), fill = species))+
  geom_boxplot()+ 
  facet_wrap( ~ tissue, nrow = 2)+
  theme_bw()

ggplot(hsc70[hsc70$species == 'mml' | hsc70$species == 'ppa',], aes(y = RPKM, x = reorder(species, RPKM), fill = species))+
  geom_boxplot()+ 
  facet_wrap( ~ tissue, nrow = 2)+
  theme_bw()

ggplot(hsc70[hsc70$species == 'mml' | hsc70$species == 'ptr',], aes(y = RPKM, x = reorder(species, RPKM), fill = species))+
  geom_boxplot()+ 
  facet_wrap( ~ tissue, nrow = 2)+
  theme_bw()


ggplot(hsc70[hsc70$species == 'mml' | hsc70$species == 'ppy',], aes(y = RPKM, x = reorder(species, RPKM), fill = species))+
  geom_boxplot()+ 
  facet_wrap( ~ tissue, nrow = 2)+
  theme_bw()







