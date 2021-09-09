rm(list = ls(all = T))

amebas <- read.table('../../Body/2_Derived/unicellular.euckariotes.ranged.expression.HSP.cluster.with.species.names.txt')
five_percent <- round(dim(amebas)[1] * 0.05)

high_exp <- as.vector(amebas[(dim(amebas)[1]-five_percent):dim(amebas)[1],'species'])
low_exp <- as.vector(amebas[1:five_percent,'species'])


write(high_exp, '../../Body/2_Derived/unicellular_eucariotes_with_highest_expression_of_hsp.txt')

write(low_exp, '../../Body/2_Derived/unicellular_eucariotes_with_lowest_expression_of_hsp.txt')

h <- read.csv('../../Body/2_Derived/unicellular_eucariotes_with_highest_expression_of_hsp.txt')
l <- read.csv('../../Body/2_Derived/unicellular_eucariotes_with_lowest_expression_of_hsp.txt')

pdf('../../Body/4_Figures/unicellular.eukariotes.supergroup.distribution.pdf')
par(mfcol=c(1,1))
barplot(table(h[,'Infrakingdom']), main = 'Unicellular marine eucariotes with the highest hsp expression level', col='cadetblue2')
barplot(table(l[,'Infrakingdom']),main = 'Unicellular marine eucariotes with the lowest hsp expression level', col = 'cadetblue2')
dev.off()



library(ggplot2)
library(gridExtra)
h_class <- as.data.frame(table(h$Class))
l_class <- as.data.frame(table(l$Class))
h_infraphylum <- as.data.frame(table(h$Infraphylum))
l_infraphylum <- as.data.frame(table(l$Infraphylum))
h_phylum <- as.data.frame(table(h$Phylum))
l_phylum <- as.data.frame(table(l$Phylum))
h_order <- as.data.frame(table(h$Order))
l_order <- as.data.frame(table(l$Order))

pdf('../../Body/4_Figures/unicellular.eukariotes.taxonomic.distribution.pdf')

h_1 <- ggplot(data=h_phylum, aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", color = 'steelblue', fill = 'steelblue')+
  theme_light()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1), plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5))+
  labs(x = 'Phylum', y = 'Frequency')+
  ggtitle('Unicellular marine eucariotes with the highest hsp expression level')

l_1 <- ggplot(data=l_phylum, aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", color = 'steelblue', fill = 'steelblue')+
  theme_light()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1), plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5))+
  labs(x = 'Phylum', y = 'Frequency')+
  ggtitle('Unicellular marine eucariotes with the lowest hsp expression level')

h_2 <- ggplot(data=h_infraphylum, aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", color = 'steelblue', fill = 'steelblue')+
  theme_light()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1), plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5))+
  labs(x = 'Infraphylum', y = 'Frequency')+
  ggtitle('Unicellular marine eucariotes with the highest hsp expression level')

l_2 <- ggplot(data=l_infraphylum, aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", color = 'steelblue', fill = 'steelblue')+
  theme_light()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1), plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5))+
  labs(x = 'Infraphylum', y = 'Frequency')+
  ggtitle('Unicellular marine eucariotes with the lowest hsp expression level')

h_3 <- ggplot(data=h_class, aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", color = 'steelblue', fill = 'steelblue')+
  theme_light()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1), plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5))+
  labs(x = 'Class', y = 'Frequency')+
  ggtitle('Unicellular marine eucariotes with the highest hsp expression level')

l_3 <- ggplot(data=l_class, aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", color = 'steelblue', fill = 'steelblue')+
  theme_light()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1), plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5))+
  labs(x = 'Class', y = 'Frequency')+
  ggtitle('Unicellular marine eucariotes with the lowest hsp expression level')

h_4 <- ggplot(data=h_order, aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", color = 'steelblue', fill = 'steelblue')+
  theme_light()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1), plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5))+
  labs(x = 'Order', y = 'Frequency')+
  ggtitle('Unicellular marine eucariotes with the highest hsp expression level')

l_4 <- ggplot(data=l_order, aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", color = 'steelblue', fill = 'steelblue')+
  theme_light()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1), plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5))+
  labs(x = 'Order', y = 'Frequency')+
  ggtitle('Unicellular marine eucariotes with the lowest hsp expression level')

grid.arrange(h_1, l_1, nrow=2)
grid.arrange(h_2, l_2, nrow=2)
grid.arrange(h_3, l_3, nrow=2)
grid.arrange(h_4, l_4, nrow=2)


dev.off()
