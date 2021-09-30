rm(list=ls(all=TRUE))

hsp <- read.table("../../Body/2_Derived/hsp.kn.kn.species.orders.txt", header = T, sep = ',')

pdf("../../Body/4_Figures/hsp.kn.kn.species.orders.pdf")

library('ggplot2')
theme_set(theme_classic())
ggplot(data = hsp, aes(x = log10(Generation_Length), y = dN.dS_ENSG00000096384))+
  geom_point(aes(color = Order))+
  geom_smooth(method = lm)+
  labs(y = "Kn/Ks")+
  ggtitle('HSP90AB1')

dev.off()
