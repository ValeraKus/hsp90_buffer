
#скрипт который делает боксплотики

rm(list=ls(all=TRUE))
hsp <- read.csv("/home/anastasia/Desktop/Hsp90/Body/2_Derived/ENSG00000096384.csv", header = TRUE)
head(hsp)

table(hsp$Species)
table(hsp$Tissue)

VecOfTissues = unique(hsp$Tissue); length(VecOfTissues)
setwd('/home/anastasia/Desktop/Hsp90/Body/4_Figures')
pdf('HspPrimatesVsMouse.pdf')
par(mfrow=c(2,3))
for (i in 1:length(VecOfTissues))
{ # i = 1
  TEMP = hsp[hsp$Tissue == VecOfTissues[i],]
  boxplot(TEMP[TEMP$Species!= 'mml',]$RPKM, TEMP[TEMP$Species == 'mml',]$RPKM, names = c('primates','macaque'), outline = FALSE, main = VecOfTissues[i])
}
dev.off()  

