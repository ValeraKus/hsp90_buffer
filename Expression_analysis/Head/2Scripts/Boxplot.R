#скрипт который делает боксплотики

rm(list=ls(all=TRUE))
setwd('/hdd/SCIENCE_PROJECTS_BODY/Hsp/1_RAW/Brawand_SM');
hsp <- read.table("Expr.txt", header = TRUE)
head(hsp)

table(hsp$Species)
table(hsp$Tissue)

VecOfTissues = unique(hsp$Tissue); length(VecOfTissues)
setwd('/hdd/SCIENCE_PROJECTS_BODY/Hsp/4_FIGURES/');
pdf('HspPrimatesVsMouse.pdf')
par(mfrow=c(2,3))
for (i in 1:length(VecOfTissues))
{ # i = 1
  TEMP = hsp[hsp$Tissue == VecOfTissues[i],]
  boxplot(TEMP[TEMP$Species!= 'mml',]$RPKM, TEMP[TEMP$Species == 'mml',]$RPKM, names = c('primates','mouse'), outline = FALSE, main = VecOfTissues[i])
}
dev.off()  

