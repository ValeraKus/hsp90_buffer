rm(list=ls(all=TRUE))

RPKM <- read.csv('D:/Konstantin/Hsp90/Body/2_Derived/ENSG00000096384_RPKM_tissues.csv') #file with RPKM for 6 species in 6 tissues
Martinez <- c( 'hsa','ggo','ptr','ppy','mml' ) # The order of the organisms by effective population size according to Prado Martinez


VecOfTissues = unique(RPKM$Tissue)
VecOfSpecies = unique(RPKM$Species)
setwd('D:/Konstantin/Hsp90/Body/4_Figures')
pdf('Mortinez_N & Kaesmann_RPKM boxplots.pdf')
for (i in 1:length(VecOfTissues))
{ # i = 1
  TEMP = RPKM[RPKM$Tissue == VecOfTissues[i]]
  boxplot(TEMP[TEMP$Species == 'hsa',]$Matinez)
}
dev.off()  

#Error in `[.data.frame`(RPKM, RPKM$Tissue == VecOfTissues[i]) : 
#undefined columns selected




