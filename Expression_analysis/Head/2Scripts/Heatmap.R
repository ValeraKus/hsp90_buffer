
rbind() # Bind by column
cbind() # Bind by row
require('ggplot2')



RPKM_Aligned_Amniote = read.table('/home/anastasia/Desktop/Hsp90/Body/1_Raw/RPKM_1to1Orthologs/NormalizedRPKM_ConstitutiveAlignedExons_Amniote1to1Orthologues.txt', sep = '\t', header=F)
RPKM_Aligned_Primate = read.table('/home/anastasia/Desktop/Hsp90/Body/1_Raw/RPKM_1to1Orthologs/NormalizedRPKM_ConstitutiveAlignedExons_Primate1to1Orthologues.txt', sep = '\t', header=F)
RPKM_Amniote= read.table('/home/anastasia/Desktop/Hsp90/Body/1_Raw/RPKM_1to1Orthologs/NormalizedRPKM_ConstitutiveExons_Amniote1to1Orthologues.txt', sep = '\t', header=F)
RPKM_Primate= read.table('/home/anastasia/Desktop/Hsp90/Body/1_Raw/RPKM_1to1Orthologs/NormalizedRPKM_ConstitutiveExons_Primate1to1Orthologues.txt', sep = '\t', header=F)

plot(RPKM_Aligned_Primate$gga[ENSG00000096384],RPKM_Aligned_Primate$hsa.br.M.1)
View(RPKM_Aligned_Primate)
newRPKM <-RPKM_Aligned_Primate[which()ggplot(RPKM_Aligned_Primate,aes(x = hsa, col=))
 geom_histogram()
 newRPKM1 <- RPKM_Aligned_Primate[which(RPKM_Aligned_Primate$hsa == 'ENSG00000096384')]
 newRPKM1 = t(newRPKM1)
 
 x = rownames(newRPKM1)[10:141]
 require(ggplot2)
 
 qplot(x = rownames(newRPKM1), as.numeric(newRPKM1[,1]))
 
 ggplot(newRPKM1[10:141,]), aes(newRPKM1[]newRPKM1$1[10:141,]))
 write.table(OnlyHSP, 'filename.txt', sep = '\t')
 names(OnlyHSP)
 getwd()
 
 ####
 
 table <- read.csv("C:/Users/User/Desktop/ENSG00000096384.csv")
 map = as.data.frame(table[c(1,2,5)])
x = c()
 for (i in unique(map$Species)){
   for (j in unique(map$Tissue)){
     x = c(x, i, j, mean(map[which(map$Species == i & map$Tissue == j),][,3]))
   }
 }

y = data.frame(x[seq(1, length(x), 3)], x[seq(2, length(x), 3)], x[seq(3, length(x), 3)]) 
colnames(y) = c('Species', 'Tissue', 'RPKM')

 mean(map[which(map$Species == "hsa" & map$Tissue == "br"),][,3])
 library(ggplot2)

  #means
ggplot(y, aes(Species, Tissue)) + geom_tile(aes(fill=as.numeric(RPKM)), colour = "white") +
   scale_fill_gradient(low = "white", high = 'steelblue')

ggplot(map, aes(Species, Tissue)) + geom_tile(aes(fill=as.numeric(RPKM)), colour = "white") +
  scale_fill_gradient(low = "white", high = 'steelblue')


 
 p =  ggplot(as.character(map$RPKM), aes(map$Species, map$Tissue))
 gencode.v27.annotation.gtf.Genes = read.csv('/home/anastasia/Desktop/Hsp90/Body/1_Raw/gencode.v27.annotation.gtf.Genes')
 