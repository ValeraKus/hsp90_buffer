#Кристина

# Объединить все:
rbind() # Bind by column - по колонкам
cbind() # Bind by row
require('ggplot2')
RPKM1 <- read.csv("C:/Users/User/Desktop/NormalizedRPKM_ConstitutiveAlignedExons_Amniote1to1Orthologues.csv")
RPKM2 <- read.csv("C:/Users/User/Desktop/NormalizedRPKM_ConstitutiveAlignedExons_Primate1to1Orthologues.csv")
RPKM3 <- read.csv("C:/Users/User/Desktop/NormalizedRPKM_ConstitutiveExons_Amniote1to1Orthologues.csv")
RPKM4 <- read.csv("C:/Users/User/Desktop/NormalizedRPKM_ConstitutiveExons_Primate1to1Orthologues.csv")
str(RPKM1)
summary(RPKM1)
plot(RPKM1$gga[ENSG00000080824],RPKM1$hsa.br.M.1)
View(RPKM1)
newRPKM1 <- RPKM1[which()
ggplot(RPKM1,aes(x = hsa, col=))
 geom_histogram()
 RPKM1 = read.csv('NormalizedRPKM_ConstitutiveAlignedExons_Amniote1to1Orthologues.csv')
 newRPKM1 <- RPKM1[which(RPKM1$hsa == 'ENSG00000100364')]
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
 