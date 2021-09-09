########################################################################
############ 28 March 2018: derive Kaesmann's expression levels for orthologous clusters from OrthoDB (one-to-many)
########################################################################

### SPECIES:
# Pan paniscus = bonobo
# Pan troglodytes = chimp
# Pongo abelii = orangutan
# Papio anubis = baboon
# Monodelphis domestica = opossum
# Ornithorhynchus anatinus = platypus

#### 1 read "OrthoDB.Mammals.EOG090A0320.txt" and save only species from Kaesmann paper: "OrthoDB.Mammals.EOG090A0320.Kaesmann.txt"
rm(list=ls(all=TRUE))
Orthologs <- read.table('/home/anastasia/Desktop/Hsp90/Body/1_Raw/OrthoDB.Mammals.EOG090A0320.csv', sep = '\t', header = TRUE)
VecOfSpecies = unique(Orthologs$organism_name); VecOfSpecies
KaesmannSpecies = c('Macaca mulatta', 'Pan paniscus', 'Pan troglodytes', 'Pongo abelii', 'Gorilla gorilla gorilla', 'Homo sapiens', 'Papio anubis', 'Mus musculus', 'Monodelphis domestica', 'Ornithorhynchus anatinus')
Orthologs = as.data.frame(Orthologs)
str(Orthologs)
Orthologs = Orthologs[Orthologs$organism_name %in% KaesmannSpecies,]
Orthologs$pub_gene_id = gsub("\\;(.*)",'',Orthologs$pub_gene_id)  # ENSMUSG00000023944;Hsp90ab1 => ENSMUSG00000023944
Orthologs[Orthologs$organism_name == 'Homo sapiens' & Orthologs$pub_gene_id == 'HSP90AA1',]$pub_gene_id = 'ENSG00000080824'  # HSP90AA1 => 'ENSG00000080824' from OrthDB
Orthologs[Orthologs$organism_name == 'Pan paniscus' & Orthologs$pub_gene_id == '100983575',]$pub_gene_id = 'ENSPTRG00000006736'  # Because Kaesmann's annotation is based on Pan troglodites
Orthologs[Orthologs$organism_name == 'Pan paniscus' & Orthologs$pub_gene_id == '100973292',]$pub_gene_id = 'ENSPTRG00000018217'  #  Because Kaesmann's annotation is based on Pan troglodites
write.table(Orthologs,'/home/anastasia/Desktop/Hsp90/Body/2_Derived/OrthoDB.Mammals.EOG090A0320.Kaesmann.txt')

##### pan paniscus - there are no Ensembl ID!!! keasman derived them from P trogl  !!!!!!!!!!!!!!!!!!!!

#### 2 For each gene from "OrthoDB.Mammals.EOG090A0320.Kaesmann.Edited.txt" look for expression level from Kaesmenn Supplementry data

rm(list=ls(all=TRUE))
Orthologs <- read.table('/home/anastasia/Desktop/Hsp90/Body/2_Derived/OrthoDB.Mammals.EOG090A0320.Kaesmann.txt', header = TRUE)

KaesmannMusMus<-read.table('/home/anastasia/Desktop/Hsp90/Body/1_Raw/Kaesmann_1 to_many/Mouse_Ensembl57_TopHat_UniqueReads.txt', sep = '\t', header = TRUE)
KaesmannMusMus = KaesmannMusMus[KaesmannMusMus$GeneID %in% Orthologs$pub_gene_id,]

KaesmannMonDom<-read.table('/home/anastasia/Desktop/Hsp90/Body/1_Raw/Kaesmann_1 to_many/Opossum_Ensembl57_TopHat_UniqueReads.txt', sep = '\t', header = TRUE)
KaesmannMonDom = KaesmannMonDom[KaesmannMonDom$GeneID %in% Orthologs$pub_gene_id,]

KaesmannOrnAna<-read.table('/home/anastasia/Desktop/Hsp90/Body/1_Raw/Kaesmann_1 to_many/Platypus_Ensembl57_TopHat_UniqueReads.txt', sep = '\t', header = TRUE)
KaesmannOrnAna = KaesmannOrnAna[KaesmannOrnAna$GeneID %in% Orthologs$pub_gene_id,]

KaesmannBonobo<- read.table('/home/anastasia/Desktop/Hsp90/Body/1_Raw/Kaesmann_1 to_many/Bonobo_Ensembl57_TopHat_UniqueReads.txt', sep = '\t', header = TRUE)
KaesmannBonobo= KaesmannBonobo[KaesmannBonobo$GeneID %in% Orthologs$pub_gene_id,]

KaesmannChicken<- read.table('/home/anastasia/Desktop/Hsp90/Body/1_Raw/Kaesmann_1 to_many/Chicken_Ensembl57_TopHat_UniqueReads.txt', sep = '\t', header = TRUE)
KaesmannChicken= KaesmannChicken[KaesmannChicken$GeneID %in% Orthologs$pub_gene_id,]

KaesmannChimpanzee<- read.table('/home/anastasia/Desktop/Hsp90/Body/1_Raw/Kaesmann_1 to_many/Chimpanzee_Ensembl57_TopHat_UniqueReads.txt', sep = '\t', header = TRUE)
KaesmannChimpanzee= KaesmannChimpanzee[KaesmannChimpanzee$GeneID %in% Orthologs$pub_gene_id,]

KaesmannGorilla<- read.table('/home/anastasia/Desktop/Hsp90/Body/1_Raw/Kaesmann_1 to_many/Gorilla_Ensembl57_TopHat_UniqueReads.txt', sep = '\t', header = TRUE)
KaesmannGorilla= KaesmannGorilla[KaesmannGorilla$GeneID %in% Orthologs$pub_gene_id,]

KaesmannHuman<- read.table('/home/anastasia/Desktop/Hsp90/Body/1_Raw/Kaesmann_1 to_many/Human_Ensembl57_TopHat_UniqueReads.txt', sep = '\t', header = TRUE)
KaesmannHuman= KaesmannHuman[KaesmannHuman$GeneID %in% Orthologs$pub_gene_id,]

KaesmannMacaque<- read.table('/home/anastasia/Desktop/Hsp90/Body/1_Raw/Kaesmann_1 to_many/Macaque_Ensembl57_TopHat_UniqueReads.txt', sep = '\t', header = TRUE)
KaesmannMacaque= KaesmannMacaque[KaesmannMacaque$GeneID %in% Orthologs$pub_gene_id,]

KaesmannOrangutan<- read.table('/home/anastasia/Desktop/Hsp90/Body/1_Raw/Kaesmann_1 to_many/Orangutan_Ensembl57_TopHat_UniqueReads.txt', sep = '\t', header = TRUE)
KaesmannOrangutan= KaesmannOrangutan[KaesmannOrangutan$GeneID %in% Orthologs$pub_gene_id,]

KaesmannOrangutan<- read.table('/home/anastasia/Desktop/Hsp90/Body/1_Raw/Kaesmann_1 to_many/Orangutan_Ensembl57_TopHat_UniqueReads.txt', sep = '\t', header = TRUE)
KaesmannOrangutan= KaesmannOrangutan[KaesmannOrangutan$GeneID %in% Orthologs$pub_gene_id,]

pdf('/home/anastasia/Desktop/Hsp90/Body/4_Figures/BoxplotKaesmanOneToMany.pdf', width = 14, height = 14)
GenderVec = c('All', 'Male', 'Female')
for (Gender in GenderVec)
#Gender = 'Female'
{ # 
Final = c()
VecOfSpecies = c('KaesmannBonobo','KaesmannMonDom','KaesmannChimpanzee','KaesmannGorilla','KaesmannHuman','KaesmannMacaque','KaesmannMusMus','KaesmannOrangutan','KaesmannOrnAna')
VecOfTissues= c('Cerebellum', "Brain",'Heart','Kidney','Liver',"Testis")
for (sp in VecOfSpecies)
{ # sp = 'KaesmannBonobo'
  if (sp == 'KaesmannBonobo') {Temp = KaesmannBonobo} 
  if (sp == 'KaesmannMonDom') {Temp = KaesmannMonDom} 
  if (sp == 'KaesmannChimpanzee') {Temp = KaesmannChimpanzee} 
  if (sp == 'KaesmannGorilla') {Temp = KaesmannGorilla}
  if (sp == 'KaesmannHuman') {Temp = KaesmannHuman}
  if (sp == 'KaesmannMacaque') {Temp = KaesmannMacaque}
  if (sp == 'KaesmannMusMus') {Temp = KaesmannMusMus}
  if (sp == 'KaesmannOrangutan') {Temp = KaesmannOrangutan}
  if (sp == 'KaesmannOrnAna') {Temp = KaesmannOrnAna}
  
  if (Gender == 'Male') {Temp = Temp[grep('Male', colnames(Temp))]}
  if (Gender == 'Female') {Temp = Temp[grep('Female', colnames(Temp))]}
  
  for (tis in VecOfTissues)
  {
    x = Temp[grep(tis, colnames(Temp))] 
    if (ncol(x) > 0)
    {
      y = unlist(x)
      Line = data.frame(sp,tis,median(y))
      Final = rbind(Final,Line)
    }
  }
}

Min = min(Final$median.y.)
Max = max(Final$median.y.)
LittleSpecies = c('KaesmannMusMus','KaesmannMacaque')
BigSpecies = c('KaesmannBonobo','KaesmannChimpanzee','KaesmannGorilla','KaesmannHuman')
par(mfrow=c(3,2))
for (tis in VecOfTissues)
{ # tis = 'Brain'
  x = Final[Final$tis == tis,]
  if (ncol(x) > 0 & nrow(x) > 0)
  {
  boxplot(x[x$sp %in% LittleSpecies,]$median.y.,x[x$sp %in% BigSpecies,]$median.y., main = paste(Gender,tis,sep = ' '), names = c('LittleSpecies','BigSpecies'), ylim=c(Min, Max))
  plot(x$sp, x$median.y.,main = paste(Gender,tis,sep = ' '), ylim=c(Min, Max))
  }
}

dev.off()


### plot segments

#Data = data.frame(cbind(VecOfSpecies,c(1:9))); names(Data) = c('sp','Number')
##Data1 = merge(Final,Data)
#Temp = Data1[Data1$tis == 'Brain',]
#Temp$Number = as.numeric(Temp$Number)
#Temp = Temp[order(Temp$Number),]
#plot(Temp$sp, Temp$median.y., type = "l", col = 'red', lwd = 2)


  

#переменная для таблички для хранения
#for (i in VecOfTissues)
#{ x = Final[Final$i == i,]
#boxplot(x[x$sp %in% LittleSpecies,]$median.y.,x[x$sp %in% BigSpecies,]$median.y., main = i, names = c('LittleSpecies','BigSpecies'), ylim=c(Min, Max))
#}

#вытащить из датафрейма кусок для ткани
#посчитать среднее
#запихать это среднее в таблицу
#табличка из ортологов один к многим
#





########################################################################
############ 
########################################################################

#rm(list=ls(all=TRUE))
#RPKM_Aligned_Amniote = read.table('/home/anastasia/Desktop/Hsp90/Body/1_Raw/RPKM_1to1Orthologs/NormalizedRPKM_ConstitutiveAlignedExons_Amniote1to1Orthologues.txt', sep = '\t', header=F)
#RPKM_Aligned_Primate = read.table('/home/anastasia/Desktop/Hsp90/Body/1_Raw/RPKM_1to1Orthologs/NormalizedRPKM_ConstitutiveAlignedExons_Primate1to1Orthologues.txt', sep = '\t', header=F)
#RPKM_Amniote= read.table('/home/anastasia/Desktop/Hsp90/Body/1_Raw/RPKM_1to1Orthologs/NormalizedRPKM_ConstitutiveExons_Amniote1to1Orthologues.txt', sep = '\t', header=F)
#RPKM_Primate= read.table('/home/anastasia/Desktop/Hsp90/Body/1_Raw/RPKM_1to1Orthologs/NormalizedRPKM_ConstitutiveExons_Primate1to1Orthologues.txt', sep = '\t', header=F)

