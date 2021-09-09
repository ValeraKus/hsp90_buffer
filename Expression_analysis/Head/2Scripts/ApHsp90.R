########################################################################
############ 28 March 2018: derive Kaesmann's expression levels for orthologous clusters from OrthoDB
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
Orthologs <- read.table('/hdd/SCIENCE_PROJECTS_BODY/Hsp/1_RAW/FromOrthoDB/OrthoDB.Mammals.EOG090A0320.txt', sep = '\t', header = TRUE)
VecOfSpecies = unique(Orthologs$organism_name); VecOfSpecies
KaesmannSpecies = c('Macaca mulatta', 'Pan paniscus', 'Pan troglodytes', 'Pongo abelii', 'Homo sapiens', 'Papio anubis', 'Mus musculus', 'Monodelphis domestica', 'Ornithorhynchus anatinus')
Orthologs = as.data.frame(Orthologs)
str(Orthologs)
Orthologs = Orthologs[Orthologs$organism_name %in% KaesmannSpecies,]
Orthologs$pub_gene_id = gsub("\\;(.*)",'',Orthologs$pub_gene_id)  # ENSMUSG00000023944;Hsp90ab1 => ENSMUSG00000023944
write.table(Orthologs,'/hdd/SCIENCE_PROJECTS_BODY/Hsp/1_RAW/FromOrthoDB/OrthoDB.Mammals.EOG090A0320.Kaesmann.txt')

# pan paniscus - there are no Ensembl ID!!! keasman derived them from P trogl 

#### 2 add by hand EnsemblID to "OrthoDB.Mammals.EOG090A0320.Kaesmann.txt" and save it as "OrthoDB.Mammals.EOG090A0320.Kaesmann.Edited.txt"

#### 3 For each gene from "OrthoDB.Mammals.EOG090A0320.Kaesmann.Edited.txt" look for expression level from Kaesmenn Supplementry data

rm(list=ls(all=TRUE))
Orthologs <- read.table('/hdd/SCIENCE_PROJECTS_BODY/Hsp/1_RAW/FromOrthoDB/OrthoDB.Mammals.EOG090A0320.Kaesmann.Edited.txt', sep = '\t', header = TRUE)

KaesmannMusMus<-read.table('/hdd/SCIENCE_PROJECTS_BODY/Hsp/1_RAW/Brawand_SM/Supplementary_Data2/Mouse_Ensembl57_TopHat_UniqueReads.txt', sep = '\t', header = TRUE)
KaesmannMusMus = KaesmannMouse[KaesmannMouse$GeneID %in% Orthologs$pub_gene_id,]

KaesmannMonDom<-read.table('/hdd/SCIENCE_PROJECTS_BODY/Hsp/1_RAW/Brawand_SM/Supplementary_Data2/Opossum_Ensembl57_TopHat_UniqueReads.txt', sep = '\t', header = TRUE)
KaesmannMonDom = KaesmannMonDom[KaesmannMonDom$GeneID %in% Orthologs$pub_gene_id,]

KaesmannOrnAna<-read.table('/hdd/SCIENCE_PROJECTS_BODY/Hsp/1_RAW/Brawand_SM/Supplementary_Data2/Platypus_Ensembl57_TopHat_UniqueReads.txt', sep = '\t', header = TRUE)
KaesmannOrnAna = KaesmannOrnAna[KaesmannOrnAna$GeneID %in% Orthologs$pub_gene_id,]







########################################################################
############ 
########################################################################

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
