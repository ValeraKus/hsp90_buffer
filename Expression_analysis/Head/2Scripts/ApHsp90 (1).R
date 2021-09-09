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
############ boxplots
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


########################################################################
############ 16 Aug 2018, extract Kn/Ks from ensembl, merge with GT and correlate Kn/Ks vs GT
########################################################################

rm(list=ls(all=TRUE))
# read Kn/Ks and GT data
KnKs = read.csv('/media/konstantinpopadin/ac45df81-e084-4d30-9653-5c57cc9b58fd/konstantinpopadin/SCIENCE_PROJECTS_BODY/Hsp/1_RAW/ComparaOrthologs/orthologues-ComparaOrthologs-Homo_sapiens_Gene_Compara_Ortholog_ENSG00000096384.csv')
GT = read.table('/media/konstantinpopadin/ac45df81-e084-4d30-9653-5c57cc9b58fd/konstantinpopadin/SCIENCE_PROJECTS_BODY/Hsp/1_RAW/GenerationLenghtforAllMammals/GenerationLenghtforMammals.xlsx.txt', sep = '\t', header = TRUE)

### keep only 1 to 1 orthologs
nrow(KnKs)
table(KnKs$Type)
KnKs = KnKs[KnKs$Type == '1-to-1View Gene Tree',]
nrow(KnKs)

### filter out raws with NA dn.ds
KnKs = KnKs[KnKs$dN.dS != 'n/a',]
str(KnKs)
KnKs$dN.dS = as.numeric(as.character(KnKs$dN.dS))
nrow(KnKs)

### edit name
KnKs$Scientific_name = gsub("(.*)\\(",'',KnKs$Species)
KnKs$Scientific_name = gsub("\\)",'',KnKs$Scientific_name)
KnKs$Scientific_name = gsub('Canis lupus familiaris','Canis lupus',KnKs$Scientific_name)
KnKs$Scientific_name = gsub('Colobus angolensis palliatus','Colobus angolensis',KnKs$Scientific_name)
KnKs$Scientific_name = gsub('Gorilla gorilla gorilla','Gorilla gorilla',KnKs$Scientific_name)
KnKs$Scientific_name = gsub('Mustela putorius furo','Mustela putorius',KnKs$Scientific_name)
KnKs$Scientific_name = gsub('Peromyscus maniculatus bairdii','Peromyscus maniculatus',KnKs$Scientific_name)
KnKs$Scientific_name = gsub('Panthera tigris altaica','Panthera tigris',KnKs$Scientific_name)

### take subset of columns
KnKs = KnKs[,grepl("Scientific_name|dN.dS", names(KnKs))]
KnKs = aggregate(KnKs$dN.dS, by = list(KnKs$Scientific_name), FUN = mean)
names(KnKs)=c('Scientific_name','dN.dS')

### keep subset of columns from GT
GT = GT[,grepl("Scientific_name|GenerationLength_d", names(GT)),]

### merge by Scientific_name
ALL = merge(KnKs,GT, by = 'Scientific_name')
cor.test(ALL$dN.dS,ALL$GenerationLength_d, method = 'spearman')
pdf('/media/konstantinpopadin/ac45df81-e084-4d30-9653-5c57cc9b58fd/konstantinpopadin/SCIENCE_PROJECTS_BODY/Hsp/4_FIGURES/KnKsVsGtMammalsHsp.pdf')
plot(log2(ALL$GenerationLength_d),log2(ALL$dN.dS), cex = 2, pch = 16, col = 'black')
A<-lm(log2(ALL$dN.dS)~log2(ALL$GenerationLength_d))
summary(A) # intercept doesn't differ from zero => make it from the zero
abline(A, col = 'red', lwd = 3)
dev.off()


## my regression from Popadin (MBE) is very similar (at least not steeper): Kn/Ks = 0.094 + 1.18*10^(-5)*GT
A<-lm(ALL$dN.dS~ALL$GenerationLength_d)
summary(A) # intercept doesn't differ from zero => make it from the zero
B<-lm(ALL$dN.dS~0+ALL$GenerationLength_d)
summary(B) 
#Call:
#  lm(formula = ALL$dN.dS ~ 0 + ALL$GenerationLength_d)
#
#Residuals:
#  Min       1Q   Median       3Q      Max 
#-0.19021 -0.07294 -0.03678 -0.01115  0.49727 
#
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#ALL$GenerationLength_d 2.158e-05  5.452e-06   3.958  0.00033 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.1403 on 37 degrees of freedom
#Multiple R-squared:  0.2975,	Adjusted R-squared:  0.2785 
#F-statistic: 15.67 on 1 and 37 DF,  p-value: 0.0003297
