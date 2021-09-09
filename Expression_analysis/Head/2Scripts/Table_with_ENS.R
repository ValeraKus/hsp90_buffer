# Creation a small table with species,tissues and EnsemblId
rm(list=ls(all=TRUE))
setwd('/hdd/SCIENCE_PROJECTS_BODY/GOOD_DATASETS/Gencode25AndGeneAnnotation/RawData/')
G = read.table("C:/Users/User/Desktop/gencode.v27.annotation.gtf", sep = '\t',header=F)

# derive EnsemblId:
G$EnsemblId = gsub("(.*)gene_id ",'',G$V9) 
G$EnsemblId = gsub(";(.*)",'',G$EnsemblId)#??????????? ?????? ? ?????????? ????? ?????
G$EnsemblId = gsub("\\.(.*)",'',G$EnsemblId)#????????????? ??? ??????????
length(unique(G$EnsemblId))  # 57 992 ? ???? 58243

GENE = G[G$V3 == 'gene',]
GENE$GeneLength = GENE$V5 - GENE$V4
GENE$GeneType = gsub("(.*)gene_type ",'',GENE$V9); GENE$GeneType = gsub(";(.*)",'',GENE$GeneType);     
GENE$GeneName = gsub("(.*)gene_name ",'',GENE$V9); GENE$GeneName = gsub(";(.*)",'',GENE$GeneName);     
GENE = data.frame(GENE$EnsemblId, GENE$V1,GENE$V4,GENE$V5,GENE$GeneLength,GENE$GeneType,GENE$GeneName); names(GENE)=c('EnsemblId','GeneChr','GeneStart','GeneEnd','GeneLength','GeneType','GeneName')
GENE$GeneChr = gsub("chr",'',GENE$GeneChr);
# derive NumberOfTranscripts, AverageNumberOfExons and AverageExonLength:
T = G[G$V3 == 'transcript',]
T$TranscriptNumber = 1
TAgg = aggregate(T$TranscriptNumber, by = list(T$EnsemblId), FUN = sum); names(TAgg) = c('EnsemblId','NumberOfTranscripts')

E = G[G$V3 == 'exon',]
E$ExonNumber = 1
E$AverageExonLength = E$V5 - E$V4
EAgg = aggregate(list(E$ExonNumber,E$AverageExonLength), by = list(E$EnsemblId), FUN = sum); names(EAgg) = c('EnsemblId','NumberOfExons','AverageExonLength')
EAgg$AverageExonLength = EAgg$AverageExonLength/EAgg$NumberOfExons

TE = merge(TAgg,EAgg, by = 'EnsemblId')

ALL = merge(GENE,TE, all.x = TRUE, by = 'EnsemblId')
ALL$AverageNumberOfExonsPerTranscript = ALL$NumberOfExons/ALL$NumberOfTranscripts
ALL= ALL[-c(9)]

write.table(ALL,'gencode.v27.annotation.gtf.Genes', sep = '\t', quote = FALSE, row.names = FALSE)
OnlyHSP = (ALL[grep("HSP", ALL$GeneName),])
ALL$GeneName

#?????

temp_data =RPKM4[RPKM4$hsa == 'ENSG00000096384',]
temp_data1 = c('hsa.br.M.1', 'hsa.br.M.2', 'hsa.br.M.3', 'hsa.br.M.4', 'hsa.br.M.5')

ggplot(RPKM4,aes(x = temp_data, col=temp_data1))
library(ggplot2)
temp_data = RPKM4[RPKM4$hsa == 'ENSG00000096384', c('hsa.br.M.2', 'hsa.br.M.2', 'hsa.br.M.3', 'hsa.br.M.4', 'hsa.br.M.5')]
newRPKM4 = t(temp_data)
qplot(row.names(newRPKM4), as.numeric(newRPKM4[,1]))

head(RPKM4)
hsa_br = c('hsa.br.M.1', 'hsa.br.M.2', 'hsa.br.M.3', 'hsa.br.M.4', 'hsa.br.M.5', 'hsa.br.F.1')
brain = RPKM4[RPKM4$hsa == 'ENSG00000096384',hsa_br] 
qplot(row.names((brain), as.numeric(brain[,1])))

a = t(temp_data)
a1= cbind
?cbind
cbind(as.data.frame(a))
class(a1)
grep(as.data.frame(a1, x, ignore.case = FALSE, perl = FALSE, value = FALSE, fixed = FALSE, useBytes = FALSE, invert = FALSE))

write.table(as.data.frame(temp_data))
?strsplit
b = strsplit(a, 'hsa' ,'br')

c = strsplit(a, 'ptr')
d = strsplit(a, 'ggo')
e = strsplit(a, 'ppy')
f = strsplit(a, 'mml')

b = as.data.frame(b)
c = as.data.frame(c)
d = as.data.frame(d)
e = as.data.frame(e)
f = as.data.frame(f)
