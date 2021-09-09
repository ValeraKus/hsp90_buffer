rm(list=ls(all=TRUE))

data <- read.csv('../../Body/1_Raw/unicelullar_expression.csv')
HSP90 <- as.vector(data$Human.gene.symbol[grepl('HSP90',data$Human.gene.symbol)])
HSP90_UniRef50 <- as.vector(data[grepl('HSP90',data$Human.gene.symbol),]$UniRef50.ID)
print(data.frame(HSP90, HSP90_UniRef50))

#outcome:
#HSP90  HSP90_UniRef50
#HSP90AA1 UniRef50_P07900
#HSP90B1 UniRef50_P14625
