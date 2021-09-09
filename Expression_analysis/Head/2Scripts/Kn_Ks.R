#Deriving table from Blast with Kn/Ks
#Opossum (Monodelphis domestica)	
#Macaque (Macaca mulatta) 
#Chicken (Gallus gallus)
#Gorilla (Gorilla gorilla gorilla)
#Mouse (Mus musculus)
#Platypus (Ornithorhynchus anatinus)
#Bonobo (Pan paniscus)
#Chimpanzee (Pan troglodytes)
#Orangutan (Pongo abelii)


rm(list=ls(all=TRUE))
df <- read.csv('D:/Konstantin/Hsp90/Body/1_Raw/orthologues-ComparaOrthologs-Homo_sapiens_Gene_Compara_Ortholog_ENSG00000096384_Placental.csv')

df = df[df$Type == '1-to-1View Gene Tree', 1:4]
df = df[,-2:-3]



L = read.csv('D:/Konstantin/Hsp90/Body/1_Raw/GenerationLenghtforMammals.csv')
        



?merge



Kn.Ks_table_for_Kaesman_species = a
setwd('D:/Konstantin/Hsp90/Body/2_Derived')
Kn.Ks_table_for_Kaesman_species = as.data.frame(Kn.Ks_table_for_Kaesman_species)

write.table(Kn.Ks_table_for_Kaesman_species, file = "Kn.Ks_table_for_Kaesman_species.csv", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")


