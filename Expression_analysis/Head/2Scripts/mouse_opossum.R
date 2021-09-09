,tab = read.csv("C:/Users/User/Desktop/Ортологи.csv")
mouse_opossum = tab[which(tab$organism_name == "Mus musculus" | tab$organism_name == "Monodelphis domestica", arr.ind = FALSE,  useNames = TRUE), ] 
write.table(a,"mouse_opossum.txt", sep = "\t")

