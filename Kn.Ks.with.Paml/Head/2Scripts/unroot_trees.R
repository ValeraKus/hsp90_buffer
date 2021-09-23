rm(list=ls(all=TRUE))
library(ggplot2)
library(ggtree)
library(ape)

trees <- list.files('../../../../Raw/clients/trees/')


for (tree in trees){
  file_tree <- read.tree(paste('../../../../Raw/clients/trees/', tree, sep = ''))
  unrooted_tree <-  unroot(file_tree)
  if (!is.rooted(unrooted_tree)){
    write.tree(unrooted_tree, paste("../../../../Raw/clients/trees_unrooted/", tree, sep = ''))
  }
  else {
    print(paste(tree, ' not unrooted!'))
    write.tree(file_tree, paste("../../../../Raw/clients/trees_unrooted/", tree, sep = ''))
  }
}



trees <- list.files('../../../../Raw/nonclients/trees/')


for (tree in trees){
  file_tree <- read.tree(paste('../../../../Raw/nonclients/trees/', tree, sep = ''))
  unrooted_tree <-  unroot(file_tree)
  if (!is.rooted(unrooted_tree)){
    write.tree(unrooted_tree, paste("../../../../Raw/nonclients/trees_unrooted/", tree, sep = ''))
  }
  else {
    print(paste(tree, ' not unrooted!'))
    write.tree(file_tree, paste("../../../../Raw/nonclients/trees_unrooted/", tree, sep = ''))
  }
}

