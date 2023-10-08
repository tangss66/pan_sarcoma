#' @title  MatchProPho
#' @description 


MatchProPho <- function(pro,pho){
  library(stringr)
  phoGene <- sapply(pho, function(x) str_split(x, pattern = '/', simplify = T)[1], simplify = T)
  index <- which(phoGene %in% pro)
  selectPho <- pho[index]
  return(selectPho)
}

TransPhoPro <- function(pho){
  library(stringr)
  phoGene <- sapply(pho, function(x) str_split(x, pattern = '/', simplify = T)[1], simplify = T)
  return(phoGene)
}
