#' @title  LevelGene
#' @description 


LevelGene<- function(ExpData, key_tf, n =3) {
  #detect the expression level of the gene based on the median 
  for (i in key_tf){
    ExpData[[i]] <- ntile(ExpData[[i]],n) 
  }
  return(ExpData)
}


