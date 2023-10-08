#' @title CalculateGroupMean

CalculateGroupMean <- function(data, group){
  cats <- unique(group)
  data.mean <- data.frame('protein' = rownames(data))
  for (i in cats){
    colIndex <- which(group == i)
    mean.sub <- apply(data[,colIndex], 1, function(x) mean(x,na.rm = T))
    data.mean <- cbind(data.mean, mean.sub)
  }
  colnames(data.mean) <- c('Protein', cats)
  return(data.mean)
}