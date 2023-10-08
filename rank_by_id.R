
rank_by_id <- function(data, col_id, id_rank) {
  indexAll <- c()
  for (i in id_rank) {
    index <- which(data[,col_id] == i)
    indexAll <- append(indexAll, index)
  }
  data[indexAll,]
}