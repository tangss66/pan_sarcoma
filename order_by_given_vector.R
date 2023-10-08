order_by_given_vector <- function(data, order_vector) {
  new_order <- c()
  for (i in order_vector) {
    index <- which(data == i)
    new_order <- append(new_order, index)
  }
  return(new_order)
}
  