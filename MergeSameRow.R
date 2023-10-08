#' @title  MergeSameRow
#' @description merge the rows with same name or ID

# data <- exp_phoGene_lms

MergeSameRow.numeric <- function(data, na.rm = T) {
  #the first columns of data is name or ID
  #the other columns except is numeric
  if (na.rm == T) {
    data <- apply(data, c(1,2), function(x) {
      if(is.na(x)){
      return(0)
      } else return(x)
      })
  }
  
  if (ncol(data) != 2) {
    data_name <- data[, 1]
    data2 <- data[, -1]
    data2 <- apply(data2, c(1,2), as.numeric)
    dup_num <- which(duplicated(data_name))
    dup_name <- unique(data_name[dup_num])
    for (i in dup_name) {
      row_num <- which(data_name == i)
      dup_data <- data2[row_num[-1], ]
      if (is.null(dim(dup_data))) {
        data2[row_num[1],] <- data2[row_num[1],] + dup_data
      } else {
        x <- data2[row_num[1],]
        dim(x) <- c(1, length(x))
        data2[row_num[1],] <- sweep(x, 2, apply(dup_data, 2, sum), FUN = '+')
        }
    }
    data3 <- cbind(data_name, data2)
    data3 <- data3[-dup_num,]
    return(data3)
  } else {
    dup_num <- which(duplicated(data[, 1]))
    if (length(dup_num) != 0) {
      dup_name <- data[dup_num, 1]
      dup_name <- unique(dup_name)
      for (i in dup_name) {
        row_num <- which(data[, 1] == i)
        for (j in seq(2, length(row_num))) {
          data[row_num[1], 2] = as.numeric(data[row_num[1], 2]) + as.numeric(data[row_num[j], 2])
        }
      }
      data <- data[-dup_num, ]
    }
    return(data)
  }
}



MergeSameRow.character <- function(data,sep = ';') {
  data_name <- data[, 1]
  data2 <- data[, -1]
  dup_num <- which(duplicated(data_name))
  dup_name <- unique(data_name[dup_num])
  if (ncol(data) !=2 ) {
    for (i in dup_name) {
      row_num <- which(data_name == i)
      dup_data <- data2[row_num, ]
      dup_data_sum <- apply(dup_data, 2, function(x) {paste(x,collapse = sep)})
      data2[row_num[1],] <- dup_data_sum
    }
    data3 <- cbind(data_name, data2)
  } else {
      for (i in dup_name) {
        row_num <- which(data_name == i)
        dup_data <- data2[row_num]
        dup_data_sum <- paste(dup_data,collapse = sep)
        data2[row_num[1]] <- dup_data_sum
      }
      data3 <- data.frame(data_name,data2)
      colnames(data3) <- colnames(data)
  }
  data3 <- data3[-dup_num,]
  return(data3)
}
