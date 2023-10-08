#' @title ChangeValues 
#' @description change all values fitting the pattern.
#' @import stringr
 


#Function==========
###-----------判断输入数据为matrix/data.frame还是vector
##-----------判断是部分匹配(仅对字符串有效)还是全匹配
##-----------全匹配
#----------判断元素是否为NA
##-----------部分匹配
#----------判断是否全替代

#循环pattern中的元素
#返回需要替换值的序号
#替换值
ChangeValues <- function(data, pattern, replacement, full.match = T, full.replace = T) {
  #update 2022.05.17
  
  #function: replace the pattern in the data by replacement
  #data : data.frame/matrix/vector
  #pattern, replacement : can be NA, string or vector, factor is invalid. 
  #the length of pattern and replacement should be the same.
  require(stringr)
  if (is.data.frame(data) | is.matrix(data)) { # detect if the data is matrix or data.frame
    if (full.match == T) { # full match  
      for (i in seq_along(pattern)){ #cycle
        for (j in seq_along(data)) {
          if (! is.na(pattern[i])){ # not null
            index <- which(data[, j] == pattern[i])
          } else index <- which(is.na(data[, j])) #null
          data[index, j] <- replacement[i]
        }
      }
    } 
    if (full.match == F) { # part match 
      for (i in seq_along(pattern)) { #cycle
        for (j in seq_along(data)) {
          index <- sapply(data[,j], function(x) str_detect(x, pattern = pattern[i]), simplify = TRUE) # require(stringr)
          if (full.replace == F) { # part replace
            data[index,j] <- sapply(data[index,j],
                                    function(x) str_replace(x, pattern = pattern[i], replacement = replacement[i]),
                                    simplify = T)
          } else data[index,j] <- replacement[i] # full replace
        }
      }
    }
    return(data)
  } else if (is.vector(data)) { #detect if the data is vector
    if (full.match == T) {  # full match
      for (i in seq_along(pattern)){ #cycle
        if (! is.na(pattern[i])){ # not null
          index <- which(data == pattern[i])
        } else index <- which(is.na(data)) #null
        data[index] <- replacement[i] 
      }
    } else { # part match
      for (i in seq_along(pattern)) { #cycle
        index <- sapply(data, function(x) str_detect(x, pattern = pattern[i]), simplify = TRUE) # require(stringr)
        if (full.replace == F) { # part replace
          data[index] <- sapply(data[index],
                                  function(x) str_replace(x, pattern = pattern[i], replacement = replacement[i]),
                                  simplify = T)
        } else data[index] <- replacement[i] # full replace
      }
    }
    return(data)
  } else {
    print("the type of data shoud be matrix, data.frame or vector")
    return(NULL)
  }
}


# #Examples==========
# score <- rnorm(30,80,20)
# score <- sapply(score, function(x) round(x,0), simplify = T)
# data <- data.frame('ID' = rep(c('A','B','C'),10), 'Num' = rep(c(1:9,NA),each = 3), 'Score' = score)
# head(data)
# #data frame single replacement
# demo1 <- ChangeValues(data, 'A', 'a')
# head(demo1)
# table(data$ID)
# table(demo1$ID)
# 
# #data frame multiple replacements
# demo2 <- ChangeValues(data, c('A','B'), c('a', 'bs'))
# demo3 <- ChangeValues(data, c('A',1), c('a', 'no data'))
# table(demo2$ID)
# head(demo3)
# 
# 
# #vector multiple replacements
# demo4 <- ChangeValues(score, c(57,89,76), c('no pass', 'excellet','good'))
# score
# demo4
# 
# #part replace
# data3 <- c('Hello World', 'Hello R', 'R studio', 'Internet')
# ChangeValues(data3, 'Hello', 'Hi', full.match = F, full.replace = F)
# ChangeValues(data3, 'Hello', 'Hi', full.match = F, full.replace = T)
# 
# #Warning
# data4 <- list('A' = c(1,2,3), 'B' ='Test Data', 'C' = c('new', 'cow', 'tween'))
# is.vector(data4)
# ChangeValues(data4, 1, 5)
# ChangeValues(list(1,2,3), 1, 5)
# ChangeValues(list(c(1,2,3), c(1,4,7), c(2,1,8)), 1, 5)

