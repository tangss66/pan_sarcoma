#' @title  scale axis
#' @update 20221110
#' @description 对ggplot的x轴或y轴进行缩放,可能只能用于缩放没有点的区域
#' @references https://zhuanlan.zhihu.com/p/358781655

library(scales)

squash_axis <- function(from, to, factor) { 
  # Args:
  #   from: left end of the axis
  #   to: right end of the axis
  #   factor: the compression factor of the range [from, to]
  
  trans <- function(x) {    
    # get indices for the relevant regions
    isq <- x > from & x < to
    ito <- x >= to
    
    # apply transformation
    x[isq] <- from + (x[isq] - from)/factor
    x[ito] <- from + (to - from)/factor + (x[ito] - to)
    
    return(x)
  }
  
  inv <- function(x) {
    # get indices for the relevant regions
    isq <- x > from & x < from + (to - from)/factor
    ito <- x >= from + (to - from)/factor
    
    # apply transformation
    x[isq] <- from + (x[isq] - from) * factor
    x[ito] <- to + (x[ito] - (from + (to - from)/factor))
    
    return(x)
  }
  
  # return the transformation
  return(trans_new("squash_axis", trans, inv))
}


# shiyanhe <- data.frame(group=rep(c('A', 'B', 'C', 'D'), each = 10), 
#                        value=c(rnorm(10), rnorm(10)+100))
# # 把 5 到 95 范围的 y 轴压缩 10倍
# ggplot(shiyanhe, aes(x=group, y = value))+
#   geom_point()+
#   coord_trans(y = squash_axis(5, 95, 10))
