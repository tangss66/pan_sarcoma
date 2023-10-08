#' @title  scatter_reg
#' @description 绘制散点图,并添加回归曲线

f_processData <- function(data, x, y, filter_x = NA, filter_y = NA, f.x = NA, f.y = NA) {
  if (!is.na(filter_x)) data <- data[which(data[[x]]> filter_x),]
  if (!is.na(filter_y)) data <- data[which(data[[y]]> filter_y),]
  if (typeof(f.x) %in% c('closure', 'builtin')) data[[x]] <- f.x(data[[x]])
  if (typeof(f.y) %in% c('closure', 'builtin')) data[[y]] <- f.y(data[[y]])
  return(data)
}

f_scatter_add_coef  <- function(gg.obj, method = 'pearson') {
  library(ggplot2)
  x <- as.character(gg.obj$mapping$x[[2]])
  y <- as.character(gg.obj$mapping$y[[2]])
  data <- gg.obj$data
  corResult <- cor.test(data[,x], data[,y], method = method)
  
  x_position <- fivenum(data[,x])[4]
  y_position <- fivenum(data[,y])[5]
  ano_text <- paste0('R = ', round(corResult$estimate, 2), '\n p = ', signif(corResult$p.value,2))
  gg.obj + annotate('text',x = x_position, y = y_position,label = ano_text)
}

f_corplot_one_step <- function(data,x, y, filter_x = NA, filter_y = NA, f.x = NA, f.y = NA, method = 'pearson',...) {
  library(ggplot2)
  library(ggpubr)
  data <- f_processData(data,x, y, filter_x = filter_x, filter_y = filter_y, f.x = f.x, f.y = f.y)
  p <- ggscatter(data, x, y,...)
  f_scatter_add_coef(p, method = method)
}

