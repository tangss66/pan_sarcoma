#' @tiltle PlotCor
#' @description plot the point plot to show the correlation of x & y

PlotCor <- function(data, x, y, color.point = '#FF63487F', color.line = 'red') {
  ggplot(data = data, aes_string(x = x ,y = y)) + geom_point(size=3,color = color.point)+
  geom_smooth(method = 'lm', se = F, color = color.line, size = 0.5) + theme_classic() + 
  ggpubr::stat_cor(method = 'spearman',  size = 6)
}

