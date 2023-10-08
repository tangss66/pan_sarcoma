#' @title  PlotHRbar
#' @description  plot the bar plot of Hazard ratio

PlotHRbar <- function(){
  ggplot(data = data_p2_left_ha, aes(x = id, y = HR)) +
    geom_errorbar(aes(ymin = HR.95L, ymax = HR.95H),width = 0.2, position = position_dodge(0.9)) +
    geom_point(size = 4,color = '#3e848c') +
    theme_classic() + geom_hline(yintercept = 1,lty='dotdash',col="black",lwd=0.8) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_y_continuous(limits = c(0.5, 3)) 
}