#' @title PlotBubble
#' @description  plot bubble plot of enrichment 

PlotBubble <- function() {
  goc1 <- c('GO:0031424','GO:0044255')
  goc2 <- c('GO:0042254','GO:0032259')
  goc3 <- c('GO:0002274','GO:0042119')
  # goc4 <- c('GO:0002443','GO:0002263')
  goc5 <- c('GO:0003012','GO:0061061')
  goc6 <- c('GO:0030334','GO:0048514','GO:0001525')
  
  gocAll <- c(goc1, goc2, goc3, goc5,goc6)
  
  data <- engo_merge[engo_merge$ID %in% gocAll,]
  
  data$Description <- factor(data$Description, levels = c(
    'cellular lipid metabolic process',
    'keratinization',
    'methylation',
    'ribosome biogenesis',
    'neutrophil activation',
    'myeloid leukocyte activation',
    'regulation of cell migration',
    'muscle system process',
    'muscle structure development',
    'angiogenesis',
    'blood vessel morphogenesis'
  ))
  
  ggplot(data = data, aes_string(x = 'hc', y = 'Description', size = 'Count')) +
    geom_point(aes(color = -log10(pvalue)))+ theme_bw() + 
    scale_colour_gradient(low = "#FEE0D2", high = "#A50F15") +
    scale_size_continuous(range = c(2,7))
  
}