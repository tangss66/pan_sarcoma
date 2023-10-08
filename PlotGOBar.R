#' @title PlotGOBar
#' @description Plot bar plot of GO enrichment terms


PlotGOBar.base <- function(data, show.cat = 10, color = '#00CED1') {
  # recommended color #00CED1, #F08080
  data <- data %>% arrange(pvalue) %>% mutate(order=-log10(pvalue)) %>% head(show.cat)
  p <- ggplot(data = data) + geom_bar(aes(x = reorder(Description,order), y = order), fill = color, stat = 'identity') +
    coord_flip() + theme_classic() + labs(x='Enrichment terms',y='-log10(p-value)') +
    theme(legend.title = element_blank()) + geom_hline(yintercept = 0)
  return(p)
}

PlotGOBar.updown <- function(dataup,datadown,show.cat.up=10,show.cat.down=10,color=c('#00CED1','#F08080')){
  #Plot the top significant GO terms in both upregulation and downregulation gene data
  dataup <- dataup %>% arrange(pvalue) %>% mutate(order=-log10(pvalue)) %>% head(show.cat.up)
  datadown <- datadown %>% mutate(order=log10(pvalue)) %>% arrange(pvalue) %>% head(show.cat.down)
  data <- rbind(dataup,datadown) %>% arrange(desc(order)) %>% mutate(updown=ifelse(order>0,'up','down'))
  plotbar <- ggplot(data=data) + geom_bar(aes(x = reorder(Description,order),y = order,fill = updown), stat = 'identity') +
    coord_flip() + scale_fill_manual(values =color)+theme_classic() + labs(x='GO terms',y='-log10(p-value)') +
    theme(legend.title = element_blank()) + geom_hline(yintercept = 0)
  #change the breaks in y axis
  #plotbar=plotbar+scale_y_continuous(breaks = c(),labels = c())
  return(plotbar)
}



# ggplot(data = clin, aes(x = ProteinCluster2)) + geom_bar(aes(fill = metastasis), position = "fill")