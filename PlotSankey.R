#' @title PlotSankey
#' @description  

# sankeyData2 <- tcga.class.p2[, c('tcga.class.p2','iCluster cluster')]
# sankeyData2$`iCluster cluster` <- paste('iCluster', sankeyData2$`iCluster cluster`, sep = '')
# sankeyData2$Value <- 1
# sankeyData2 <- sankeyData2 %>% arrange(tcga.class.p2)


PlotSankey <- function(data, axis1, axis2, axis1.name, axis2.name) {
  require(ggalluvial)
  data[['Value']] <- 1
  data <- data %>% arrange(data[[axis1]], data[[axis2]])
  ggplot(data, aes_string(y= 'Value', axis1 = axis1, axis2 = axis2)) +
    geom_alluvium(aes_string(fill = axis1), width= 0.1)+ theme_classic()+
    geom_stratum(width = 1/6, aes_string(fill = axis1)) +
    geom_label(stat = "stratum", infer.label = TRUE) +
    scale_x_discrete(limits= c(axis1.name, axis2.name), expand = c(.05, .05)) 
}


# ggsankey
# colors <- c("#3B9AB2","#78B7C5","#EBCC2A","#E1AF00","#F21A00",
#             "#ECCBAE","#046C9A","#D69C4E","#ABDDDE","#000000")
# 
# df <- mtcars %>%
#   make_long(cyl, vs, am, gear, carb)
# # simple plot
# ggplot(df, aes(x = x, 
#                next_x = next_x, 
#                node = node, 
#                next_node = next_node,
#                fill = factor(node))) +
#   geom_sankey()+
#   scale_fill_manual(values = colors)+
#   theme_minimal()+labs(x=NULL)
# #
# ggplot(df, aes(x = x, next_x = next_x,
#                node = node, next_node = next_node,
#                fill = factor(node), label = node)) +
#   geom_sankey(flow.alpha = .6,
#               node.color = "gray30") +
#   geom_sankey_text(size = 3, color = "black") +
#   scale_fill_manual(values = colors) +
#   theme_sankey(base_size = 18) +
#   labs(x = NULL) +
#   theme(legend.position = "none",
#         plot.title = element_text(hjust = .5)) +
#   ggtitle("Car features")
order_his <- c('AS','ES','WDLPS','MLPS','DDLPS','MFS','otherFS','MPNST','SS','UPS','RMS','LMS')
order_hc <- c('HC1','HC2','HC3','HC4','HC5','HC6')
order_pc2 <- c('PC1','PC2','PC3')
orderAll <- c(order_his, order_hc, order_pc2)

col_his <- c('#B92025','#DDA0DD','#F0E68C','#A79930','#D2862B','#9BC799','#B3CF41','#A8D2F4','#BFBFD1','#DBEBF9','#DBA0A4','#E8C7D0')
col_hc <- c('#8C1D27','#253669','#B78AB8','#B16E29','#61668E','#53978A')
col_pc2 <- c('#DE8C3E','#D3605E','#6583B5')
colAll <- c(col_his, col_hc, col_pc2)

df <- clin2 %>% make_long(his_short, his_hc, ProteinCluster2)
df$node <- factor(df$node, levels = orderAll)

ggplot(df, aes(x = x, next_x = next_x,
               node = factor(node), next_node = next_node,
               fill = node, label = node)) +
      geom_sankey(flow.alpha = .6,
                  node.color = "gray30", node.fill = colAll, 
                  width =0.1, shift =1) +
      geom_sankey_text(size = 3) +
      #scale_fill_manual(values = colAll) +
      theme_sankey(base_size = 18) +
      labs(x = NULL) +
      theme(legend.position = "none",
            plot.title = element_text(hjust = .5)) +
      ggtitle("")
