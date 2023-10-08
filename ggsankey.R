library(tidyverse)
library(ggsankey)
data(diamonds)
View(diamonds)

# 处理数据
df <- diamonds %>%
  make_long(cut, color, clarity) #需要绘制桑基图的节点


psankey <- ggplot(df, aes(x = x, next_x = next_x,
               node = node, next_node = next_node,
               fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = .6,
              node.color = "gray30") +
  geom_sankey_text(size = 3, color = "black") +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5)) +
  ggtitle("sankey of all clusters")

psankey
# 自定义颜色 
col_cut <- c('#8C1D27','#253669','#B78AB8','#B16E29','#61668E')
names(col_cut) <- c('Very Good', 'Premium','Ideal','Good', 'Fair')

col_color <- c('#B92025','#DDA0DD','#F0E68C','#A79930','#D2862B','#9BC799','#B3CF41')
names(col_color) <- c('J','I','H','G','F','E','D')

col_clarity <- c('#A8D2F4','#BFBFD1','#DBEBF9','#DBA0A4','#E8C7D0','#377EB7','#E31A1C','#FFFF33')
names(col_clarity) <- c('I1','SI2','SI1','VS2','VS1','VVS2','VVS1','IF')

col = c(col_cut, col_color, col_clarity)


psankey + scale_fill_manual(values = col)


