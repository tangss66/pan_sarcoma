##南丁格尔图
library(ggplot2)
a <- data.frame(read.csv('STE3.csv'))
pdf()
ggplot(a,aes(Symbol,Freq,colour=Subtype))+
  geom_bar(stat="identity",
           #aes(fill=Subtype),
           position="dodge", # 普通柱形图
           fill="transparent",# 填充透明度
           size=1)+
  coord_polar()+ #极坐标转换
  #scale_x_continuous(breaks = 1:5,labels = LETTERS[1:5])+
  facet_wrap(~Subtype,nrow=2) # 水平分割，分割后的图呈两行排列
dev.off()

