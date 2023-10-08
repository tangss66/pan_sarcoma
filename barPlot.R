# 百分比柱状图

# 建立数据
mydata <- data.frame(id = 1:98,
                     ORR = c(rep("SD", 54), rep("CR", 44)),
                     PD_1 = c(rep("TC0", 8), rep("TC1", 19), rep("TC2", 27),
                              rep("TC0", 6), rep("TC1", 13), rep("TC2", 25)))

head(mydata) # 查看前6行

# 卡方检验
## 将数据整理为四格表
ka<-xtabs(~mydata$PD_1+mydata$ORR,data=mydata)
ka
# 卡方检验
chisq.test(ka)
# 百分比柱状图
library(scales)
library(ggplot2)
library(RColorBrewer)
p <- ggplot(mydata, aes(x = ORR, fill = PD_1)) +
  geom_bar(width = 0.9, position = "fill") + # 百分比柱状图
  scale_fill_brewer(palette = "Greens")  +
  scale_y_continuous(labels = percent) +
  labs(title = "mUC",
       y = "Patients Percentage")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"))
p

p1 <- p + coord_cartesian(clip = 'off', ylim = c(0,1))+ 
  theme(plot.margin = margin(0.5,0.5,1.2,0.5,'cm'))+ #自定义图片上左下右的边框宽度
  annotate( "text",
            cex=5,
            x=1.5, y=-0.15, # 根据自己的数据调节p value的位置
            label='p = 0.760', # 添加P值
            color="black")+
  annotate("rect", xmin = 0.55, xmax = 1.45, ymin = -0.1, ymax = -0.02, 
           fill = "#40a1cf")+
  annotate("rect", xmin = 1.55, xmax = 2.45, ymin = -0.1, ymax = -0.02, 
           fill = "#dd816d")
p1
## 同理话另外两个 且不显示y轴
p2 <- p1+theme(
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.title.y = element_blank()
)
p2

## 拼图
#install.packages("patchwork")
library(patchwork)
p1+ p2+ p2 + plot_layout(guides = 'collect')
ggsave('stack_barplot.pdf',width = 8,height = 6)