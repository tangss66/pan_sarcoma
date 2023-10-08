#' @title PlotBoxPlot
#' @importFrom ggplot2

PlotViolin <- function(data, x, y, x.lab,y.lab, compares = list(c(''),c('')), color ) {
  require(ggplot2)
  require(ggsignif)
  ggplot(data=data, aes_string(x = x ,y = y))+ geom_violin(aes_string(color = x)) +   theme_classic() + theme(legend.position = '') + labs(x = x.lab, y= y.lab) +
  geom_signif(comparisons = compares,step_increase = 0.1) +geom_boxplot(aes_string(color = x),width=0.1,position = position_identity(),fill="white") +
  stat_summary(fun="mean",geom="point",shape=23, size=2,aes_string(fill = x)) +
  theme_classic() + theme(legend.position = "") +
  scale_color_manual(values = color)+ scale_fill_manual(values = color)
}

PlotBox <- function(data, x, y, x.lab,y.lab, compares= list(c(''),c('')), color) {
  require(ggplot2)
  require(ggsignif)
  ggplot(data=data, aes_string(x = x, y = y, fill = x))+
    geom_boxplot() +   theme_classic() + theme(legend.position = '') + labs(y = y.lab, x = x.lab) +
    geom_jitter(aes_string(group = x)) + scale_fill_manual(values = color) +
    geom_signif(comparisons = compares, test = t.test,step_increase = 0.1)
  
}

PlotViolin2 <- function(data, x, y, x.lab,y.lab, compares = list(c(''),c('')), color ) {
  require(ggplot2)
  require(ggsignif)
  ggplot(data=data, aes_string(x = x ,y = y))+ geom_violin(aes_string(color = x)) +   theme_classic() + theme(legend.position = '') + labs(x = x.lab, y= y.lab) +
    geom_signif(comparisons = compares, test = t.test,step_increase = 0.1) +geom_boxplot(aes_string(color = x),width=0.1,position = position_identity(),fill="white") +
    stat_summary(fun="mean",geom="point",shape=23, size=2,aes_string(fill = x)) +
    theme_classic() + theme(legend.position = "") + geom_jitter(aes_string(group = x)) + 
    scale_color_manual(values = color)+ scale_fill_manual(values = color)
}


# PlotBox2 <- function() {
#   P1 <- ggplot(Data1,aes(x=Group,y=Value,fill=Group))+ #”fill=“设置填充颜色
#     stat_boxplot(geom = "errorbar",width=0.15,aes(color="black"))+ #由于自带的箱形图没有胡须末端没有短横线，使用误差条的方式补上
#     geom_boxplot(size=0.5,fill="white",outlier.fill="white",outlier.color="white")+ #size设置箱线图的边框线和胡须的线宽度，fill设置填充颜色，outlier.fill和outlier.color设置异常点的属性
#     geom_jitter(aes(fill=Group),width =0.2,shape = 21,size=2.5)+ #设置为向水平方向抖动的散点图，width指定了向水平方向抖动，不改变纵轴的值
#     scale_fill_manual(values = c("#E69F00", "#0072B2","#F0E442"))+  #设置填充的颜色
#     scale_color_manual(values=c("black","black","black"))+ #设置散点图的圆圈的颜色为黑色
#     ggtitle("Car Milleage Data")+ #设置总的标题
#     theme_bw()+ #背景变为白色
#     theme(legend.position="none", #不需要图例
#           axis.text.x=element_text(colour="black",family="Times",size=14), #设置x轴刻度标签的字体属性
#           axis.text.y=element_text(family="Times",size=14,face="plain"), #设置x轴刻度标签的字体属性
#           axis.title.y=element_text(family="Times",size = 14,face="plain"), #设置y轴的标题的字体属性
#           axis.title.x=elesment_text(family="Times",size = 14,face="plain"), #设置x轴的标题的字体属性
#           plot.title = element_text(family="Times",size=15,face="bold",hjust = 0.5), #设置总标题的字体属性
#           panel.grid.major = element_blank(), #不显示网格线
#           panel.grid.minor = element_blank())+
#     ylab("Miles Per Gallon")+xlab("Number of Cylinders") #设置x轴和y轴的标题
# 
# }


# p1<-ggboxplot(ROC,'outcome','a',color = 'outcome',palette = 'lancet',add = 'jitter',legend='none',ggtheme = theme_bw())+
#   stat_compare_means(comparisons = list(c('group1','group2')),method = 't.test')


f_boxplot <- function(data, x = 'ProteinCluster2', y,col_pal = col_pc, filter_y = NA, f = NA, plot = 'b',compare = list(c('PC1','PC2'),c('PC2','PC3'),c('PC1','PC3')),...) {
  library(ggpubr)
  plot <- match.arg(plot, c('b','v'))
  if (!is.na(filter_y)) data <- data[which(data[[y]]> filter_y),]
  if (typeof(f) %in% c('closure', 'builtin')) data[[y]] <- f(data[[y]])
  data[[x]] <- factor(data[[x]], levels = names(col_pc))
  
  sum_mean <- aggregate(data[[y]], by = list(data[[x]]), mean)
  colnames(sum_mean) <- c('group',paste0(y,' mean'))  
  print(sum_mean)
  
  if (plot == 'b') {
    ggboxplot(data, x, y, xlab = '', 
              fill = x, add ="jitter", add.params = list(color = 'grey', alpha = 0.8), palette = col_pal,...) +
      stat_compare_means(comparisons = compare, method = 't.test') +theme(legend.position = '' )
  } else if (plot == 'v') {
    ggviolin(data, x, y, xlab = '', 
             fill = x, add ="boxplot", palette = col_pal,...) +
      stat_compare_means(comparisons = compare, method = 't.test')+theme(legend.position = '' )
  }
}

#' @title f_boxplot2
#' @param data 
#' @param x 设置为x轴的列名
#' @param y 设置为y轴的列名
#' @param col_pal 设置颜色和颜色对应的组名, 默认为NA, 自动根据x轴变量数n选择rainbow(n)作为颜色
#' @param filter_y 输入一个值，筛选基因表达水平大于这个值的样本， 不能为NA
#' @param f 输入函数对y轴进行处理，常用log2, log10等, 默认为NA，即不处理
#' @param plot 'b'表示boxplot, 'v'表示violin, 默认为b
#' @param compare 对哪些组进行统计学检验，输入为列表
#' @param p.method 统计检验方法，默认为t.test, 其余选项有wilcox.test, anova, kruskal.test 
f_boxplot2 <- function(data, x, y,col_pal = NA, filter_y = NA, f = NA, plot = 'b',p.method = 't.test',compare = NA,verbose = T,...) {
  plot <- match.arg(plot, c('b','v'))
  p.method <- match.arg(p.method,c('wilcox.test','anova','t.test','kruskal.test'))
  library(ggpubr)
  
  if (all(is.na(col_pal))) {
    col_pal <- rainbow(length(unique(data[[x]])))
    names(col_pal) <- sort(unique(data[[x]]))
  }
  if (!is.na(filter_y)) data <- data[which(data[[y]]> filter_y),]
  if (typeof(f) %in% c('closure', 'builtin')) data[[y]] <- f(data[[y]])
  data[[x]] <- factor(data[[x]], levels = names(col_pal))
  
  if (verbose == T) {
    sum_mean <- aggregate(data[[y]], by = list(data[[x]]), mean)
    colnames(sum_mean) <- c('group','mean')  
    print(sum_mean)
  }
  if (plot == 'b') {
    p <- ggboxplot(data, x, y, xlab = '', 
                   fill = x, add ="jitter", add.params = list(color = 'grey', alpha = 0.8), palette = col_pal,...) +
      theme(legend.position = '' )
  } else if (plot == 'v') {
    p <- ggviolin(data, x, y, xlab = '', 
                  fill = x, add ="boxplot", palette = col_pal,...) +
      theme(legend.position = '' )
  }
  
  if (p.method %in% c('wilcox.test','t.test')) p <- p + stat_compare_means(comparisons = compare, method = p.method)
  if (p.method %in% c('anova','kruskal.test')) p <- p + stat_compare_means(method = p.method)
  return(p)
}


f_patch <- function(items, f_in, ncol =3, ...) {
  library(patchwork)
  for (i in 1:length(items)) {
    assign(paste0('p',i), f_in(items[i],...))
  }
  exp_text <- paste(paste0('p', 1:length(items)), collapse = '+') 
  eval(parse(text = exp_text)) + plot_layout(ncol = ncol)
}




