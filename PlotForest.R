#' @title PlotForest
#' @description generate a forest plot of HR
#' @author shi-nian
#' @update 20221013
#' 
#' @param data 表达数据, 纵坐标为基因名，横坐标为样本名称
#' @param meta 生存数据矩阵
#' @param HRtable HR统计表
#' @param time 生存时间
#' @param status 生存状态
#' @param gene 基因名称
#' @param hr.low 95% CI 下限
#' @param hr.high 95% CI 上限
#' @param hr hazard ratio
#' @param pvalue hazard ratio' pvalue
#' @param col.pal 设定颜色, 包含三个元素的向量,顺序为good, no,poor


PlotForest.hr <- function(data, meta, id, time, status){
  require(survival)
  meta <- meta[,c(id,time,status)]
  data <- merge(meta,data,by.x = id, by.y = 0)
  data <- data[,-1]
  genes <- colnames(data)[-which(colnames(data) %in% c(time, status))]
  HRtable <- data.frame() #generate a blank data.frame to save the result
  for(i in genes) {
    cox <- coxph(Surv(data[,time], data[,status]) ~ data[,i])
    coxSummary = summary(cox)
    HRtable =rbind(HRtable ,
                   cbind(id=i,
                         HR = coxSummary$conf.int[,"exp(coef)"],
                         HR.95L = coxSummary$conf.int[,"lower .95"],
                         HR.95H = coxSummary$conf.int[,"upper .95"],
                         pvalue = coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
  HRtable[,2:5] <- apply(HRtable[,2:5], c(1,2), as.numeric)
  return(HRtable)
}

PlotForest.sortGene <- function(HRtable, gene = 'id',hr = 'HR', pvalue = 'pvalue', 
                                sort.method = c('na','p','p-re','hr','hr-re')) {
  require(dplyr)
  sort.method <- match.arg(sort.method, c('na','p','p-re','hr','hr-re'))
  if (sort.method == 'p')  index = order(HRtable[,pvalue])
  if (sort.method == 'p-re') index = order(desc(HRtable[,pvalue]))
  if (sort.method == 'hr') index = order(HRtable[,hr])
  if (sort.method == 'hr-re') index = order(desc(HRtable[,hr]))
  if (sort.method == 'na') index = seq(1:nrow(HRtable))
  HRtable <- HRtable[index,]
  HRtable[, gene] <- factor(HRtable[, gene], levels = HRtable[, gene]) 
  return(HRtable)
}

PlotForest.plt1 <- function(HRtable, gene = 'id', hr = 'HR', hr.low = 'HR.95L', hr.high = 'HR.95H', 
                            pvalue = 'pvalue',sort.method = 'na') {
  require(ggplot2)
  require(viridis)
  HRtable <- PlotForest.sortGene(HRtable, gene, hr, pvalue, sort.method)
  ggplot(data = HRtable,aes_string(x = hr, y= gene, color=pvalue)) +
    geom_vline(xintercept = 1,linetype='dashed',linewidth=1) +
    geom_errorbar(aes_string(xmax = hr.high, xmin = hr.low),color="black",width = 0.3,size=1)+
    geom_point(aes_string(x = hr, y = gene),size=6,shape=18) + #prismatic
    coord_trans(x='log2') +
    labs(y = "", x = 'Hazard Ratio')+  #标签
    labs(color="P value",title ="Univariate Cox regression analysis" )+
    scale_color_viridis() +
    theme_classic()
}


PlotForest.plt2 <- function(HRtable, gene = 'id', hr = 'HR', hr.low = 'HR.95L', hr.high = 'HR.95H', 
                            pvalue = 'pvalue', col.pal = c('#346D9C','grey','#BD0026'),sort.method = 'na') {
  require(ggplot2)
  require(viridis)
  HRtable <- PlotForest.sortGene(HRtable, gene, hr, pvalue, sort.method)
  names(col.pal) <- c('good','no','poor')
  HRtable$change <- ifelse(HRtable[,pvalue] <= 0.05, ifelse(HRtable[,hr] >1, 'poor', 'good'), 'no')
  uniLab <- sort(unique(HRtable$change))  
  col.pal <- col.pal[uniLab]
  ggplot(data = HRtable,aes_string(x = hr, y= gene, color= 'change')) +
    geom_vline(xintercept = 1,linetype='dashed',size=1) +
    geom_errorbar(aes_string(xmax = hr.high, xmin = hr.low),color="black",width = 0.3,linewidth=1)+
    geom_point(aes_string(x = hr, y = gene),size=4,shape=18) + #prismatic
    coord_trans(x='log2') +
    labs(y = "", x = 'Hazard Ratio')+  #标签
    labs(color="P value",title ="Univariate Cox regression analysis" )+
    scale_color_manual(values = col.pal) +
    theme_classic()
}


# example --------
# data 蛋白表达矩阵，纵坐标为基因名，横坐标为样本名称
# meta  生存信息
# id  样本编号
# time meta中的生存日期
# status meta中的生存状态
# proData <- read.csv('D:/tang/proData.csv')
# clin <- read.csv('D:/tang/clin.csv')
# HRtable <- PlotForest.hr(data = t(proData), meta = clin, id = 'firmianaT', time = 'os.time.month', status = 'os.livingstatus')
# HRtable1 <- HRtable %>% arrange(pvalue) %>% head(10)
# p <- PlotForest.plt1(HRtable1)
# p
# 
# 
# 
# 
# # 更改配色方案
# library(colorspace)
# p + colorspace:: scale_color_continuous_sequential(palette = 'Burg')
# library(viridis)
# p + viridis::scale_color_viridis(option = "D")
# p + xlim(-3,3)
# 
# # # 对基因进行排序
# PlotForest.plt1(HRtable1)
# PlotForest.plt2(HRtable1)
# PlotForest.plt1(HRtable1, sort.method = 'hr')
# PlotForest.plt2(HRtable1, sort.method = 'hr-r')

