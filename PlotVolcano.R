#' @title  PlotValcano
#' @author shi-nian
#' @update 20221013

# old version ---------
PlotVolcano_cor <-  function(data, col.name, col.r, col.p, cut.off.p = 0.05, cut.off.r = 0,label = T, label.text, 
                             x.name = 'r', color.set=c("#546de5", "#d2dae2","#ff4757"), max_overlap = 20) {
  #Plot correlation volcano figure
  require(ggplot2)
  require(ggrepel)
  require(stringr)
  
  data[['change']] <- ifelse(data[[col.p]] <= cut.off.p, ifelse(data[[col.r]] > cut.off.r, 'Pos', 'Neg'), 'No')
  data[[col.p]] <- sapply(data[[col.p]], function(x) -log10(x), simplify = T)
  
  
  p<-ggplot(data, aes_string(x = col.r, y = col.p, color = 'change')) +
    geom_jitter(alpha=1,size=2) +
    scale_color_manual(values=color.set) +
    geom_hline(yintercept = -log10(cut.off.p), lty='dotdash',col="black",lwd=0.8) +
    labs(x= x.name, y="-log10 (p-value)") + scale_x_continuous(limits = c(-1, 1)) + theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), 
          legend.position="right", 
          legend.title = element_blank())
  if(label == T){
    data[['label']] = ifelse(data[[col.p]] > -log10(cut.off.p), as.character(data[[col.name]]), "")
  } else{
    data[['label']] <- NA
    # index <- sapply(label.text, function(x) which(data[[col.name]] == x), simplify = T)
    index <- which(data[[col.name]] %in% label.text)
    data$label[index]=data[[col.name]][index]
  }
  p <- p + geom_text_repel(data = data, aes_string(x = col.r, 
                                                   y = col.p, 
                                                   label = 'label'),
                           size = 3,box.padding = unit(0.5, "lines"),
                           point.padding = unit(0.8, "lines"), 
                           segment.color = "black", 
                           show.legend = FALSE,max.overlaps = max_overlap)
  p
}


PlotVolcano_DEG <- function(data,col.name='Symbol',col.p='p.value',col.FC='FC',doFCLog=F,label=T,label.text,
                        cut.off.pvalue=0.05,cut.off.FC=2,color.set=c("#546de5", "#d2dae2","#ff4757"), 
                        point.size = 2, x.labs = "log2(fold change)", max_overlap = 20){
  ###function: vacano_plot
  ###parameter:
  ##data: DataFrame formation
  ##col_name: original row names
  ##col_p: the column names of the p value
  ##col_FC: the column names of fold change
  ##doFCLog: whether transform FC to log2FC
  #cut_off_pvalue: the threshold value of the significance
  #cut_off_FC: the threshold value to identify different expression genes
  #color_set : the color of down,stable and up
  library(ggplot2)
  library(ggrepel)
  library(stringr)
  cut_off_logFC<-log2(cut.off.FC)
  #rename
  colnames(data)[which(colnames(data)==col.name)]="Symbol"
  colnames(data)[which(colnames(data)==col.p)]="p.value"
  colnames(data)[which(colnames(data)==col.FC)]="logFC"
  # log2 translate
  if(doFCLog==T){
    data$logFC=sapply(data$logFC,log2,simplify = T)
  }
  data$change<-ifelse(data$p.value<=cut.off.pvalue &abs(data$logFC)>=cut_off_logFC, ifelse(data$logFC>=cut_off_logFC,"Up","Down"), "Stable")
  
  p<-ggplot(data,aes(x=logFC,y=-log10(p.value),color=change))+
    geom_point(alpha=1,size= point.size)+
    scale_color_manual(values=color.set)+
    geom_vline(xintercept=c(-cut_off_logFC,cut_off_logFC),lty='dotdash',col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(cut.off.pvalue),lty='dotdash',col="black",lwd=0.8) +
    labs(x= x.labs, y="-log10 (p-value)")+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5), 
          legend.position="right", 
          legend.title = element_blank())
  {if(label==T){
    data$label = ifelse(data$p.value < cut.off.pvalue & abs(data$logFC) >= 1, as.character(data$Symbol),"")}
    else{
      data$label=NA
      index = sapply(label.text,function(x) which(data$Symbol==x),simplify = T)
      index = unlist(index)
      data$label[index]=data$Symbol[index]}}
  p=p+geom_text_repel(data = data, aes(x = data$logFC, 
                                       y = -log10(data$p.value), 
                                       label = label),
                      size = 3,box.padding = unit(0.5, "lines"),
                      point.padding = unit(0.8, "lines"), 
                      segment.color = "black", 
                      show.legend = FALSE,max.overlaps = max_overlap)
  p
  
}

PlotVolcano_HD <- function(data= HD, col.name='id', col.r ='HR', col.p='pvalue', cut.off.p = 0.05, cut.off.r = 1, label = T,label.text =NA, 
                           x.name = 'HR', color.set=c("#546de5", "#d2dae2","#ff4757"), max_overlap = 20) {
  
  require(ggplot2)
  require(ggrepel)
  require(stringr)
  
  data[['change']] <- ifelse(data[[col.p]] <= cut.off.p, ifelse(data[[col.r]] > cut.off.r, 'Pos', 'Neg'), 'No')
  data[[col.p]] <- sapply(data[[col.p]], as.numeric, simplify = T)
  data[[col.r]] <- sapply(data[[col.r]], as.numeric, simplify = T)
  data[[col.p]] <- sapply(data[[col.p]], function(x) -log10(x), simplify = T)
  
  
  p<-ggplot(data, aes_string(x = col.r, y = col.p, color = 'change')) +
    geom_jitter(height = 0.2,alpha=1,size=4) +
    scale_color_manual(values=color.set) +
    geom_hline(yintercept = -log10(cut.off.p), lty='dotdash',col="black",lwd=0.8) + geom_vline(xintercept = 1,lty='dotdash',col="black",lwd=0.8)+
    labs(x= x.name, y="-log10 (p-value)") + theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), 
          legend.position="right", 
          legend.title = element_blank())
  
  
  if(label == T){
    data[['label']] = ifelse(data[[col.p]] > -log10(cut.off.p), as.character(data[[col.name]]), "")
  } else{
    data[['label']] <- ''
    index <- which(data[[col.name]] %in% label.text)
    data$label[index]=data[[col.name]][index]
  }
  p <- p + geom_text_repel(data = data, aes_string(x = col.r, 
                                                   y = col.p, 
                                                   label = 'label'),
                           size = 3,box.padding = unit(0.5, "lines"),
                           point.padding = unit(0.8, "lines"), 
                           segment.color = "black", 
                           show.legend = FALSE,max.overlaps = max_overlap)
  p
  
}

PlotVolcano_HD2 <- function(hd_result) {
  data_vol <- hd_result
  
  color.set <- c( '#4dbbd5b2',"#d2dae2",'#f39b7fb2',"#3c5488b2",'#91d1c2b2','#dc0000b2')
  names(color.set) <- c('Neg','No','Pos','kinase','tf', 'oncogene')
  
  data_vol[['pvalue']] <- sapply(data_vol[['pvalue']], as.numeric, simplify = T)
  data_vol[['HR']] <- sapply(data_vol[['HR']], as.numeric, simplify = T)
  data_vol[['pvalue']] <- sapply(data_vol[['pvalue']], function(x) -log10(x), simplify = T)
  
  
  data_vol[['change']] <- ifelse(data_vol[['pvalue']] >= -log10(0.05), ifelse(data_vol[['HR']] > 1, 'Pos', 'Neg'), 'No')
  data_vol$change[which(data_vol$id %in% hd_kinase_result$id)] <- 'kinase'
  data_vol$change[which(data_vol$id %in% hd_tf_result$id)] <- 'tf'
  data_vol$change[which(data_vol$id %in% hd_onco_result$id)] <- 'oncogene'
  
  data_vol$change <- factor(data_vol$change, levels = c('Neg','No','Pos','kinase','tf', 'oncogene'))
  
  data_vol <- data_vol %>% arrange(change, by.group <- c('Neg','No','Pos','kinase','tf','oncogene'))
  p <-  ggplot(data_vol, aes_string(x = 'HR', y = 'pvalue', color = 'change')) +
    geom_jitter(alpha =1, size=3) +
    scale_color_manual(values=color.set) +
    geom_hline(yintercept = -log10(0.05), lty='dotdash',col="black",lwd=0.8) + geom_vline(xintercept = 1)+
    labs(x= 'HR', y="-log10 (p-value)") + theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), 
          legend.position="right", 
          legend.title = element_blank())
  
  data_vol[['label']] = ifelse(data_vol[['change']] %in% c('kinase','oncogene','tf'), as.character(data_vol[['id']]), "")
  p <- p + geom_text_repel(data = data_vol, aes_string(x = 'HR', 
                                                       y = 'pvalue', 
                                                       label = 'label'),
                           size = 3,box.padding = unit(0.5, "lines"),
                           point.padding = unit(0.8, "lines"), 
                           segment.color = "black", 
                           show.legend = FALSE,max.overlaps = 20)
  
  data_vol[['label2']] = ifelse(data_vol[['pvalue']] > -log10(0.05), as.character(data_vol[['id']]), "")
  
  p + geom_text_repel(data = data_vol, aes_string(x = 'HR', 
                                                  y = 'pvalue', 
                                                  label = 'label2'),
                      size = 3,box.padding = unit(0.5, "lines"),
                      point.padding = unit(0.8, "lines"), 
                      segment.color = "black", 
                      show.legend = FALSE,max.overlaps = 5)
  
}


# new version 20221013 -----------
#' @name PlotVolcano.processData
#' @description prepare data for volcano plot
#' @param data 
#' @param id the column name of gene
#' @param col.x the column name of x-axis
#' @param col.y the column name of y-axis
#' @param cut.x the threshold value of x to separate data
#' @param cut.x.dbl bool, default True, True will separate data by cut.x and -cut.x
#' @param cut.y the threshold value of y to separate data
#' @param group.auto bool, default True, True will group data by cut.x and cut.y
#' @param group.manual bool, default False,True will add new group manually
#' @param group.manual.text a list of new groups and ids them contains
#' @param label.auto bool, default True, add label based on cut.x and label.auto.cut.x
#' @param label.auto.cut.x 
#' @param label.manual bool, default False, True will add new group manually
#' @param label.manual.text 
PlotVolcano.processData <- function(data, id, col.x, col.y, cut.x, cut.x.dbl = T,cut.y = 0.05, 
                                   group.auto = T, group.manual = F, group.manual.text, 
                                   label.auto = T, label.auto.cut.x = NA,label.manual = F, label.manual.text = group.manual.text){
  data[[col.y]] <- sapply(data[[col.y]], function(x) -log10(x), simplify = T)
  
  data[['group']] <- 1
  if (group.auto == T) {
    if (cut.x.dbl == T) {
      data[['group']] <- ifelse(data[[col.y]] >= -log10(cut.y) & abs(data[[col.x]]) > cut.x , ifelse(data[[col.x]] > cut.x, 2, 0), 1)
    } else data[['group']] <- ifelse(data[[col.y]] >= -log10(cut.y), ifelse(data[[col.x]] > cut.x, 2, 0), 1)
  }
  if (group.manual == T) {
    if (!is.list(group.manual.text)) {
      stop('group.manual.text should be a list')
    } else {
      for (i in names(group.manual.text)) {
        index <- which(data[,id] %in% group.manual.text[[i]])
        data[index,'group'] <- i
      }
    }
  }
  
  data[['label']] <- ''
  if (is.na(label.auto.cut.x)) label.auto.cut.x <- cut.x
  if (label.auto == T) {
    if (cut.x.dbl == T) {
      data[['label']] <- ifelse(data[[col.y]] >= -log10(cut.y) & abs(data[[col.x]]) >= label.auto.cut.x, as.character(data[[id]]), '')
    } else data[['label']] <- ifelse(data[[col.y]] >= -log10(cut.y), as.character(data[[id]]), '')
  }
  if (label.manual == T) {
    for (i in names(label.manual.text)) {
      index <- which(data[,id] %in% label.manual.text[[i]])
      data[index,'label'] <- data[index, id]
    }
  }
  return(data)
}


#' @name PlotVolcano.processData
#' @param vline.dbl bool, default True , True will plot two vertical line at cut.x and -cut.x 
PlotVolcano.plt1 <- function(data, id, col.x, col.y, cut.x, cut.y = 0.05, col.group = 'group', 
                             col.label = 'label',point.size = 2, max_overlap = 20, vline.dbl = T, 
                             col.pal = c("#00AFBB", "#999999", "#FC4E07")) {
  require(ggplot2)
  require(ggrepel)
  data[,col.group] <- factor(data[,col.group])
  xint <- cut.x
  if (vline.dbl == T) xint <- c(-cut.x, cut.x)
  p <- ggplot(data, aes_string(x = col.x, y = col.y, color = col.group)) +
    geom_point(size = point.size) + geom_hline(yintercept = -log10(cut.y), lty = 'dashed') +
    geom_vline(xintercept = xint, lty = 'dashed') +
    theme_classic() + theme(legend.position="right",legend.title = element_blank()) + labs(y="-log10 (p-value)") +
    geom_text_repel(aes_string(label = col.label), size = 3,box.padding = unit(0.5, "lines"),
                    point.padding = unit(0.8, "lines"),segment.color = "black", show.legend = FALSE,
                    max.overlaps = max_overlap) 
  if (length(unique(data[,'group'])) == length(col.pal)){
    p <- p + scale_color_manual(values=col.pal)
  } else print("the number of colors don't equal th number of group")
  return(p)
} 


PlotVolcano.processPlt <- function(data, id, col.x, col.y, cut.x, cut.x.dbl = T, cut.y = 0.05,
                                   group.auto = T, group.manual = F, group.manual.text, 
                                   label.auto = T, label.auto.cut.x = NA,label.manual = F, label.manual.text = group.manual.text,
                                   point.size = 2, max_overlap = 20, vline.dbl = T,
                                   col.pal = c("#00AFBB", "#999999", "#FC4E07")) {
  processData <- PlotVolcano.processData(data, id, col.x, col.y, cut.x, cut.x.dbl, cut.y, group.auto, group.manual, group.manual.text, 
                          label.auto, label.auto.cut.x,label.manual, label.manual.text)
  PlotVolcano.plt1(processData, id, col.x, col.y, cut.x, cut.y, col.group = 'group', 
                   col.label = 'label',point.size,max_overlap,vline.dbl,col.pal)
}



# example ------------
# DEG(different expression protein) ------
# ## prepare data & parameters
# library(ggVolcano)
# data("deg_data", package = "ggVolcano")
# head(deg_data)
# 
# ## plot
# args1 <- list(data = deg_data, id = 'row', col.x ='log2FoldChange', col.y ='padj', cut.x = 1, label.auto.cut.x = 2)
# processData <- do.call(PlotVolcano.processData, args1)
# PlotVolcano.plt1(data = processData, id = 'row', col.x ='log2FoldChange', col.y ='padj', cut.x = 1)
# do.call(PlotVolcano.processPlt, args1)
# 
# # correlation data ---------
# ## prepare data & parameters
# data('GeneExp', package = "CancerSubtypes")
# CalculateCorr.gene <- function(data, geneName1, geneName2, corr.method = 'pearson') {
#   CorrData <- as.data.frame(matrix(rep(NA, length(geneName1)*length(geneName2)*4), ncol = 4))
#   colnames(CorrData) <- c('gene1','gene2', 'corr.r', 'corr.p')
#   CorrData[['gene1']] <- rep(geneName1, each = length(geneName2))
#   CorrData[['gene2']] <- rep(geneName2, length(geneName1))
#   options(warn = -1)
#   for (i in seq(nrow(CorrData))) {
#     gene1 <- CorrData[['gene1']][i]
#     gene2 <- CorrData[['gene2']][i]
#     corr <- cor.test(as.vector(t(data[gene1, ])), as.vector(t(data[gene2, ])), method = corr.method)
#     CorrData[['corr.r']][i] <- corr$estimate
#     CorrData[['corr.p']][i] <- corr$p.value
#   }
#   options(warn = 0)
#   return(CorrData)
# }
# dataCorr <- CalculateCorr.gene(GeneExp, geneName1 = 'CCK', geneName2 = rownames(GeneExp)[-1])
# 
# ## plot
# args2 <- list(data = dataCorr, id = 'gene2', col.x ='corr.r', col.y ='corr.p', cut.x = 0, cut.x.dbl =F)
# processdDataCorr <- do.call(PlotVolcano.processData, args2)
# PlotVolcano.plt1(data = processdDataCorr, id = 'gene2', col.x ='corr.r', col.y ='corr.p', cut.x = 0)
# do.call(PlotVolcano.processPlt, c(args2,list(vline.dbl = F, max_overlap = 100)))
# 
# 
# # hazard ratio ----------
# ## prepare data
# CalculateHR <- function(data, time = 'os.time', status = 'os.livingstatus') {
#   require(survival)
#   # require(survminer)
#   genes <- colnames(data)[-which(colnames(data) %in% c(time, status))]
#   HRtable <- data.frame() #建立空白数据框用于存结果
#   for(i in genes) {
#     cox <- coxph(Surv(data[,time], data[,status]) ~ data[,i])
#     coxSummary = summary(cox)
#     HRtable =rbind(HRtable ,
#                    cbind(id=i,
#                          HR = coxSummary$conf.int[,"exp(coef)"],
#                          HR.95L = coxSummary$conf.int[,"lower .95"],
#                          HR.95H = coxSummary$conf.int[,"upper .95"],
#                          pvalue = coxSummary$coefficients[,"Pr(>|z|)"])
#     )
#   }
#   HRtable[,2:5] <- apply(HRtable[,2:5], c(1,2), as.numeric)
#   return(HRtable)
# }
# 
# 
# data('myeloma', package = "survminer")
# data <- myeloma
# data1 <- t(GeneExp[1:100,])
# data <- data[1:nrow(data1),]
# rownames(data) <- c(1:nrow(data1)) 
# rownames(data1) <- c(1:nrow(data1)) 
# data2 <- cbind(data,data1)
# data2 <- data2[,-c(1:3)]
# data2$time <- as.numeric(data2$time)
# dataHR <- CalculateHR(data2, time = 'time', status = 'event')
# 
# 
# ## plot
# args3 <- list(data = dataHR, id = 'id', col.x ='HR', col.y ='pvalue', cut.x = 1, cut.x.dbl =F)
# processdDataHR <- do.call(PlotVolcano.processData, args3)
# PlotVolcano.plt1(data = processdDataHR, id = 'id', col.x ='HR', col.y ='pvalue', cut.x = 1, vline.dbl = F)
# do.call(PlotVolcano.processPlt, c(args3,list(vline.dbl = F, col.pal = c("#00AFBB", "#999999") )))
# 
# 
# # 自定义组别(group)和标签(label)
# ## sample1
# data('cc.genes', package = 'Seurat')
# gene.s <- intersect(cc.genes$s.genes, deg_data$row)
# gene.g2m <- intersect(cc.genes$g2m.genes, deg_data$row)
# group1 <- list('cell_cycle_s' = gene.s, 'cell_cycle_g2m' = gene.g2m)
# 
# args4 <- list(group.auto =F,group.manual = T, group.manual.text = group1, label.auto =F, 
#               label.manual =T, label.manual.text = group1)
# 
# do.call(PlotVolcano.processPlt, c(args1, args4,list(col.pal = c("#99999960", '#BEBADA', '#FB8072'), max_overlap = 1000) ))
# result1 <- do.call(PlotVolcano.processData, c(args1, args4))
# 
# ## sample2
# library(dplyr)
# deg_data <- deg_data %>% arrange(padj)
# group2 <- list('term1' = deg_data[1:20,'row'], 'term2' = deg_data[21:50,'row'])
# args5 <- list(group.manual = T, group.manual.text = group2, label.auto =F, 
#               label.manual =T, label.manual.text = group2, col.pal = c("#00AFBB", "#99999960", "#FC4E07",'#BEBADA', '#FB8072'))
# do.call(PlotVolcano.processPlt, c(args1, args5))
# 
