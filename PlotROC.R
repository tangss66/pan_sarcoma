#' @title  PlotROC

library(pROC)
library(ggplot2)

# 生成ROC表达式
plotROC.expression <- function(data, value, outcome, smooth = F) {
  require(pROC)
  data_str <- deparse(substitute(data))
  rocobj_exp <- paste0('roc(', outcome, '~')
  rocobj_exp <- paste0(rocobj_exp, value[1])
  if (length(value) != 1) {
    for(i in value[-1]) rocobj_exp <- paste0(rocobj_exp, '+',i)
  } 
  rocobj_exp <- paste0(rocobj_exp, ',data=', data_str, ',ci=T ')
  if (smooth == T) rocobj_exp <- paste0(rocobj_exp, ', smooth = T')
  rocobj_exp <- paste0(rocobj_exp,')')
  return(rocobj_exp)
}


# 执行ROC表达式
ROC <- function(...){
  rocobj_exp <- plotROC.expression(...)
  eval(parse(text = rocobj_exp))
}


# 批量计算ROC
calculateROC <- function(data, items, outcome, smooth = F){
  n <- length(items)
  result <- data.frame(items = items, outcome = rep(outcome,n), 
                       auc = rep(NA, n))
  rownames(result) <- items
  for (value in items){
    rocobj <- tryCatch(
      ROC(data = data, value = value, outcome = outcome),
      error = function(e) list('auc' = NA)
    )
    rocobj_auc <- round(as.numeric(rocobj$auc),3)
    result[value,'auc'] <- rocobj_auc
  }
  return(result)
}

# 提取ROC曲线的AUC和95%CI
plotROC.auc <- function(rocobj) {
  if (names(rocobj)[1] == 'percent') {
    print('only one  variant')
    ci_value <- as.numeric(rocobj$ci)
    result <- ci_value
  } else {
    result <- data.frame()
    for (i in names(rocobj)) {
      ci_value <- as.numeric(rocobj[[i]]$ci)
      result <- rbind(result, ci_value)
    }
  }
  names(result) <- c('95.low','auc','95.up')
  return(result)
}


# 生成图形注释，包括AUC和95%CI
plotROC.ano <- function(rocobj) {
  anos <- ''
  if (names(rocobj)[1] == 'percent') {
    ci_value <- as.numeric(rocobj$ci)
    ci_value <- round(ci_value,2)
    anos <- paste0(' AUC = ', ci_value[2], '\n ', '95% CI:',
                  ci_value[1] , '-',  ci_value[3])
  } else {
    for (i in names(rocobj)) {
      ci_value <- as.numeric(rocobj[[i]]$ci)
      ci_value <- round(ci_value,2)
      ano <- paste0(i, ': AUC = ', ci_value[2], '\n ', '95% CI:',
                    ci_value[1] , '-',  ci_value[3])
      anos <- paste(anos, sep = '\n \n', ano)
    }
  }
  return(anos)
}



# 绘图
plotROC.plt <- function(rocobj, ano = T) {
  require(ggplot2)
  p <- ggroc(rocobj, aes=c("color"))+theme_classic() + 
       ggsci::scale_color_lancet() +
       theme(legend.title = element_blank()) +
       geom_abline(intercept = 1, slope = 1,lty='longdash')
  if (ano == T) {
    anos <- plotROC.ano(rocobj)
    p <- p + annotate("text", x = 0.25, y = 0.3, label = anos, size =3)
  }
  return(p)
}


# ggroc(rocobj2, color = 'red')+theme_classic() +
#   annotate("text", x = 0.25, y = 0.3, label = ano) + 
#   ggsci::scale_color_lancet() +
#   theme(legend.title = element_blank()) +
#   geom_abline(intercept = 1, slope = 1,lty='longdash')



# 一键完成所有操作
plotROC.onestep <- function(...) {
  rocobj <-  ROC(...)
  tb.auc <- plotROC.auc(rocobj)
  p <- plotROC.plt(rocobj)
  return(list('auc' = tb.auc, 'plot' = p))
}



# example -----
# data(aSAH)
# x <- plotROC.onestep(aSAH, c('age', 's100b'), 'outcome')
# x
# x1 <- plotROC.onestep(aSAH, c('age', 's100b'), 'outcome', smooth =T)
# x1
# plotROC.onestep(aSAH, c('age'), 'outcome')
# plotROC.onestep(aSAH, c('age'), 'outcome') + ggsci::scale_color_aaas()
# 
# data(aSAH)
# rocobj2 <- ROC(aSAH, 'age', 'outcome')
# ano <- plotROC.ano(rocobj2)
# ggroc(rocobj2, color = 'red')+theme_classic() +
#   annotate("text", x = 0.25, y = 0.3, label = ano) +
#   theme(legend.title = element_blank()) +
#   geom_abline(intercept = 1, slope = 1,lty='longdash')
