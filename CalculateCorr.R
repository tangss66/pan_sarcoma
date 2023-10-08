#' @title  CalculateCorr
#' @description 



# CalculateCorr <- function() {
#   
# }

CalculateCorr.KinaseSub <- function(kinase,subsite, corr.method = 'pearson') {
  both_colnames <- intersect(colnames(kinase),colnames(subsite))
  kinase <- kinase[,both_colnames]
  subsite <- subsite[,both_colnames]
  CorrData <- as.data.frame(matrix(rep(NA, nrow(kinase)*nrow(subsite)*4), ncol = 4))
  colnames(CorrData) <- c('kinase','sites','corr.r','corr.p')
  CorrData[['kinase']] <- rep(rownames(kinase), each = nrow(subsite))
  CorrData[['sites']] <- rep(rownames(subsite), nrow(kinase))
  options(warn = -1)
  for (i in seq(nrow(CorrData))) {
    gene <- CorrData[['kinase']][i]
    site <- CorrData[['sites']][i]
    corr <- cor.test(as.vector(t(kinase[gene, ])), as.vector(t(subsite[site, colnames(kinase)])), method = corr.method)
    CorrData[['corr.r']][i] <- corr$estimate
    CorrData[['corr.p']][i] <- corr$p.value
  }
  options(warn = 0)
  return(CorrData)
}




CalculateCorr.gene <- function(data, geneName1, geneName2, corr.method = 'pearson', f = NA,...) {
  CorrData <- as.data.frame(matrix(rep(NA, length(geneName1)*length(geneName2)*4), ncol = 4))
  colnames(CorrData) <- c('gene1','gene2', 'corr.r', 'corr.p')
  CorrData[['gene1']] <- rep(geneName1, each = length(geneName2))
  CorrData[['gene2']] <- rep(geneName2, length(geneName1))
  options(warn = -1)
  dataC <- data
  for (i in seq(nrow(CorrData))) {
    gene1 <- CorrData[['gene1']][i]
    gene2 <- CorrData[['gene2']][i]
    if (typeof(f) %in% c('closure', 'builtin')) data <- f(dataC,x =gene1 ,y =gene2 ,...)
    corr <- cor.test(as.vector(t(data[gene1, ])), as.vector(t(data[gene2, ])), method = corr.method)
    CorrData[['corr.r']][i] <- corr$estimate
    CorrData[['corr.p']][i] <- corr$p.value
  }
  options(warn = 0)
  return(CorrData)
}

CalculateCorr.gene2 <- function(data, geneName1, geneName2, corr.method = 'pearson', f = NA,...) {
  CorrData <- as.data.frame(matrix(rep(NA, length(geneName1)*length(geneName2)*4), ncol = 4))
  colnames(CorrData) <- c('gene1','gene2', 'corr.r', 'corr.p')
  CorrData[['gene1']] <- rep(geneName1, each = length(geneName2))
  CorrData[['gene2']] <- rep(geneName2, length(geneName1))
  options(warn = -1)
  dataC <- data
  for (i in seq(nrow(CorrData))) {
    gene1 <- CorrData[['gene1']][i]
    gene2 <- CorrData[['gene2']][i]
    if (typeof(f) %in% c('closure', 'builtin')) data <- f(dataC,x =gene1 ,y =gene2 ,...)
    corr <- cor.test(data[,gene1], data[,gene2], method = corr.method)
    CorrData[['corr.r']][i] <- corr$estimate
    CorrData[['corr.p']][i] <- corr$p.value
  }
  options(warn = 0)
  return(CorrData)
}



#' @title f_scatter_add_coef
#' @description 用于给散点图的回归曲线添加回归系数和显著水平，由于现在的函数添加的回归系数在导入Adobe Illustrator时
#'              会发生字体格式的变化，因此通过annotate手动添加回归系数。
#' @param gg.obj ggplot对象
#' @param method 相关性计算方法
f_scatter_add_coef  <- function(gg.obj, method = 'pearson') {
  library(ggplot2)
  x <- as.character(gg.obj$mapping$x[[2]])
  y <- as.character(gg.obj$mapping$y[[2]])
  data <- gg.obj$data
  corResult <- cor.test(data[,x], data[,y], method = method)
  
  x_position <- fivenum(data[,x])[4]
  y_position <- fivenum(data[,y])[5]
  ano_text <- paste0('R = ', round(corResult$estimate, 2), '\n p = ', signif(corResult$p.value,2))
  gg.obj + annotate('text',x = x_position, y = y_position,label = ano_text)
}



