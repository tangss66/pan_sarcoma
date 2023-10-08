#' @title  StatToolKit
#' @description  a toolkit for data summary and process
#' @author shi-nian
#' @update 2022-10-15



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


CalculateCorr.gene <- function(data, geneName1, geneName2, corr.method = 'pearson') {
  CorrData <- as.data.frame(matrix(rep(NA, length(geneName1)*length(geneName2)*4), ncol = 4))
  colnames(CorrData) <- c('gene1','gene2', 'corr.r', 'corr.p')
  CorrData[['gene1']] <- rep(geneName1, each = length(geneName2))
  CorrData[['gene2']] <- rep(geneName2, length(geneName1))
  options(warn = -1)
  for (i in seq(nrow(CorrData))) {
    gene1 <- CorrData[['gene1']][i]
    gene2 <- CorrData[['gene2']][i]
    corr <- cor.test(as.vector(t(data[gene1, ])), as.vector(t(data[gene2, ])), method = corr.method)
    CorrData[['corr.r']][i] <- corr$estimate
    CorrData[['corr.p']][i] <- corr$p.value
  }
  options(warn = 0)
  return(CorrData)
}
