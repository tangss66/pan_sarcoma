#' @description calculate  AUC of proteins

CalDEGauc <- function(DEG, ExpData, clinfo = clinfo, n1, c1,c2) {
  library(pROC)
  DEG_auc <- data.frame(DEG$Gene, rep(NA, nrow(DEG)))
  colnames(DEG_auc) <- c('Gene', 'AUC')
  rownames(DEG_auc) <- DEG_auc$Gene
  for (i in rownames(DEG_auc)) {
    control <- clinfo$firmianaT[(clinfo$ccp.group_3 == c1 |clinfo$ccp.group_3 == c2)]
    case <- clinfo$firmianaT[clinfo$ccp.group_3 == n1]
    rocobj <- roc(controls = c(t(ExpData[i,control])), cases = c(t(ExpData[i, case])),
                  smooth = F, direction = '<')
    DEG_auc[i,'AUC'] <- as.numeric(auc(rocobj))
  }
  return(DEG_auc)
}