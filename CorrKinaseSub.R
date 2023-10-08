#' @title  CorrKinaseSub
#' @description calculate the correlation between kinase and substrates

CorrKinaseSub <- function(kinase_exp, sub_exp) {
  require(dplyr)
  expNum <- intersect(colnames(kinase_exp), colnames(sub_exp))
  kinase_exp <- t(kinase_exp[,expNum])
  sub_exp <- t(sub_exp[,expNum])
  result <- list()
  for (k in colnames(kinase_exp)) {
    data <- apply(sub_exp, 2, function(x, y = kinase_exp[,k]) cor.test(x, y, method = 'spearman'))
    pvalue <- lapply(data, function(x) return(x$p.value)) %>% unlist()
    estimate <- lapply(data, function(x) return(x$estimate)) %>% unlist()
    result[[k]] <- data.frame('substrate' = colnames(sub_exp),'pvalue' = pvalue, 'estimate' = estimate) %>% 
      mutate('kinase' = k, FDR = p.adjust(pvalue)) %>% filter(pvalue <= 0.05) %>% dplyr::select(kinase, everything()) 
  }
  return(result)
}

