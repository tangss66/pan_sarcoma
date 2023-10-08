#' @title steamnessScore
#' @description calculate steamness score
#' https://www.jianshu.com/p/efbd952e13e2

mRNAsi <- readRDS('D:/Projects/Pan-sarcoma/AnalysisBranch/mRNAsi.rds')

predict.mRNAsi <- function(exp, model) {
  
  common <- intersect(model$HUGO, rownames(exp))
  X <- exp[common, ]
  w <- model$Weight
  names(w) <- model$HUGO
  w <- w[common]
  w <- as.numeric(w)
  
  score <- apply(X, 2, function(z) {cor(z, w, method="sp", use="complete.obs")})
  score <- score - min(score)
  score <- score / max(score)
}

score_steamness <- predict.mRNAsi(ExpData, model = mRNAsi)