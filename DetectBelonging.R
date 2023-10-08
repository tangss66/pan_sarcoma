#' @title DetectBelonging
#' @description detect if values in one set belong to another

DetectBelonging.list <- function(current_data, target_data, by.current, by.target) {
  data <- current_data[[by.current]] %in% target_data[[by.target]]
  data[which(isTRUE(data))] = 1
  data[which(isFALSE(data))] = 0
  return(data)
}


BelongingBatch <- function(data, by = 'id') {
  # load('D:/Projects/Pan-sarcoma/AnalysisBranch/RData/DataBase.RData')
  load('D:/tang/projects/pan_sarcoma/workflow/rda/DataBase.RData')
  data$kinase <- DetectBelonging.list(data, db_kinase, by,'HGNC.Name')
  data$tf <- DetectBelonging.list(data, db_tf, by,'Symbol')
  data$tf_cofactor <- DetectBelonging.list(data, db_tf_cofactor, by,'Symbol')
  data$oncogene <- DetectBelonging.list(data, db_oncogene, by,'OncogeneName')
  # data$secreted <- DetectBelonging.list(data, db_secretedPro, by, 'Gene')
  # data$ctAntigen <- DetectBelonging.list(data, db_ctAnti, by, 'Family.member')
  # data <- merge(data, db_drug[,c('TARGTYPE','SYMBOL')], by.x = by, by.y = 'SYMBOL', all.x = T)
  return(data)
}

# DetectSubstrate <- function(data, by = ) {
#   #find substrates and phosphorylation sites for kinase
#   load('D:/Projects/Pan-sarcoma/Analysis/R_function/DataBase.RData')
#   
# }
# 
# DetectTG <- function(data, by = ) {
#   #find target genes for transcription factors
#   
# }