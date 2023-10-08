#' @title  CalculateDEG
#' @description calculate the FC and p-value of different expression proteins between groups
#' @importFrom limma
#' @importFrom dplyr

CalculateDEG.wilcox.paired <- function(data, id1, id2) {
  result <- data.frame(symbol = rownames(data), difference = NA, p_value = NA, FDR = NA)
  rownames(result) <- rownames(data)
  result$difference <- apply(data[,id1], 1, median) - apply(data[, id2], 1, median)
  
  for (i in rownames(data)) {
    test1 <- wilcox.test(x = t(data[i, id1]), y = t(data_tn[i, id2]), paired=TRUE, alternative = "two.sided")
    result[i, 'p_value'] <- test1$p.value
  }
  result$FDR <- p.adjust(result$p_value)
  return(result)
}

# data = exp_pro_lms; rowname = 0;group = clinfo1; 
# group.exp = 'firmiana_T'; group.level = 'ccp.group_3';
# class1 = 2; class2 = c(1,3); class1_name = '2'; class2_name = '13'

CalculateDEG.t <- function(data,rowname=0,group,group.exp=1,group.level=2,class1, class2,  class1_name = NA,class2_name = NA) {
  ###function: calculate the FC between class2 & class1(class2/class1) and p.value
  
  ###parameter:
  ##data:a DataFrame or a matrix
  ##rowname: the number of the column which should be regarded as rownames
  #rowname=0 means row indexes of the data are the real rownames
  ##group: a dataframe with 2 columns, 
  #the first one is Experiment names matching column names of data
  #the second is the group information
  ##group_exp,group_level: the column name or number of experiment number and factor
  ##class1,class2: the selected level in the group
  require(dplyr)
  data <- as.data.frame(data)
  group2 <- rep(NA,ncol(data))
  for (i in seq_along(data)) {
    num <- which(group[,group.exp] == colnames(data)[i])
    if(length(num)) {
      group2[i] <- group[num,group.level]
    }
  } 
  
  {if (length(class1) == 1) index1 <- which(group2==class1)
    else index1 <- which(group2 %in% class1)
  }
  {if (length(class2) == 1) index2 <- which(group2==class2)
    else index2 <-  which(group2 %in% class2)
  }
  
  
  if (rowname == 0){
    Symbol <- rownames(data)
  } else Symbol <- data[,rowname]
  
  if (is.na(class1_name)) class1_name <- class1
  if (is.na(class2_name)) class2_name <- class2
  
  
  
  subdata <- data[,c(index1,index2)]

  #t.test, two tails
  pvalue <- apply(subdata, 1,function(x) {
    group_f <- factor(c(rep('A',length(index1)), rep('B',length(index2))))
    fvalue <-  var.test(x~group_f) #F test
    if (!is.na(fvalue$p.value)) {
      if (fvalue$p.value >0.05) {
        t.test(x~group_f,var.equal =T)
      } else t.test(x~group_f,var.equal = F)
    } else  return(NA)
  }
  )
  
  # calucutae FC 
  gene.na <- which(is.na(pvalue))
  if(length(gene.na) >0){
    pvalue <- pvalue[-gene.na]
    data <- data[-gene.na,]
    Symbol <- Symbol[-gene.na]
  }
  
  result_ttest <- data%>%transmute(mean1 = rowMeans(dplyr::select(.,all_of(index1))),
                         mean2 = rowMeans(dplyr::select(.,all_of(index2))),
                         FC = mean1/mean2,
                         pvalue = as.numeric(unlist(lapply(pvalue, function(x) x$p.value)))
                         )
  
  result_ttest$Symbol <- Symbol
  result_ttest <- result_ttest %>% dplyr::select(Symbol,everything())
  colnames(result_ttest) <- c('Gene',
                              paste('mean',class1_name,sep=''),
                              paste('mean',class2_name,sep = ""),
                              paste('FC',paste(class1_name,'/',class2_name,sep = ''),sep = "_"),
                              'pvalue')
  result_ttest$FDR <- p.adjust(result_ttest$pvalue)
  return(result_ttest)
}



# data1 <- CalculateDEG.t(data = ExpData, rowname = 0, group = clinfo1,
#                     group.exp = 'firmiana_T',group.level = 'ccp.group_3',
#                     class1 = 1, class2 = 2)
# 
# data1_filt <- data1 %>% filter(t.test.p <=0.05)

CalculateDEG.limma <- function(data, group, group.exp=1,group.level=2,class1, class2,class1_name=NA,class2_name =NA,log2T = T) {
  require(dplyr)
  require(limma)
  data <- as.data.frame(data)
  group2 <- rep(NA,ncol(data))
  for (i in seq_along(data)) {
    num <- which(group[,group.exp] == colnames(data)[i])
    if(length(num)) {
      group2[i] <- group[num,group.level]
    }
  } 
  
  {if (length(class1) == 1) index1 <- which(group2==class1)
    else index1 <- which(group2 %in% class1)
  }
  {if (length(class2) == 1) index2 <- which(group2==class2)
    else index2 <-  which(group2 %in% class2)
  }
  if (is.na(class1_name)) class1_name <- class1
  if (is.na(class2_name)) class2_name <- class2
  
  
  data <- data[,c(index1,index2)]

  if (log2T == T ) data <- log2(data + 1)
  group_list <- factor(c(rep(class2_name,length(index2)), rep(class1_name,length(index1)))) 
  design <- model.matrix(~ group_list)
  colnames(design) <- levels(group_list)
  rownames(design) <- colnames(data)
  fit <- lmFit(data,design)
  fit <- eBayes(fit, trend = T)
  result_limma <- topTable(fit,coef = 2, n=Inf)
  return(result_limma)
  
}
#   
# data2 <- CalculateDEG.limma(data = ExpData,group = clinfo1,
#                     group.exp = 'firmiana_T',group.level = 'ccp.group_3',
#                     class1 = 1, class2 = 2,log2T = T)
# 
# data2_filt <- data2 %>% filter(P.Value <=0.05)

CalculateDEG.t.multi <- function(data,rowname=0,group,group.exp=1,group.level=2,class1, class2,  class1_name = NA,class2_name = NA) {
  result <- list()
  for (i in class2) {
    DEG <- CalculateDEG.t(data, rowname, group, group.exp, group.level, class1, class2 = i, class1_name, class1_name ) 
    result[[i]] <- DEG
  }
  return(result)
}

CalculateDEG.multi <- function(data,rowname=0,group,group.exp=1,group.level=2,class1, class2,  class1_name = NA,
                               class2_name = NA, method = 't', FC.method = '/') {
  result <- list()
  if (method == 't') {
    for (i in class2) {
      DEG <- CalculateDEG.t(data, rowname, group, group.exp, group.level, class1, class2 = i, class1_name, class2_name ) 
      result[[i]] <- DEG
    }
  }
  if (method == 'wilcox') {
    for (i in class2) {
      DEG <- CalculateDEG.wilcox(data, rowname, group, group.exp, group.level, class1, class2 = i, class1_name, class2_name, FC.method) 
      result[[i]] <- DEG
    }
  }
  return(result)
}



CalculateDEG.multiCycle <- function(data,rowname=0,group,group.exp=1,group.level=2,class = NA,
                                    class1_name = 'A',  class2_name = 'C',method = 't', FC.method = '/') {
  if (is.na(class)) class <- unique(group[,group.level])
  result <- list()
  for (i in class) {
    class1 <- i
    class2 <- setdiff(class,class1)
    if (method == 't') {
      DEG <- CalculateDEG.t(data, rowname, group, group.exp, group.level, class1, class2, class1_name, class2_name ) 
    }
    if (method == 'wilcox') {
      DEG <- CalculateDEG.wilcox(data, rowname, group, group.exp, group.level, class1, class2, class1_name, class2_name, FC.method) 
    }
    result[[i]] <- DEG
  }
  return(result)
}

CalculateDEG.multiCycleFilter <- function(DEP, cut.p = 0.05, cut.fc = 1.2) {
  for (i in names(DEP) ) {
    DEP[[i]] <- DEP[[i]][which( (DEP[[i]][,4] >= cut.fc) &(DEP[[i]][,5] <= cut.p)  ),]
  }
  return(DEP)
}

CalculateDEG.wilcox <- function(data,rowname=0,group,group.exp=1,group.level=2,
                                class1, class2,  class1_name = NA,class2_name = NA,
                                FC.method = c('/','-')) {
  FC.method <- match.arg(FC.method)
  require(dplyr)
  data <- as.data.frame(data)
  group2 <- rep(NA,ncol(data))
  for (i in seq_along(data)) {
    num <- which(group[,group.exp] == colnames(data)[i])
    if(length(num)) {
      group2[i] <- group[num,group.level]
    }
  } 
  
  {if (length(class1) == 1) index1 <- which(group2==class1)
    else index1 <- which(group2 %in% class1)
  }
  {if (length(class2) == 1) index2 <- which(group2==class2)
    else index2 <-  which(group2 %in% class2)
  }
  
  
  if (rowname == 0){
    Symbol <- rownames(data)
  } else Symbol <- data[,rowname]
  
  if (is.na(class1_name)) class1_name <- class1
  if (is.na(class2_name)) class2_name <- class2
  
  
  
  # subdata <- data[,c(index1,index2)]
  
  #wilcox.test, two tails
  p.value <- apply(data, 1,function(x) {
    # group_f <- factor(c(rep('A',length(index1)), rep('B',length(index2))))
    tryCatch(wilcox.test(x[index1], x[index2], exact = FALSE)$p.value,
             error = function(x){NA}
             )
  }
  )
  
  gene.na <- which(is.na(p.value))
  if(length(gene.na) >0){
    p.value <- p.value[-gene.na]
    data <- data[-gene.na,]
    Symbol <- Symbol[-gene.na]
  }
  
  # calucutae FC 
  if (FC.method == '/') {
     result_wilcox <- data%>%transmute(mean1 = rowMeans(dplyr::select(.,all_of(index1))),
                                      mean2 = rowMeans(dplyr::select(.,all_of(index2))),
                                      FC = mean1/mean2,
                                      pvalue = p.value
    )
  } else if(FC.method == '-'){
    result_wilcox <- data%>%transmute(mean1 = rowMeans(dplyr::select(.,all_of(index1))),
                                      mean2 = rowMeans(dplyr::select(.,all_of(index2))),
                                      FC = mean1 - mean2,
                                      pvalue = p.value
    )
  } else stop("please select FC.method from '/' or '-' ")
  
  result_wilcox$Symbol <- Symbol
  result_wilcox <- result_wilcox %>% dplyr::select(Symbol,everything())
  colnames(result_wilcox) <- c('Gene',
                               paste('mean',class1_name,sep=''),
                               paste('mean',class2_name,sep = ""),
                               paste('FC',paste(class1_name,'/',class2_name,sep = ''),sep = "_"),
                               'pvalue')
  result_wilcox$FDR <- p.adjust(result_wilcox$pvalue)
  return(result_wilcox)
}


CalculateDEG.t.norm <- function(data,rowname=0,group,group.exp=1,group.level=2,class1, class2,  class1_name = NA,class2_name = NA) {
  ###function: calculate the FC between class2 & class1(class2/class1) and p.value
  
  ###parameter:
  ##data:a DataFrame or a matrix
  ##rowname: the number of the column which should be regarded as rownames
  #rowname=0 means row indexes of the data are the real rownames
  ##group: a dataframe with 2 columns, 
  #the first one is Experiment names matching column names of data
  #the second is the group information
  ##group_exp,group_level: the column name or number of experiment number and factor
  ##class1,class2: the selected level in the group
  require(dplyr)
  data <- as.data.frame(data)
  group2 <- rep(NA,ncol(data))
  for (i in seq_along(data)) {
    num <- which(group[,group.exp] == colnames(data)[i])
    if(length(num)) {
      group2[i] <- group[num,group.level]
    }
  } 
  
  {if (length(class1) == 1) index1 <- which(group2==class1)
    else index1 <- which(group2 %in% class1)
  }
  {if (length(class2) == 1) index2 <- which(group2==class2)
    else index2 <-  which(group2 %in% class2)
  }
  
  
  if (rowname == 0){
    Symbol <- rownames(data)
  } else Symbol <- data[,rowname]
  
  if (is.na(class1_name)) class1_name <- class1
  if (is.na(class2_name)) class2_name <- class2
  
  
  
  subdata <- data[,c(index1,index2)]
  subdata <- apply(subdata, c(1,2),log2)
  norm.test <- rep(NA, nrow(subdata))
  names(norm.test) <- rownames(subdata)
  for (i in seq_along(norm.test)) {
    nt <- tryCatch(shapiro.test(as.numeric(subdata[i,]))$p.value, error = function(x){NA})
    norm.test[i] <- nt
  }
  norm.test <- apply(subdata,1, function(x) shapiro.test(as.numeric(x))$p.value)
  #t.test, two tails
  pvalue <- apply(subdata, 1,function(x) {
    group_f <- factor(c(rep('A',length(index1)), rep('B',length(index2))))
    fvalue <-  var.test(x~group_f) #F test
    if (!is.na(fvalue$p.value)) {
      if (fvalue$p.value >0.05) {
        t.test(x~group_f,var.equal =T)
      } else t.test(x~group_f,var.equal = F)
    } else  return(NA)
  }
  )
  
  # calucutae FC 
  gene.na <- which(is.na(pvalue))
  if(length(gene.na) >0){
    pvalue <- pvalue[-gene.na]
    data <- data[-gene.na,]
    Symbol <- Symbol[-gene.na]
  }
  
  result_ttest <- data%>%transmute(mean1 = rowMeans(dplyr::select(.,all_of(index1))),
                                   mean2 = rowMeans(dplyr::select(.,all_of(index2))),
                                   FC = mean1/mean2,
                                   t.test.p = as.numeric(unlist(lapply(pvalue, function(x) x$p.value)))
  )
  
  result_ttest$Symbol <- Symbol
  result_ttest <- result_ttest %>% dplyr::select(Symbol,everything())
  colnames(result_ttest) <- c('Gene',
                              paste('mean',class1_name,sep=''),
                              paste('mean',class2_name,sep = ""),
                              paste('FC',paste(class1_name,'/',class2_name,sep = ''),sep = "_"),
                              't.test.p')
  result_ttest$FDR <- p.adjust(result_ttest$t.test.p)
  return(result_ttest)
}



# CalculateEPhc <- function(data =proData, group =clin, class1) {
#   class2 <-  c('HC1','HC2','HC3','HC4','HC5','HC6')
#   class2 <- setdiff(class2, class1)
#   EPhc <- CalculateDEG.wilcox(data, group = group, 
#                               group.exp = 'firmianaT', group.level = 'his_hc',
#                               class1 = class1, class2 = class2, 
#                               class1_name = 'HC', class2_name = 'Ctl', FC.method = '-')
#   return(EPhc)
# }
# 
# # 按照FC和FDR筛选富集通路
# FilterEPhc <- function(EPhc){
#   EPhc <- EPhc %>% filter(`FC_HC/Ctl` > 0 & FDR <= 0.05)
#   return(EPhc)
# }
# # EPhc1 <- FilterEPhc(EPhc1)
# # EPhc2 <- FilterEPhc(EPhc2)
# # EPhc3 <- FilterEPhc(EPhc3)
# # EPhc4 <- FilterEPhc(EPhc4)
# # EPhc5 <- FilterEPhc(EPhc5)
# # EPhc6 <- FilterEPhc(EPhc6)
# 
# # 计算数据的分组均值
# CalculateGroupRowMean <- function(data, group, group.exp, group.level){
#   level <- unique(group[,group.level])
#   meanTab <- matrix(rep(NA, nrow(data)*length(level)), nrow = nrow(data) )
#   rownames(meanTab) <- rownames(data); colnames(meanTab) <- level
#   for (i in level) {
#     index <- group[[group.exp]][group[[group.level]] == i]
#     meanOneClass <- rowMeans(data[,index])
#     meanTab[,i] <- meanOneClass
#   }
#   return(as.data.frame(meanTab))
# }