#' @title Imputation
#' @description  several methods to imputate null values of the matrix

Imputation.knn <- function(data,k=10,scale=T,meth='weighAvg',distData=NULL){
  ###function: interpolation null values using knn method
  
  ###parameter:
  ##data:a matrix with numeric values 
  #the rows without null values in the matrix should be more than k
  ##scale,meth,distData: parameters of knnImputation
  
  ###warning: this function will delete the columns without any non empty values 
  #and the rows with less than 2 non empty values 
  require(DMwR2)
  # iszero=sapply(data,FUN=function(x){return(x==0)},simplify = T)
  # data[iszero]=NA
  #删除全为空值或只有一个值的行
  nacount_row=apply(data,1,function(x) length(which(is.na(x))))
  data=data[-which(nacount_row>=ncol(data)-1),]
  #删除全为空值的列
  nacount_col=apply(data,2,function(x) length(which(is.na(x))))
  if (length(which(nacount_col>=(nrow(data)-1)))){
    data=data[,-which(nacount_col>=(nrow(data)-1))]}
  k_max=length(which(complete.cases(as.matrix(data))))
  if (k_max<k&k_max>1){
    k=k_max
    print(paste('the max k is',k,sep = ' '))}
  if (k_max<=1){
    print('the max k is less than 1')}
  data_knn=DMwR2::knnImputation(data,k=k,scale=scale,meth=meth,distData = distData)
  return(data_knn)}


Imputation.quantile <- function(data){
  ###function: normalize the data using the quantile method
  ###parameter:
  ##data:a matrix or dataframe with numeric values
  require(preprocessCore)
  data2 <- normalize.quantiles(as.matrix(data))
  colnames(data2) <- colnames(data)
  rownames(data2) <- rownames(data)
  return(data2)
}