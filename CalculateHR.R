#' @title CalculateHD
#' @description Calculate hazard ratio
#' @importFrom survival
#' @importFrom survminer

# HDvolcano_fib<- clinfo_fib[,c('dfs.time', 'dfs.livingstatus')]
# rownames(HDvolcano_fib) <- HDvolcano_fib$firmiana_T
# HDvolcano_fib <- cbind(HDvolcano_fib, t(ExpData_fib))


CalculateHR <- function(data, time = 'os.time', status = 'os.livingstatus') {
  require(survival)
  # require(survminer)
  genes <- colnames(data)[-which(colnames(data) %in% c(time, status))]
  HRtable <- data.frame() #建立空白数据框用于存结果
  for(i in genes) {
    cox <- coxph(Surv(data[,time], data[,status]) ~ data[,i])
    coxSummary = summary(cox)
    HRtable =rbind(HRtable ,
                 cbind(id=i,
                       HR = coxSummary$conf.int[,"exp(coef)"],
                       HR.95L = coxSummary$conf.int[,"lower .95"],
                       HR.95H = coxSummary$conf.int[,"upper .95"],
                       pvalue = coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
  HRtable[,2:5] <- apply(HRtable[,2:5], c(1,2), as.numeric)
  return(HRtable)
}
