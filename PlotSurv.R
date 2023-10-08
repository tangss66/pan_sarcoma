#' @title PlotSurv
#' @description generate a plot of KM survival curve
#' @author shi-nian
#' @update 20221013 增加多重比较，添加注释等功能
#' @update 20221028 解决cox检验的p < 0.001 显示为0的问题; 给多重比较添加log-rank检验和Cox检验注释; 修该了错误的分组信息
#'
#' @param data 生存数据矩阵
#' @param time 生存时间
#' @param status 生存状态
#' @param yitem 分组数据
#' @param n  用n分位数分组, 仅当cut.method = 'ntile'时有用
#' @param cut.method 分组方式(best/ntile)
#' @param addCoxAnotation 是否添加cox检验值作为注释
#' @param col.pal 设定颜色


# fit <- survfit(Surv(os.time, os.livingstatus) ~ ccp.group_3, data = data)
# cutoff <- surv_cutpoint(shc1tf, time = 'dfs.time', event = 'dfs.livingstatus', variables = 'SHC1')
# cutcat <- surv_categorize(cutoff)
# fit<- survfit(Surv(dfs.time, dfs.livingstatus) ~ SHC1, data = cutcat)
# ggsurvplot(fit, data = cutcat, pval = T, pval.method = T, palette = col_hl)

# KM surv function------------------
PlotSurv <- function(fit, data = clinfo1, col.pal = "npg") {
  require(survival)
  require(survminer)
  ps <- ggsurvplot(fit,
    data = data,
    surv.median.line = "hv",
    pval = T, pval.size = 8,
    pval.method = T, pval.method.size = 8,
    conf.int = F,
    risk.table = F,
    palette = col.pal,
    legend = "right"
  )
  return(ps$plot)
}

PlotSurv2 <- function(data, time, status, yitem, col.pal = "npg", ...) {
  library(survival)
  library(survminer)
  tb <- table(data[[yitem]])
  tb_label <- paste(names(tb), tb, sep = "=")
  tb_label <- paste(tb_label, collapse = ";\n ")

  fit.text <- sprintf("surv_fit(Surv(%s, %s) ~ %s, data = data)", time, status, yitem)
  fit <- eval(parse(text = fit.text))

  ano <- surv_pvalue(fit)
  ano <- paste0(ano$method, ":", ano$pval.txt, "\n ", tb_label)

  psurv <- ggsurvplot(fit,
    data = data,
    surv.median.line = "hv",
    # pval = T,pval.size = 8,
    # pval.method = T, pval.method.size = 8,
    legend.title = yitem,
    # legend.labs=tb_lab,
    conf.int = F,
    risk.table = F,
    palette = unname(col.pal),
    legend = "right", ...
  )
  x.position <- quantile(data[, time], na.rm = T)[2]
  x.position <- unname(x.position)
  psurv$plot + annotate("text", x = x.position, y = 0.25, label = ano)
}


PlotSurv.groupBestPoint <- function(data, time, status, yitem) {
  require(survminer)
  cutoff <- surv_cutpoint(data, time, status, variables = yitem)
  cutcat <- surv_categorize(cutoff)
  cutcat[which(cutcat[, yitem] == "low"), yitem] <- 0
  cutcat[which(cutcat[, yitem] == "high"), yitem] <- 1
  return(cutcat)
}

PlotSurv.groupBestPoint2 <- function(..., yitem, yitem2) {
  require(survminer)
  cutcat <- PlotSurv.groupBestPoint(..., yitem)
  cutcat2 <- PlotSurv.groupBestPoint(..., yitem2)
  cutcat <- cbind(cutcat, cutcat2[, 3])
  colnames(cutcat)[4] <- colnames(cutcat2)[3]
  return(cutcat)
}

PlotSurv.groupNtile <- function(data, time, status, yitem, n = 2) {
  require(dplyr)
  cutcat <- data[, c(time, status, yitem)]
  cutcat[, yitem] <- ntile(cutcat[, yitem], n)
  cutcat <- cutcat[which(cutcat[, yitem] %in% c(1, n)), ]
  cutcat[which(cutcat[, yitem] == 1), yitem] <- 0
  cutcat[which(cutcat[, yitem] == n), yitem] <- 1
  return(cutcat)
}

PlotSurv.groupNtile2 <- function(data, time, status, yitem, yitem2, n = 2) {
  require(dplyr)
  cutcat <- data[, c(time, status, yitem, yitem2)]
  cutcat[, yitem] <- ntile(cutcat[, yitem], n)
  cutcat[, yitem2] <- ntile(cutcat[, yitem2], n)

  cutcat <- cutcat[which(cutcat[, yitem] %in% c(1, n)), ]
  cutcat <- cutcat[which(cutcat[, yitem2] %in% c(1, n)), ]

  cutcat[which(cutcat[, yitem] == 1), yitem] <- 0
  cutcat[which(cutcat[, yitem] == n), yitem] <- 1
  cutcat[which(cutcat[, yitem2] == 1), yitem2] <- 0
  cutcat[which(cutcat[, yitem2] == n), yitem2] <- 1
  return(cutcat)
}


PlotSurv.mergeGroup <- function(cutcat) {
  cutcat$group <- paste0(cutcat[, 3], cutcat[, 4])
  cutcat$group <- ifelse(cutcat$group == "00", "1", ifelse(cutcat$group == "10", "2", ifelse(cutcat$group == "01", "3", "4")))
  return(cutcat)
}

PlotSurv.pairwiseCom <- function(cutcat, time, status, group = "group") {
  require(survminer)
  require(survival)
  # fit.text <-  sprintf('pairwise_survdiff(Surv(%s, %s) ~  %s, data = cutcat)', time, status, group)
  # comp <- eval(parse(text = fit.text))
  colnames(cutcat)[1:2] <- c("time", "status")
  comp <- pairwise_survdiff(Surv(time, status) ~ group, data = cutcat)
  pvalue <- as.vector(unlist(comp$p.value))
  name <- array(dim = c(3, 3))
  for (i in 1:3) {
    for (j in 2:4) {
      name[i, j - 1] <- paste(i, j, sep = " vs ")
    }
  }
  pvalue_name <- as.vector(t(name))
  logrank <- data.frame(pvalue_name, pvalue)
  # drop data with p-value larger than 0.05
  logrank_sig <- subset(logrank, pvalue < 0.05)
  # if p < 0.0001, replace p-value by <0.0001
  logrank_sig$pvalue <- lapply(logrank_sig$pvalue, function(i) {
    ifelse(i < 0.0001, "<0.0001", round(i, 4))
  })

  return(logrank_sig)
}

PlotSurv.cox <- function(data, time, status, yitem) {
  cox <- coxph(Surv(data[, time], data[, status]) ~ data[, yitem])
  coxSummary <- summary(cox)
  return(coxSummary)
}

PlotSurv.addCoxAnotation <- function(signif_num, ...) {
  require(dplyr)
  coxSummary <- PlotSurv.cox(...)
  HR <- coxSummary$conf.int[, "exp(coef)"] %>%
    as.numeric() %>%
    signif(signif_num)
  HR.95L <- coxSummary$conf.int[, "lower .95"] %>%
    as.numeric() %>%
    signif(signif_num)
  HR.95H <- coxSummary$conf.int[, "upper .95"] %>%
    as.numeric() %>%
    signif(signif_num)
  pvalue <- coxSummary$coefficients[, "Pr(>|z|)"] %>%
    as.numeric() %>%
    signif(signif_num)
  anoCox <- paste0("HR: ", HR, "; 95%CI: ", HR.95L, "-", HR.95H, "\n pvalue = ", pvalue)
  return(anoCox)
}

PlotSurv.splt <- function(data, time, status, yitem, n = 2, # special_name = F,
                          cut.method = c("best", "ntile"), addCoxAnotation = F, col.pal = c("#1D91C0", "#F79EC0"), ...) {
  require(survival)
  require(survminer)
  cut.method <- match.arg(cut.method, c("best", "ntile"))
  if (cut.method == "best") {
    cutcat <- PlotSurv.groupBestPoint(data, time, status, yitem)
  } else {
    cutcat <- PlotSurv.groupNtile(data, time, status, yitem, n)
  }
  # define legend text
  tb <- table(cutcat[[yitem]])
  tb_lab <- paste0(c("low = ", "high = "), c(tb[1], tb[2]))
  # survival model
  # if (!special_name) {
  #  fit.text <- sprintf('surv_fit(Surv(%s, %s) ~ %s, data = cutcat)',time, status, yitem)
  # } else fit.text <- sprintf("surv_fit(Surv(%s, %s) ~ cutcat[,'%s'], data = cutcat)",time, status, yitem)
  # fit<- eval(parse(text = fit.text))
  fit <- surv_fit(Surv(cutcat[, time], cutcat[, status]) ~ cutcat[, yitem], data = cutcat)

  # annotation
  ano <- surv_pvalue(fit)
  ano <- paste0(ano$method, ":", ano$pval.txt)
  if (addCoxAnotation == T) {
    anoCox <- PlotSurv.addCoxAnotation(data, time, status, yitem)
    ano <- paste0(anoCox, "\n ", ano)
  }
  psurv <- ggsurvplot(fit,
    data = cutcat,
    surv.median.line = "hv",
    # pval = T,pval.size = 8,
    # pval.method = T, pval.method.size = 8,
    legend.title = yitem, legend.labs = tb_lab,
    conf.int = F,
    risk.table = F,
    palette = col.pal,
    legend = "right", ...
  )
  x.position <- quantile(cutcat[, time], na.rm = T)[2]
  x.position <- unname(x.position)
  psurv$plot + annotate("text", x = x.position, y = 0.25, label = ano)
}


## nSurv ----------
PlotSurv.mplt <- function(data, time, status, yitem, yitem2, n = 2,
                          cut.method = c("best", "ntile"), col.pal = "npg", addCoxAnotation = T, signif_num = 2, ...) {
  require(survival)
  require(survminer)
  cut.method <- match.arg(cut.method, c("best", "ntile"))
  if (cut.method == "best") {
    cutcat <- PlotSurv.groupBestPoint2(data, time, status, yitem = yitem, yitem2 = yitem2)
  } else {
    cutcat <- PlotSurv.groupNtile2(data, time, status, yitem, yitem2, n)
  }
  cutcat <- PlotSurv.mergeGroup(cutcat)
  colnames(cutcat)[1:2] <- c("time", "status")
  time <- "time"
  status <- "status"
  # fit
  # fit.text <- sprintf('surv_fit(Surv(%s, %s) ~ group, data = cutcat)',time, status)
  # fit<- eval(parse(text = fit.text))
  fit <- surv_fit(Surv(time, status) ~ group, data = cutcat)

  # annotation
  logrank_sig <- PlotSurv.pairwiseCom(cutcat, time, status)
  cutcat[, "group"] <- as.factor(cutcat[, "group"])
  coxph.fit <- coxph(Surv(time, status) ~ group, data = cutcat)
  coxSummary <- summary(coxph.fit)
  hr <- round(coef(coxSummary)[, 2], 5)
  HR <- c(1, as.vector(unlist(hr)))

  A <- c("Low; ", "High;", "Low; ", "High;")
  B <- c("Low; ", "Low; ", "High;", "High;")
  Group <- c(1, 2, 3, 4)
  groupNum <- table(cutcat$group)
  text.legend1 <- paste0(yitem, " = ", A, yitem2, " = ", B, " Group = ", Group, ", n = ", groupNum, ", HR = ", HR)
  text.legend2 <- paste0(logrank_sig$pvalue_name, " : ", logrank_sig$pvalue)
  ano <- ""
  for (i in text.legend2) {
    ano <- paste(ano, i, sep = "\n ")
  }

  ano2 <- surv_pvalue(fit)
  ano2 <- paste0(ano2$method, ":", ano2$pval.txt)
  if (addCoxAnotation == T) {
    anoCox <- PlotSurv.addCoxAnotation(signif_num, data, time, status, yitem)
    ano2 <- paste0(anoCox, "\n ", ano2)
  }

  # plot
  psurv <- ggsurvplot(fit,
    data = cutcat,
    surv.median.line = "hv",
    # pval = T,pval.size = 8,
    # pval.method = T, pval.method.size = 8,
    legend.title = yitem, legend.labs = text.legend1,
    conf.int = F,
    risk.table = F,
    palette = col.pal,
    legend = c(0.35, 0.15), ...
  )
  x.position <- mean(fivenum(cutcat[, "time"], na.rm = T)[c(3, 4)])

  x.position2 <- mean(fivenum(cutcat[, "time"], na.rm = T)[c(1, 2)])
  psurv$plot + annotate("text", x = x.position, y = 1, label = ano) + annotate("text", x = x.position2, y = 0.4, label = ano2)
}




# example--------
# library(survival)
# library(survminer)
# data(cancer, package = "survival")
# colnames(cancer)
# #
# # # #单组比较
# args1 <- list(data = cancer, time = "time", status = "status", yitem = "wt.loss")
# args2 <- list(cut.method = "best", addCoxAnotation = T)
# args3 <- list(n = 3, cut.method = "ntile", addCoxAnotation = T)
# 
# cutcat1 <- PlotSurv.groupBestPoint(data = cancer, time = "time", status = "status", yitem = "wt.loss")
# cutcat1 <- do.call(PlotSurv.groupBestPoint, args1) # 根据生存差异最显著的阈值点分组
# cutcat2 <- do.call(function(...) PlotSurv.groupNtile(..., n = 3), args1) # 根据中位值或1/n位值分组
# do.call(PlotSurv.cox, args1) # 单因素cox检验
# do.call(PlotSurv.addCoxAnotation, args1) # 向KM plot添加cox注释信息

# do.call(PlotSurv.splt, c(args1, args2)) # 绘制KM曲线，根据最显著阈值分组
# do.call(PlotSurv.splt, c(args1, args3)) # 绘制KM曲线，根据1/n位值分组, n =3
# do.call(function(...)
#   PlotSurv.splt(..., xlab = 'overall survival', title = 'KM-survival curve of wt.loss'),
#   c(args1, args3)) #设置ggsurvplot的其他参数
#
# # # #多组比较
# data(rotterdam, package = 'survival')
# argsm <- list(data = rotterdam, time = 'dtime', status = 'death', yitem = 'age', yitem2 = 'er')
#
# cutcat3 <- do.call(PlotSurv.groupBestPoint2, argsm) #根据生存差异最显著的阈值点分组
# cutcat4 <- do.call(function(...) PlotSurv.groupNtile2(..., n =3), argsm) # 根据中位值或1/n位值分组
# cutcat <- PlotSurv.mergeGroup(cutcat4) #合并两个组别
# pairwiseCom <- PlotSurv.pairwiseCom(cutcat, 'dtime', 'death') #进行两两组别的比较
# do.call(PlotSurv.mplt, c(argsm, list(cut.method = 'ntile'))) #二变量KM生存曲线
#


# # onestep --------------
# library(survival)
# library(survminer)
# data(cancer, package = 'survival')
# head(cancer)
# # 连续变量单因素生存分析
# PlotSurv.splt(data = cancer, time = 'time', status = 'status', yitem = 'wt.loss',
#               cut.method = 'best',  addCoxAnotation = T,  col.pal = c('#1D91C0','#F79EC0')) #cut.method = 'best', 用差异最显著的值分组
#
# PlotSurv.splt(data = cancer, time = 'time', status = 'status', yitem = 'wt.loss',
#               cut.method = 'ntile', n =2, addCoxAnotation = T,  col.pal = c('#1D91C0','#F79EC0')) #cut.method = 'ntile', n =2, 用均值分组
#
# PlotSurv.splt(data = cancer, time = 'time', status = 'status', yitem = 'wt.loss',
#               cut.method = 'ntile', n =3, addCoxAnotation = T,  col.pal = c('#1D91C0','#F79EC0')) #cut.method = 'ntile', n =3, 用25%和75%的值分组
#
# # 连续变量双因素生存分析
# PlotSurv.mplt(data = cancer, time = 'time', status = 'status',yitem = 'wt.loss', yitem2 = 'meal.cal',
#               cut.method = 'ntile', n =2, addCoxAnotation = T)


DataSurv.splt <- function(data, time, status, yitem, n = 2, # special_name = F,
                          cut.method = c("best", "ntile"), col.pal = c("#1D91C0", "#F79EC0"), ...) {
  require(survival)
  require(survminer)
  cut.method <- match.arg(cut.method, c("best", "ntile"))
  if (cut.method == "best") {
    cutcat <- PlotSurv.groupBestPoint(data, time, status, yitem)
  } else {
    cutcat <- PlotSurv.groupNtile(data, time, status, yitem, n)
  }
  # define legend text
  tb <- table(cutcat[[yitem]])
  # tb_lab <- paste0(c("low = ", "high = "), c(tb[1], tb[2]))
  # survival model
  # if (!special_name) {
  #  fit.text <- sprintf('surv_fit(Surv(%s, %s) ~ %s, data = cutcat)',time, status, yitem)
  # } else fit.text <- sprintf("surv_fit(Surv(%s, %s) ~ cutcat[,'%s'], data = cutcat)",time, status, yitem)
  # fit<- eval(parse(text = fit.text))
  fit <- surv_fit(Surv(cutcat[, time], cutcat[, status]) ~ cutcat[, yitem], data = cutcat)
  
  # annotation
  ano <- surv_pvalue(fit)
  # ano <- paste0(ano$method, ":", ano$pval.txt)
  return(c(tb[1], tb[2], ano$pval.txt))

}


cancer <- read.csv('D:/MM86FORXCELL-1.csv', row.names = 1)
os.time <- 'OS'
os.status <- 'Status'
yitems <- setdiff(colnames(cancer), c(os.time, os.status))
low <- c()
high <- c()
p_value <- c()
for (i in yitems) {
  result <- tryCatch(DataSurv.splt(data = cancer, time = os.time, status = os.status,
                                    yitem = i, cut.method = 'best'),
                     error = function(e) {c(NA, NA, NA)}
                     )
  # result <- DataSurv.splt(data = cancer, time = os.time, status = os.status,
  #                        yitem = i, cut.method = 'best')
  low <- append(low, result[1])
  high <- append(high, result[2])
  p_value <- append(p_value, result[3])
}

tb1 <- data.frame(Symbol = yitems, low = low, high = high, p_value = p_value)

write.csv(tb1, 'D:/tb1.csv', row.names = FALSE)

# PlotSurv.splt(data = cancer, time = 'time', status = 'status', yitem = 'wt.loss',
#               cut.method = 'best',  addCoxAnotation = T,  col.pal = c('#1D91C0','#F79EC0')) #cut.method = 'best'