## Function to export results
##
export_results <- function(modresults, Year, modelname,datatype, cells,result_folder, outcomeVar){
  # collapse list of data.frames back to a data.table
  # rename
  modresults <- data.frame(modresults)
  colnames(modresults)[1:4] = c("coef","se", "pvalue","samplesize")
  modresults$samplesize = as.integer(modresults$samplesize)
  save(modresults, file = paste0(result_folder, paste0(Year, "_", outcomeVar,"_",modelname,"_",datatype,"_",cells,"_",Sys.Date(),".RData"))) 
  return(modresults)
}

lambda <- function(p) median(qchisq(p, df=1, lower.tail=FALSE), na.rm=TRUE) / qchisq(0.5, df=1)

publishFormat<-function(res, rounddigit = 3){
  res$lower<-res[,1]-1.96*res[,2]
  res$upper<-res[,1]+1.96*res[,2]
  res$beta <- round(res[,1], rounddigit)
  res$CI <- paste0("(",round(res$lower, rounddigit),", ",round(res$upper, rounddigit),")")
  res$p <- round(res[,3], rounddigit)
  return(res)
}

sigResults <- function(results, psigcut = $psigcut, rounddigit = 3){
  results <- publishFormat(results)
  results <- na.omit(results)
  results$p.FDR<-p.adjust(results$pvalue,"fdr")
  results$qvalue<-qvalue(results$pvalue)$qvalues
  results<-results[results$pvalue<psigcut,]
  results<-results[order(results$pvalue),]
  return(results)
}

statsummary <- function(d,type){
  res = c(min(d),quantile(d,0.25),median(d),mean(d),quantile(d,0.75),max(d),IQR(d),sd(d))
  if(type == "beta")
    names(res) = c("beta_Min","beta_1stQuartile", "beta_Median","beta_Mean","beta_3rdQuartile","beta_Max","beta_IQR","beta_SD")
  if(type == "M")
    names(res) = c("M_Min","M_1stQuartile", "M_Median","M_Mean","M_3rdQuartile","M_Max","M_IQR","M_SD")
  return(res)
}