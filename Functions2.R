## Function to export results
##
suppressWarnings(rm(export_results))
suppressWarnings(rm(lambda))
suppressWarnings(rm(publishFormat))
suppressWarnings(rm(sigResults))
suppressWarnings(rm(statsummary))

suppressWarnings(rm(f.RLM.Adjusted.Robust.par))

### Formating function:

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

sigResults <- function(results, psigcut = psigcut, rounddigit = 3){
  results <- publishFormat(results)
  results <- na.omit(results)
  results$p.FDR<-p.adjust(results$pvalue,"fdr")
  results$qvalue<-qvalue(results$pvalue)$qvalues
  results<-results[results$pvalue<psigcut,]
  results<-results[order(results$pvalue),]
  return(results)
}

# Add summary of statistics of the tested CpG sites (will add 17 columns)
statsummary <- function(bigdata){
  samplesize <- nrow(bigdata)
  Mval <- bigdata$methy
  betaVal <- 2^Mval/(2^Mval + 1)
  
  res = c(samplesize, min(betaVal),quantile(betaVal,0.25),median(betaVal),mean(betaVal),quantile(betaVal,0.75),max(betaVal),IQR(betaVal),sd(betaVal),
          min(Mval),quantile(Mval,0.25),median(Mval),mean(Mval),quantile(Mval,0.75),max(Mval),IQR(Mval),sd(Mval))
  names(res) = c("Sample_Size","beta_Min","beta_1stQuartile", "beta_Median","beta_Mean","beta_3rdQuartile","beta_Max","beta_IQR","beta_SD",
                 "M_Min","M_1stQuartile", "M_Median","M_Mean","M_3rdQuartile","M_Max","M_IQR","M_SD")
  return(res)
}

### Modeling functions:

## Debug
methcol = setNames(seq_len(ncol(tdatRUN)), dimnames(tdatRUN)[[2]])[1]
COV = Covariates


## RLM
f.RLM.Adjusted.Robust.par <- function(methcol, VAR, model_statement, tdatRUN) { 
  
  bigdata <- data.frame(na.omit(cbind(VAR = eval(parse(text = paste0("df$", VAR))),methy = tdatRUN[, methcol], COV)))

  mod <- try(rlm(model_statement, bigdata, maxit=200))
  # pull out a data.frame with results
  if("try-error" %in% class(mod)){
    b <- rep(NA, 21)
  } else {
    cf <- try(coeftest(mod, vcov=vcovHC(mod, type="HC0")))
    if(class(cf) == "try-error"){
      b <- rep(NA, 21)
    } else {b <- c(cf[2,], statsummary(bigdata))
    }
  }
  invisible(b)
}

## LM
# TODO

## LOGISTIC
# TODO

message("Function2.R loaded!")