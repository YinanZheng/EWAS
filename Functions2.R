### Functions for EWAS Pipeline-2

suppressWarnings(rm(export_results))
suppressWarnings(rm(publishFormat))
suppressWarnings(rm(sigResults))
suppressWarnings(rm(statsummary))

suppressWarnings(rm(f.RLM.Adjusted.Robust.par))

## Function to export results
export_results <- function(modresults, NAMES_LIST, result_folder, rounddigit = rounddigit){
  # collapse list of data.frames back to a data.table
  # rename
  modresults <- data.frame(modresults)
  colnames(modresults)[1:4] = c("Estimates","StdErr", "Stat","Pvalue")
  modresults$Sample_Size = as.integer(modresults$Sample_Size)
  modresults <- publishFormat(modresults, rounddigit = rounddigit)
  saveRDS(modresults, file = paste0(result_folder, 
                                    paste0(NAMES_LIST$cohortname, "_", NAMES_LIST$Year, "_", NAMES_LIST$VAR,"_",NAMES_LIST$modelname,"_",
                                           NAMES_LIST$datatype,"_",NAMES_LIST$cells,"_", NAMES_LIST$nPC,"PC_", NAMES_LIST$tag,"_",Sys.Date(),".RDS"))) 
  message("EWAS results exported!")
  return(modresults)
}

publishFormat<-function(res, rounddigit = 3){
  res$lower<-res[,"Estimates"]-1.96*res[,"StdErr"]
  res$upper<-res[,"Estimates"]+1.96*res[,"StdErr"]
  res$beta <- round(res[,"Estimates"], rounddigit)
  res$CI <- paste0("(",round(res$lower, rounddigit),", ",round(res$upper, rounddigit),")")
  res$p <- round(res[,"Pvalue"], rounddigit)
  return(res)
}

splitAutosomal <- function(res, annot)
{
  cpg_auto <- as.character(annot$Name[!annot$chr %in% c("chrX", "chrY")])
  cpg_X <- as.character(annot$Name[annot$chr %in% c("chrX")])
  cpg_Y <- as.character(annot$Name[annot$chr %in% c("chrY")])
  
  length(cpg_auto)
  length(cpg_X)
  length(cpg_Y)
  
  results_auto <- res[which(rownames(res) %in% cpg_auto),]
  results_X <- res[which(rownames(res) %in% cpg_X),]
  results_Y <- res[which(rownames(res) %in% cpg_Y),]
  
  return(list(auto = results_auto, X = results_X, Y = results_Y))
}

sigResults <- function(results, annotcord, NAMES_LIST, psigcut = psigcut, rounddigit = rounddigit){
  results <- na.omit(results)
  results$p.FDR<-p.adjust(results$Pvalue,"fdr")
  results$qvalue<-qvalue(results$Pvalue)$qvalues
  results<-results[results$Pvalue<psigcut,]
  results<-results[order(results$Pvalue),]
  
  # Add annotation
  results = cbind(results,annotcord[match(rownames(results),annotcord$Name),])
  write.csv(results, paste0(result_folder, paste0(NAMES_LIST$cohortname, "_", NAMES_LIST$Year, "_", NAMES_LIST$VAR,"_", NAMES_LIST$modelname,"_",
                                                  NAMES_LIST$datatype,"_", NAMES_LIST$cells,"_", NAMES_LIST$nPC,"PC_", NAMES_LIST$tag,"_",Sys.Date(),".csv")))
  message("Signficant results exported!")
}

# Add summary of statistics of the tested CpG sites (will add 17 columns)
statsummary <- function(bigdata, type){
  samplesize <- nrow(bigdata)
  if(type == "Mval")
  {
    Mval <- bigdata$methy
    betaVal <- 2^Mval/(2^Mval + 1)
  }

  if(type == "beta")
  {
    betaVal <- bigdata$methy
    Mval <- log2(betaVal/(1-betaVal))
  }
  
  res = c(samplesize, min(betaVal),quantile(betaVal,0.25),median(betaVal),mean(betaVal),quantile(betaVal,0.75),max(betaVal),IQR(betaVal),sd(betaVal),
          min(Mval),quantile(Mval,0.25),median(Mval),mean(Mval),quantile(Mval,0.75),max(Mval),IQR(Mval),sd(Mval))
  names(res) = c("Sample_Size","beta_Min","beta_1stQuartile", "beta_Median","beta_Mean","beta_3rdQuartile","beta_Max","beta_IQR","beta_SD",
                 "M_Min","M_1stQuartile", "M_Median","M_Mean","M_3rdQuartile","M_Max","M_IQR","M_SD")
  return(res)
}

### Modeling functions:

## Debug
# methcol = setNames(seq_len(ncol(tdatRUN)), dimnames(tdatRUN)[[2]])[1]

## RLM
f.RLM.par <- function(methcol, VAR, COV, model_statement, datatype, tdatRUN) { 
  bigdata <- data.frame(na.omit(cbind(VAR = eval(parse(text = paste0("df$", VAR))),methy = tdatRUN[, methcol], COV)))
  mod <- try(rlm(model_statement, bigdata, maxit=200))
  # pull out a data.frame with results
  if("try-error" %in% class(mod)){
    b <- rep(NA, 21)
  } else {
    cf <- try(coeftest(mod, vcov=vcovHC(mod, type="HC0")))
    if(class(cf) == "try-error"){
      b <- rep(NA, 21)
    } else {b <- c(cf[2,], statsummary(bigdata, datatype))
    }
  }
  invisible(b)
}

## LM
f.LM.par <- function(methcol, VAR, COV, model_statement, datatype, tdatRUN) { 
  bigdata <- data.frame(na.omit(cbind(VAR = eval(parse(text = paste0("df$", VAR))),methy = tdatRUN[, methcol], COV)))
  mod <- try(lm(model_statement, bigdata))
  if("try-error" %in% class(mod)){
    b <- rep(NA, 21)
  } else {
    cf <- summary(mod)$coefficients
    b <- c(cf[2,], statsummary(bigdata, datatype))
  }
  invisible(b)
}

## LOGISTIC
f.LOGISTIC.par <- function(methcol, VAR, COV, model_statement, datatype, tdatRUN) { 
  bigdata <- data.frame(na.omit(cbind(VAR = eval(parse(text = paste0("df$", VAR))),methy = tdatRUN[, methcol], COV)))
  mod <- try(glm(model_statement, bigdata, family = binomial))
  if("try-error" %in% class(mod)){
    b <- rep(NA, 21)
  } else {
    cf <- summary(mod)$coefficients
    b <- c(cf[2,], statsummary(bigdata, datatype))
  }
  invisible(b)
}

message("Function2.R loaded!")