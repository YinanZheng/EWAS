### Functions for EWAS Pipeline-1

suppressWarnings(rm(addlogtransform))
suppressWarnings(rm(loadMyData))
suppressWarnings(rm(outlierRemove))
suppressWarnings(rm(dataQualityPlot))
suppressWarnings(rm(dataQualityPlot_out_RM))
suppressWarnings(rm(unix2dos))
suppressWarnings(rm(dos2unix))

# function to add log-transformed variables (optional)
addlogtransform <- function(pheno, var.list, base = exp(1)){
  if(!'data.frame' %in% class(pheno))
    stop("pheno type data must be in data.frame!")
  for(var in var.list)
  {
    eval(parse(text = paste0("pheno$log_", var, "=log(pheno$",var,", base = base)")))
  }
  return(pheno)
}

# Load custom dataset
loadMyData <- function(MyDataFileName, sheetName = NULL){
  suffix <- tolower(unlist(strsplit(MyDataFileName,"\\."))[2])
  
  if(suffix == "csv")
  {
    data <- read.csv(file.path(data_folder, MyDataFileName))
    warning("You are loading .csv file. Make sure the ID column is not truncated!")
  }
  
  if(suffix == "xls"))
  {
    message("Excel file detected.")
    data <- read_xls(file.path(data_folder, MyDataFileName))
  } 
  
  if(suffix == "xlsx"))
  {
    message("Excel file detected.")
    data <- read_xlsx(file.path(data_folder, MyDataFileName))
  } 
  
  if(suffix == "rds")
  {
    message(".rds file detected.")
    data <- readRDS(file.path(data_folder, MyDataFileName))
  }
  
  if(suffix == "sas7bdat")
  {
    message(".sas7bdat file detected.")
    data <- read.sas7bdat(file.path(data_folder, MyDataFileName))
  }
  
  return(data)
}
  
# Remove outliers using 3*IQR
outlierRemove <- function(df, variable_interest_outlierRemove){
  quantile2575 <- apply(df[, variable_interest_outlierRemove], 2, function(x) quantile(x, probs = c(0.25, 0.75), na.rm = T))
  IQR <- quantile2575[2,] - quantile2575[1,]
  extremeUpper <- quantile2575[2,] + 3 * IQR
  extremeLower <- quantile2575[1,] - 3 * IQR
  ind <- t((t(df[, variable_interest_outlierRemove]) > extremeUpper | t(df[, variable_interest_outlierRemove]) < extremeLower))
  i = 1
  for(var in variable_interest_outlierRemove)
  {
    message(var, ": ", sum(ind[,i], na.rm = TRUE), " outlier(s) detected!")
    i = i + 1
  }
  df[, variable_interest_outlierRemove][ind] <- NA  
  return(df)
}

# QC plot for variable of intrest in raw data
dataQualityPlot <- function(pheno, var.list, Tag, groupVar = NULL, stackratio = 0.8, width = 9, height = 3){
  if(!'data.frame' %in% class(pheno))
    stop("pheno type data must be in data.frame!")
  dir.create(file.path(result_folder, "QC_plot_Raw_data"))
  for(var in var.list){
    if (is.null(groupVar))
    {
      pheno$groupVar = var
      pheno$groupVar = as.factor(pheno$groupVar)
    } else {
      pheno$groupVar = eval(parse(text = paste0("pheno[,",groupVar,"]")))
    }
    png(paste0(result_folder, "/QC_plot_Raw_data/",Tag, "_DataQualityPlot_", var, ".png"), width = width, height = height, unit = "in", res = 400)
    par(mfrow=c(1,3), mar = c(4,4,2,2), oma=c(0,0,0,0))
    
    d = pheno[,var]
    
    eval(parse(text = paste0(var, "=na.omit(d)")))
    eval(parse(text = paste0("hist(", var, ",col='grey')")))
    eval(parse(text = paste0("qqPlot(", var, ")")))
    
    plot.new()              ## suggested by @Josh
    vps <- baseViewports()
    pushViewport(vps$figure) ##   I am in the space of the autocorrelation plot
    vp1 <-plotViewport(c(0,0,0,0)) ## create new vp with margins, you play with this values 
    
    eval(parse(text = paste0("b = ggplot(pheno, aes(x = groupVar, y = ", var, ")) + 
                             geom_boxplot( colour = '#3366FF') +
                             geom_dotplot(binaxis='y', stackdir='center',
                             stackratio=stackratio, dotsize=0.5, fill = 'grey') +
                             ggtitle('Box plot with dots') +
                             theme(plot.title=element_text(face='bold', size=8))
                             ")))
    
    print(b,vp = vp1)        ## suggested by @bpatiste
    dev.off()
  }
}

# QC plot for variable of intrest after outliers removed
dataQualityPlot_out_RM <- function(pheno, var.list,  Tag, groupVar = NULL, stackratio = 0.8, width = 9, height = 3){
  if(!'data.frame' %in% class(pheno))
    stop("pheno type data must be in data.frame!")
  dir.create(file.path(result_folder, "QC_plot_Out_RM"))
  for(var in var.list){
    if (is.null(groupVar))
    {
      pheno$groupVar = var
      pheno$groupVar = as.factor(pheno$groupVar)
    } else {
      pheno$groupVar = eval(parse(text = paste0("pheno[,",groupVar,"]")))
    }
    png(paste0(result_folder, "/QC_plot_Out_RM/", Tag, "_DataQualityPlot_out_RM_", var, ".png"), width = width, height = height, unit = "in", res = 400)
    par(mfrow=c(1,3), mar = c(4,4,2,2), oma=c(0,0,0,0))
    
    d = pheno[,var]
    
    eval(parse(text = paste0(var, "=na.omit(d)")))
    eval(parse(text = paste0("hist(", var, ",col='grey')")))
    eval(parse(text = paste0("qqPlot(", var, ")")))
    
    plot.new()              ## suggested by @Josh
    vps <- baseViewports()
    pushViewport(vps$figure) ##   I am in the space of the autocorrelation plot
    vp1 <-plotViewport(c(0,0,0,0)) ## create new vp with margins, you play with this values 
    
    eval(parse(text = paste0("b = ggplot(pheno, aes(x = groupVar, y = ", var, ")) + 
                             geom_boxplot( colour = '#3366FF') +
                             geom_dotplot(binaxis='y', stackdir='center',
                             stackratio=stackratio, dotsize=0.5, fill = 'grey') +
                             ggtitle('Box plot with dots') +
                             theme(plot.title=element_text(face='bold', size=8))
                             ")))
    print(b,vp = vp1)        ## suggested by @bpatiste
    dev.off()
  }
}

unix2dos <- function(file){
  txt <- readLines(file)
  con <- file(file, open="wb")
  writeLines(txt, con, sep="\r\n")
  close(con)
}

dos2unix <- function(file){
  txt <- readLines(file)
  con <- file(file, open="wb")
  writeLines(txt, con, sep="\n")
  close(con)
}

message("Function1.R loaded!")
