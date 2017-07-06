
installPackages <- function(packages, update = FALSE)
{
  installedPackages <- installed.packages()[,"Package"]
  source("http://bioconductor.org/biocLite.R")
  
  if(update)
  {
    biocLite(packages, suppressUpdates = TRUE)
  } else {
    new.packages <- packages[!(packages %in% installedPackages)]
    if(length(new.packages)) 
    {
      biocLite(new.packages, suppressUpdates = TRUE)
    } else {
      message("All required packages have been installed.")
    }
  }
}

list.of.packages <- c("ggplot2", "gridBase", "qualityTools", "XLConnectJars", "XLConnect", "Hmisc", "MASS", "sandwich", "lmtest",
                      "qvalue", "IlluminaHumanMethylationEPICanno.ilm10b2.hg19",
                      "gplots","scales")

installPackages(list.of.packages)
