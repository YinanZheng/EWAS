
suppressWarnings(rm(installPackages))

installPackages <- function(packages, update = FALSE)
{
  installedPackages <- installed.packages()[,"Package"]
  source("http://bioconductor.org/biocLite.R")
  
  if(update)
  {
    message("Reinstalling all required packages.")
    biocLite(packages, suppressUpdates = TRUE)
  } else {
    new.packages <- packages[!(packages %in% installedPackages)]
    if(length(new.packages)) 
    {
      message(paste0(new.packages, collapse = ", "), " are missing...")
      biocLite(new.packages, suppressUpdates = TRUE)
    } else {
      message("All required packages have been installed.")
    }
  }
}

## End