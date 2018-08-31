
suppressWarnings(rm(installPackages))

installPackages <- function(packages, update = FALSE)
{
  installedPackages <- installed.packages()[,"Package"]

  if(update)
  {
    message("Reinstalling all required packages.")
    BiocManager::install(packages, ask = FALSE)
  } else {
    new.packages <- packages[!(packages %in% installedPackages)]
    if(length(new.packages)) 
    {
      message(paste0(new.packages, collapse = ", "), " are missing...")
      BiocManager::install(new.packages, ask = FALSE)
    } else {
      message("All required packages have been installed.")
    }
  }
}

## End