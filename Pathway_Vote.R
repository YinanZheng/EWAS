## Pathway vote

pathwayVote <- function(pathwayResultsWB)
{
  analysisName <- pathwayResultsWB@filename
  analysisName <- substring(analysisName, gregexpr("CARDIA", analysisName)[[1]], gregexpr(".xlsx", analysisName)[[1]]-1)
  
  message("Processing ", analysisName, "...")
  
  wb <- loadWorkbook(paste0("./Pathway_Vote_", analysisName, ".xlsx"), create = TRUE)
  
  sheetnames <- getSheets(pathwayResultsWB)
  
  GOsheets <- sheetnames[grepl("GO", sheetnames)]
  GOcuts <- gsub("GO_Pcut_", "", GOsheets)
  
  KEGGsheets <- sheetnames[grepl("KEGG", sheetnames)]
  KEGGcuts <- gsub("KEGG_Pcut_", "", KEGGsheets)
  
  UNION_GO <- NULL
  UNION_KEGG <- NULL
  
  # Get the union set first
  for(GS in GOsheets)
  {
    temp <- readWorksheet(pathwayResultsWB, sheet = GS)
    temp <- subset(temp, qvalue < 0.05)
    UNION_GO <- rbind(UNION_GO, temp[, c(1,2)])
  }
  UNION_GO <- unique(UNION_GO)
  message("   # of union GO terms: ", nrow(UNION_GO))
  
  for(KS in KEGGsheets)
  {
    temp <- readWorksheet(pathwayResultsWB, sheet = KS)
    temp <- subset(temp, qvalue < 0.05)
    UNION_KEGG <- rbind(UNION_KEGG, temp[, c(1,2)])
  }
  UNION_KEGG <- unique(UNION_KEGG)
  message("   # of union KEGG pathways: ", nrow(UNION_KEGG))
  
  ## Count
  if(nrow(UNION_GO) > 0)
  {
    for(i in seq_len(length(GOsheets)))
    {
      temp <- readWorksheet(pathwayResultsWB, sheet = GOsheets[i])
      temp <- subset(temp, qvalue < 0.05)
      
      eval(parse(text = paste0("UNION_GO$'P", GOcuts[i], "' <- NA")))
      eval(parse(text = paste0("UNION_GO$'P", GOcuts[i], "'[match(temp$ID, UNION_GO$ID)] <- seq_len(nrow(temp))")))
    }
    UNION_GO$Mean_Rank <- rowMeans(UNION_GO[, -c(1,2)], na.rm = TRUE)
    UNION_GO$SD_Rank <- apply(UNION_GO[, -c(1,2)], 1, function(x) sd(x, na.rm = TRUE))
    UNION_GO <- UNION_GO[order(UNION_GO$Mean_Rank),]
    
    createSheet(wb, name = "GO_Vote")
    writeWorksheet(wb, UNION_GO, sheet = "GO_Vote")
  } else {
    message("   No significant GO terms, skipped.")
  }
  
  if(nrow(UNION_KEGG) > 0)
  {
    for(i in seq_len(length(KEGGsheets)))
    {
      temp <- readWorksheet(pathwayResultsWB, sheet = KEGGsheets[i])
      temp <- subset(temp, qvalue < 0.05)
      
      eval(parse(text = paste0("UNION_KEGG$'P", KEGGcuts[i], "' <- NA")))
      eval(parse(text = paste0("UNION_KEGG$'P", KEGGcuts[i], "'[match(temp$ID, UNION_KEGG$ID)] <- seq_len(nrow(temp))")))
    }
    UNION_KEGG$Mean_Rank <- rowMeans(UNION_KEGG[, -c(1,2)], na.rm = TRUE)
    UNION_KEGG$SD_Rank <- apply(UNION_KEGG[, -c(1,2)], 1, function(x) sd(x, na.rm = TRUE))
    UNION_KEGG <- UNION_KEGG[order(UNION_KEGG$Mean_Rank),]
    
    createSheet(wb, name = "KEGG_Vote")
    writeWorksheet(wb, UNION_KEGG, sheet = "KEGG_Vote")
  } else {
    message("   No significant KEGG pathways, skipped.")
  }
  
  saveWorkbook(wb)
  message("DONE!")
}

message("Pathway vote functions loaded!")
