## Subset EWAS results by gene list (separated by "/")

GenesubEWAS <- function(EWAS_Path, Genelist)
{
  analysisName <- substring(EWAS_Path, gregexpr("CARDIA", EWAS_Path)[[1]], gregexpr(".csv", EWAS_Path)[[1]]-1)
  
  EWAS = read.csv(EWAS_Path, row.names = 1)
  Genelist <- unlist(strsplit(Genelist, "/"))
  ### get alias:
  # set up your query genes
  queryGeneNames <- Genelist
  # use sql to get alias table and gene_info table (contains the symbols)
  # first open the database connection
  dbCon <- org.Hs.eg_dbconn()
  # write your SQL query
  sqlQuery <- 'SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;'
  # execute the query on the database
  aliasSymbol <- dbGetQuery(dbCon, sqlQuery)
  # subset to get your results
  
  oldGeneNames.ind = which(is.na(match(queryGeneNames, aliasSymbol[,5])))
  oldGeneNames = queryGeneNames[oldGeneNames.ind ]
  if(length(oldGeneNames) == 0){
    print("All gene symbols are in the latest version!")
  } else {
    queryGeneNames[oldGeneNames.ind] <- aliasSymbol[match(oldGeneNames, aliasSymbol[,2]),5]
  }
  
  result = queryGeneNames
  
  x <- org.Hs.egSYMBOL2EG
  xx <- as.character(x[result])
  
  y <-org.Hs.egREFSEQ
  yy <- unlist(as.list(y[xx]))
  
  
  yy = yy[substring(yy,1,2) %in% c("NM","NR")]
  
  EWAS.sub = EWAS[grep(paste0("\\b",paste(yy,collapse="\\b|\\b"),"\\b"), EWAS$UCSC_RefGene_Accession),]
  
  EWAS.sub = EWAS.sub[order(EWAS.sub$Pvalue),]
  
  write.csv(EWAS.sub, paste0("EWAS_subset_gene_", analysisName, ".csv"))
  # return(list(oldgenelist = Genelist, newgenelist = result, EWAS = EWAS.sub ))
}

message("EWAS results subset function loaded!")
