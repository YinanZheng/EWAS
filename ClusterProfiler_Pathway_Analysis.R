##### ClusterProfiler

message("Expanding java memeory to 1GB...")
while("package:XLConnect" %in% search())
{
  detach("package:XLConnect", unload = TRUE)
}
options(java.parameters = "-Xmx1024m")
library(XLConnect)

entrez2symbol_internal <- function(entrezid)
{
  paste0(getSYMBOL(unlist(strsplit(entrezid, "/")), data='org.Hs.eg'), collapse = "/")
}

RunClusterProfiler <- function(path, result_name_list, EWAS_PcutList)
{
  for(result_name in result_name_list)
  {
    message(result_name)
    
    wb <- loadWorkbook(paste0("./ClusterProfiler_Pathway_", gsub(".RDS", "", result_name), ".xlsx"), create = TRUE)
    res <- readRDS(file.path(path, result_name))
    res <- res[order(res$Pvalue),]
    
    for(Pcut in EWAS_PcutList)
    {
      
      sig_res <- subset(res, Pvalue < Pcut)
      message("  # of significant CpGs at p<", Pcut, ": ", nrow(sig_res))
      
      EntrezID <- getMappedEntrezIDs(sig.cpg = rownames(sig_res),
                                     all.cpg = rownames(res),
                                     array.type = "EPIC")
      
      ego <- enrichGO(gene          = EntrezID$sig.eg,
                      universe      = EntrezID$universe,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      minGSSize = 0,
                      maxGSSize = 999999,
                      pvalueCutoff  = 1,
                      qvalueCutoff = 1,
                      readable      = TRUE)
      
      ego_result <- ego@result
      if(sum(ego_result$qvalue < 0.05) < 10)
      {
        ego_result_sub <- ego_result[1:10,]
      } else {
        ego_result_sub <- subset(ego_result, qvalue < 0.05)
      }
      
      createSheet(wb, name = paste0("GO_Pcut_", Pcut))
      writeWorksheet(wb, ego_result_sub, sheet = paste0("GO_Pcut_", Pcut))
      
      kk <- enrichKEGG(gene         = EntrezID$sig.eg,
                       universe     = EntrezID$universe,
                       organism     = 'hsa',
                       minGSSize = 0,
                       maxGSSize = 99999,
                       pvalueCutoff = 1,
                       qvalueCutoff = 1)
      
      kk_result <- kk@result
      if(sum(kk_result$qvalue < 0.05) < 10)
      {
        kk_result_sub <- kk_result[1:10,]
      } else {
        kk_result_sub <- subset(kk_result, qvalue < 0.05)
      }
      
      
      kk_result_sub$symbol <- sapply(kk_result_sub$geneID, entrez2symbol_internal)
      
      createSheet(wb, name = paste0("KEGG_Pcut_", Pcut))
      writeWorksheet(wb, kk_result_sub, sheet = paste0("KEGG_Pcut_", Pcut))
    }
    saveWorkbook(wb)
  }
}

message("ClusterProfiler functions loaded!")

