### Functions for EWAS Pipeline-3

suppressWarnings(rm(lambda))
suppressWarnings(rm(myqqplot))
suppressWarnings(rm(ewas_diagPlot))
suppressWarnings(rm(heatmap_function))

## Lambda
lambda <- function(p) median(qchisq(p, df=1, lower.tail=FALSE), na.rm=TRUE) / qchisq(0.5, df=1)

# Function for QQ plot
myqqplot <- function(pvector, p0 = -8, col=c("#A0A0A0", "#000000"),showCI = T, ...) {
  p_order <- order(pvector,decreasing=FALSE)
  if (any(pvector == 0)) {
    pvector[pvector == 0] <- .Machine$double.xmin
  }
  o <- -log10(pvector[p_order])
  n <- length(o)
  e <- -log10 (( 1:n - 0.5)/n )
  b <- o >= p0;
  
  plot(e[!b],o[!b],pch=19,cex=0.7, ...,
       xlab=expression(Expected~~-log[10](italic(p))),
       ylab=expression(Observed~~-log[10](italic(p))),
       xlim=c(0,e[1]), ylim=c(0,o[1]),col=col[1])
  ## plot the 95% confidence interval (pointwise)
  if(showCI){
    ## the jth order statistic from a
    ## uniform(0,1) sample
    ## has a beta(j,n-j+1) distribution
    ## (Casella & Berger, 2002,
    ## 2nd edition, pg 230, Duxbury)
    c97.5 <- qbeta(0.975,1:n,n-1:n+1)
    c02.5 <- qbeta(0.025,1:n,n-1:n+1)
    polygon(c(e, rev(e)), -log10(c(c97.5, rev(c02.5))), density=NA, col="gray90")
  }
  points(e[b],o[b],pch=19,cex=0.7,col=col[2])
  abline(a=0,b=1,col=rgb(1,0.65,0),lty=1)
}

# # legend bottom-margin P-value adjust
# my.legend = function(...) {
#   opar <- par(fig=c(0,1,0,1), oma=c(0,0,0,0),
#               mar=c(0,0,0,0), new=TRUE)
#   on.exit(par(opar))
#   plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
#   legend(...)
# }

ewas_diagPlot <- function(modresults, NAMES_LIST, width = 7, height = 7){
  png(paste0(result_folder, "DiagnosticPlot_", NAMES_LIST$cohortname, "_", NAMES_LIST$Year, "_", 
             NAMES_LIST$VAR, "_", NAMES_LIST$modelname, "_", NAMES_LIST$datatype, "_", NAMES_LIST$cells, "_", NAMES_LIST$nPC, "PC_", NAMES_LIST$tag, ".png"), width = width, height = height, unit = "in", res = 300)
  
  par(mfrow=c(2,2), mar = c(5,5,2,2), oma=c(0,0,0,0))
  
  ##############################################
  ###      Histograms of sample size
  ##############################################
  hist(modresults$Sample_Size,
       main="Sample Size Histogram",
       xlab="Sample Size", col="grey")
  
  # max_sample <- max(modresults$Sample_Size)
  # modresults <- modresults[modresults$Sample_Size==max_sample, ]
 
  pval <- modresults$Pvalue
  coef <- modresults$Estimates
  nCpG <- length(pval)
  
  ##############################################
  ###      qq-plot of P-values
  ##############################################
  lambd=round(lambda(pval),3)
  myqqplot(pval, main="QQ-plot EWAS")
  legend("topleft",legend=paste0("lambda = ",lambd),bty="n")
  
  ##############################################
  ###      Histograms of P-values
  ##############################################
  hist(pval,
       main="P-value Histogram",
       xlab="P-value", col="grey")
  
  
  ##############################################
  ###      Volcano Plot 
  ##############################################
  
  ## Uncomment if you need a single volcano plot
  # tiff("VolcanoPlot.tiff", width = 5, height = 5, units = 'in', res = 600, compression = "lzw")
  
  # Color Points with gradient based on P-val
  log.pvalues <- -log10(pval)
  bonThresh = -log10(0.05/nCpG)
  
  if(max(log.pvalues) < bonThresh) color.log.pvalues = c(log.pvalues, bonThresh) else color.log.pvalues = log.pvalues
  
  colorGradient <- colorRampPalette(c("grey","orange"), alpha=F)
  colors = colorGradient(30)[as.numeric(cut(color.log.pvalues,breaks = 30))]
  
  xrg=range(coef)
  yrg=range(color.log.pvalues)
  
  #BonFerroni Cutoff
  # -log10(0.05/nCpG)
  
  qval = qvalue(pval)
  
  if(sum(qval$qvalues<=0.05) == 0) {cat("None of the results has FDR < 0.05\n"); qThresh = 10} else qThresh = max(pval[qval$qvalues<=0.05])
  
  plot(coef,log.pvalues, col=alpha(colors,0.60),
       xlim=xrg, ylim=yrg, pch=16,main="Volcano Plot",
       xlab=expression("Reg. Coeff (M-value)"~(beta)), 
       ylab=expression("-log"[10]~"(P-value)"))
  abline(h=-log10(0.05/nCpG), lty=1, col="black", lwd=2)
  abline(h=-log10(qThresh), lty=2, col="grey", lwd=2)
  legend("topright",legend=c("Bonferroni", "FDR"),lty=c(1,2),
         lwd=c(2,2),col=c("black", "dark grey"),
         horiz=FALSE, bty='n')
  
  # my.legend("topleft", c("Bonferroni", "FDR"), lty=c(1,2),
  #           lwd=c(2,2),col=c("black", "dark grey"),
  #           horiz=FALSE, bty='n')
  
  dev.off()
}

heatmap_function <- function(m_sub, modresultsCSV, NAMES_LIST, HM_CONFIG, df_var){
  CpG_top <- rownames(modresultsCSV)[1:NtopCpG]
  beta <- m_sub[CpG_top, ]
  df_var <- df_var[match(colnames(beta), df_var$ID),]
  stopifnot(identical(colnames(beta), df_var$ID))
  dat <- df_var[,NAMES_LIST$VAR]
  print(table(dat))
  col.factor <- factor(dat, levels=c("High", "Low"))
  
  na.ind <- which(is.na(col.factor))
  if(length(na.ind)>0){
    col.factor <- col.factor[-na.ind]
    beta <- beta[, -na.ind]
  }
  
  message("Distant method: ", HM_CONFIG$dist_method)
  if(HM_CONFIG$dist_method == "pearson")  
  {
    data.dist <- hyperSpec::pearson.dist(beta)
    data.dist.g <- hyperSpec::pearson.dist(t(beta))
  } else {
    data.dist <- dist(beta, method = HM_CONFIG$dist_method)
    data.dist.g <- dist(t(beta), method = HM_CONFIG$dist_method)
  }

  message("Cluster method: ", HM_CONFIG$clust_method)
  ## Row cluster
  row.clus <- hclust(data.dist, method = HM_CONFIG$clust_method)
  
  ## Column cluster
  col.clus <- hclust(data.dist.g,  method = HM_CONFIG$clust_method)
  
  ColSideColors= factor(c("orange", "blue"))[col.factor]
  colors = c(seq(-4,-0.535353535,length=100),seq(-0.5,0.5,length=100),seq(0.535353535,4,length=100))
  my_palette <- colorRampPalette(c("red", "black", "green"))(n = 299)
  
  lab.factor <- col.factor
  png(file = paste0(result_folder, "Heatmap_", NAMES_LIST$cohortname, "_", NAMES_LIST$Year, "_", 
                    NAMES_LIST$VAR, "_", NAMES_LIST$modelname, "_", NAMES_LIST$datatype, "_", NAMES_LIST$cells, "_", NAMES_LIST$nPC, "PC_", NAMES_LIST$tag, 
                    "_CpGtop", NAMES_LIST$NtopCpG,"_", HM_CONFIG$dist_method, "_", HM_CONFIG$clust_method,".png"), 
      width = HM_CONFIG$hmWidth, height = HM_CONFIG$hmHeight, res=300, unit="in")
  heatmap.2(as.matrix(beta), 
            main = paste0(NAMES_LIST$VAR),
            Rowv = as.dendrogram(row.clus),
            Colv = as.dendrogram(col.clus),
            dendrogram="both",
            scale="row",
            col=my_palette,
            cex.main = 2,
            cexRow = 1.4,
            cexCol =1,
            trace = "none", 
            density.info = "none", 
            xlab = "", ylab = "", cex.lab=2,
            ColSideColors = as.character(ColSideColors),
            labRow=rownames(beta),
            labCol=colnames(beta),
            margins = c(10, 8),
            lhei = c(2, 10)
  )
  legend("topright", levels(lab.factor), xpd = TRUE, horiz = TRUE,
         title = "",
         inset = c(-0.01,-0.04), bty = "n",
         fill=c("orange", "blue"), col = c("orange", "blue"),
         border = F,
         cex = 1)
  dev.off()
}

message("Function3.R loaded!")