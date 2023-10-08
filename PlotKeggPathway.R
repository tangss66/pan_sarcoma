#' @title  plot KEGG pathway

library(pathview)
library(org.Hs.eg.db)

load("D:/Projects/Pan-sarcoma/Analysis/272samples/proteinCluster_DEP.RData")
geneList <- DEP2
PlotKeggPath <- function(DEP, pathway.id, out.suffix) {
  entrezid <- AnnotationDbi::select(org.Hs.eg.db, keys = DEP$Gene, columns = 'ENTREZID', keytype = 'SYMBOL')
  entrezid <- entrezid[-duplicated(entrezid$SYMBOL),]
  entrezid <- entrezid[-which(is.na(entrezid$ENTREZID)),]
  
  geneList <-  merge(DEP, entrezid, by.x = 'Gene', by.y = 'SYMBOL')
  
  geneList2 <- geneList$`FC_PC/ctrl`; names(geneList2) <- geneList$ENTREZID
  geneList2 <- log2(geneList2)
  keggpath <- pathview(gene.data  = geneList2,
                       pathway.id = pathway.id,
                       species    = "hsa",
                       out.suffix = out.suffix,
                       #kegg.native = FALSE
                       #limit      = list(gene=max(abs(geneList)), cpd=1)
  )
  return(NA)
}

PlotKeggPath3 <- function(pathway.id, suffix){
  PlotKeggPath(DEP1, pathway.id, paste('DEP1', suffix, sep = '.') )
  PlotKeggPath(DEP2, pathway.id, paste('DEP2', suffix, sep = '.') )
  PlotKeggPath(DEP3, pathway.id, paste('DEP3', suffix, sep = '.') )
}

PlotKeggPath3('hsa04110', 'cell_cycle')