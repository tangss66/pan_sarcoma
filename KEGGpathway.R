#' @title KEGGpathway
#' @description 绘制KEGG通路图


PlotKeggPath <- function(DEP, pathway.id, out.suffix) {
  entrezid <- AnnotationDbi::select(org.Hs.eg.db, keys = DEP$Gene, columns = 'ENTREZID', keytype = 'SYMBOL')
  entrezid <- entrezid[-duplicated(entrezid$SYMBOL),]
  entrezid <- entrezid[-which(is.na(entrezid$ENTREZID)),]
  
  geneList <-  merge(DEP, entrezid, by.x = 'Gene', by.y = 'SYMBOL')
  
  geneList2 <- geneList$`FC_A/C`; names(geneList2) <- geneList$ENTREZID
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
