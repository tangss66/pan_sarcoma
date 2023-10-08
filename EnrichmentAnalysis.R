#' @title EnrichmentAnalysis
#' @description  



EnrichmentAnalysis.go <- function(gene, species = 'hsa', fromtype = 'SYMBOL', golevel = 'BP', 
                                  minGSSize = 1, maxGSSize =500, simple = F, sim.cut =0.7, sim.measure = 'Wang') {
  selectDb <- function(species) {
    if (species == 'hsa') {
      library(org.Hs.eg.db)
      return(org.Hs.eg.db)
    }
    if (species == 'mmu') {
      library(org.Mm.eg.db)
      return(org.Mm.eg.db)
    }
    if (species == 'rat') {
      library(org.Rn.eg.db)
      return(org.Rn.eg.db)
    }
  }
  gene=clusterProfiler::bitr(t(gene),fromType= fromtype, toType = "ENTREZID",OrgDb = selectDb(species))$ENTREZID
  ego<-clusterProfiler::enrichGO(gene =gene,
                                 #universe = names(geneList),
                                 OrgDb = selectDb(species),
                                 ont= golevel,
                                 pAdjustMethod = "BH",
                                 pvalueCutoff = 0.2,
                                 qvalueCutoff = 0.05,
                                 minGSSize = minGSSize,
                                 maxGSSize = maxGSSize,
                                 readable = T)
  if (simple == T) ego <- clusterProfiler::simplify(ego, cutoff = sim.cut, measure = sim.measure)
  return(ego)
}


EnrichmentAnalysis.kegg <- function(gene, species = 'hsa', fromtype = 'SYMBOL',  
                                  minGSSize = 1, maxGSSize =500) {
  require(clusterProfiler)
  selectDb <- function(species) {
    if (species == 'hsa') {
      library(org.Hs.eg.db)
      return(org.Hs.eg.db)
    }
    if (species == 'mmu') {
      library(org.Mm.eg.db)
      return(org.Mm.eg.db)
    }
    if (species == 'rat') {
      library(org.Rn.eg.db)
      return(org.Rn.eg.db)
    }
  }
  gene=clusterProfiler::bitr(t(gene),fromType= fromtype, toType = "ENTREZID",OrgDb = selectDb(species))$ENTREZID
  kkD<-enrichKEGG(gene = gene,organism = species, minGSSize = minGSSize, maxGSSize = maxGSSize, use_internal_data = T)
  kkD <- setReadable(kkD, selectDb(species), 'ENTREZID') 
  return(kkD)
}


EnrichmentAnalysis.base <- function(gene, species = 'hsa', fromtype = 'SYMBOL', golevel = 'BP', 
                                    minGSSize = 1, maxGSSize =500, simple = F, sim.cut =0.7, sim.measure = 'Wang') {
  #gene : a list of gene names
  require(clusterProfiler)
  
  selectDb <- function(species) {
    if (species == 'hsa') {
      library(org.Hs.eg.db)
      return(org.Hs.eg.db)
    }
    if (species == 'mmu') {
      library(org.Mm.eg.db)
      return(org.Mm.eg.db)
    }
    if (species == 'rat') {
      library(org.Rn.eg.db)
      return(org.Rn.eg.db)
    }
  }
  gene=clusterProfiler::bitr(t(gene),fromType= fromtype, toType = "ENTREZID",OrgDb = selectDb(species))$ENTREZID
  ego<-clusterProfiler::enrichGO(gene =gene,
                                 #universe = names(geneList),
                                 OrgDb = selectDb(species),
                                 ont= golevel,
                                 pAdjustMethod = "BH",
                                 pvalueCutoff = 0.2,
                                 qvalueCutoff = 0.05,
                                 minGSSize = minGSSize,
                                 maxGSSize = maxGSSize,
                                 readable = T)
  if (simple == T) ego <- clusterProfiler::simplify(ego, cutoff = sim.cut, measure = sim.measure)
  kkD<-enrichKEGG(gene = gene,organism = species, minGSSize = minGSSize, maxGSSize = maxGSSize, use_internal_data = T)
  kkD <- setReadable(kkD, selectDb(species), 'ENTREZID')
  return(list('GO'=ego,'KEGG'=kkD))
}


EnrichmentAnalysis.cut<- function(data, col.name, col.FC='FC', col.p='p.value', cut.FC=2, cut.p=0.05, positive=T,
                                fromtype='SYMBOL', golevel='BP', species = 'hsa', minGSSize = 1, maxGSSize =500,enrich = 'base', simple = F){
  # selectDb <- function(species) {
  #   if (species == 'hsa') {
  #     library(org.Hs.eg.db)
  #     return(org.Hs.eg.db)
  #   }
  #   if (species == 'mmu') {
  #     library(org.Mm.eg.db)
  #     return(org.Mm.eg.db)
  #   }
  # }
  
  data <- data[which(data[,col.p]<=cut.p),]
  {if(positive == T){
    data <- data[which(data[,col.FC] >= cut.FC),]}
    else{data <-data[which(data[,col.FC]<=cut.FC),]}}
  gene <- data[,col.name]
  if (enrich == 'base') {
    EnrichmentAnalysis.base(gene, species = species, fromtype = fromtype, golevel = golevel, 
                          minGSSize = minGSSize, maxGSSize = maxGSSize, simple = simple)
  } else if (enrich == 'kegg') {
    EnrichmentAnalysis.kegg(gene, species = species, fromtype = fromtype, golevel = golevel, 
                            minGSSize = minGSSize, maxGSSize = maxGSSize)
  } else {
    EnrichmentAnalysis.go(gene, species = species, fromtype = fromtype, golevel = golevel, 
                            minGSSize = minGSSize, maxGSSize = maxGSSize, simple = simple)
  }
}


EnrichmentAnalysis.extractGene <- function(en_result, rownames, aslist = T) {
  require(dplyr)
  require(stringr)
  result <- list()
  for (i in rownames) {
    gene <- en_result[i,'geneID'] %>%  str_split('/',simplify = T) %>% as.vector()
    result[[i]] <- gene
  }
  if (aslist == F) result <- unique(unlist(result))
  return(result)
}


EnrichmentAnalysis.extractphosite <- function(en_result, rownames, aslist = T) {
  require(dplyr)
  require(stringr)
  result <- list()
  for (i in rownames) {
    gene <- en_result[i,'phosite'] %>%  str_split(';',simplify = T) %>% as.vector()
    result[[i]] <- gene
  }
  if (aslist == F) result <- unique(unlist(result))
  return(result)
}

EnrichmentAnalysis.phosite <- function(phosite, sep = '/',species = 'hsa', fromtype = 'SYMBOL', golevel = 'BP', 
                                       minGSSize = 1, maxGSSize =500, simple = F, sim.cut =0.7, sim.measure = 'Wang') {
  gene <- str_split(phosite, sep, simplify = T)[, 1]
  site <- str_split(phosite, sep, simplify = T)[, 2]
  data <- data.frame('gene' = gene, 'mod' = site, 'phosite' = phosite)
  en_result <-  EnrichmentAnalysis.base(unique(gene), species = 'hsa', fromtype = 'SYMBOL', golevel = 'BP', 
                                        minGSSize = 1, maxGSSize =500, simple = F, sim.cut =0.7, sim.measure = 'Wang')
  en_result$GO <- en_result$GO@result; en_result$KEGG <- en_result$KEGG@result
  data_go <- en_result$GO; data_kegg <- en_result$KEGG
  data_go$phosite <- NA; data_go$site_count <- NA
  data_kegg$phosite <- NA; data_kegg$site_count <- NA
  for (i in seq(nrow(data_go))) {
    geneid <- data_go$geneID[i]
    geneid <- strsplit(geneid,'/')[[1]]
    site <- data$phosite[which(data$gene %in% geneid)]
    data_go$site_count[i] <- length(site)
    site <- paste(site, collapse = ';')
    data_go$phosite[i] <- site
  }
  for (i in seq(nrow(data_kegg))) {
    geneid <- data_kegg$geneID[i]
    geneid <- strsplit(geneid,'/')[[1]]
    site <- data$phosite[which(data$gene %in% geneid)]
    data_kegg$site_count[i] <- length(site)
    site <- paste(site, collapse = ';')
    data_kegg$phosite[i] <- site
  }
  en_result$GO <- data_go ; en_result$KEGG <- data_kegg
  return(en_result)
}
 


EnrichmentAnalysis.msigdbr <- function(gene, fromtype = 'SYMBOL', species = 'Homo sapiens', category = 'C2', 
                                       subcategory = NULL, pvalue.cut = 1, qvalue.cut = 1) {
  require(msigdbr)
  require(clusterProfiler)
  if (fromtype != 'SYMBOL') gene=clusterProfiler::bitr(t(gene),fromType= fromtype, toType = "SYMBOL",OrgDb = selectDb(species))$ENTREZID
  mdb <- msigdbr(species, category, subcategory)
  geneset <- mdb[,c('gs_name','gene_symbol')]
  
  en_result <- enricher(gene, pvalueCutoff = pvalue.cut, qvalueCutoff = qvalue.cut, TERM2GENE = geneset)
  return(en_result)
}



