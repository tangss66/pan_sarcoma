ExtractGenefromEnrich <- function(pathwayCol,enrichdata) {
  pathwayGene <- c()
  for (i in pathwayCol) {
    gene <- enrichdata[pathwayCol, 'geneID']
    gene <- str_split(gene, '/', simplify = T)
    pathwayGene <- append(pathwayGene, gene)
  }
  if (length(which(is.na(pathwayGene))) != 0) pathwayGene <- pathwayGene[-which(is.na(pathwayGene))]
  if (length(which(pathwayGene == '')) != 0) pathwayGene <- pathwayGene[-which(pathwayGene == '')]
  #pathwayGene <- pathwayGene[-which(pathwayGene == ' ')]
  return(unique(pathwayGene))
}
