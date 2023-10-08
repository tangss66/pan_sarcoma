#' @title TranslationGeneID


TranslationGeneID <- function(data, col.id, targetID, keytype = 'GI', species = 'human', rm.na = F, delID = F){
  #function : gene id translation
  #col_id : the colnames of gene id
  #targetID : the new gene id to be translated to (e.g. SYMBOL )
  #keytype : the formation of current gene id (e.g. ENSEMBL )
  #species : three species(human, mouse, rat) are supported 
  #delna : if delete the row couldn't be matched
  #delID : if delete previous gene id 
  
  require(org.Hs.eg.db)
  require(org.Mm.eg.db)
  require(org.Rn.eg.db)
  require(AnnotationDbi)
  #keytypes(org.Hs.eg.db)
  # [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"     "ENZYME"      
  # [8] "EVIDENCE"     "EVIDENCEALL"  "GENENAME"     "GO"           "GOALL"        "IPI"          "MAP"         
  # [15] "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"         "PROSITE"     
  # [22] "REFSEQ"       "SYMBOL"       "UCSCKG"       "UNIGENE"      "UNIPROT"     
  ids <- data[, col.id]
  
  if (species == 'human') id.list <- AnnotationDbi::select(org.Hs.eg.db, keys=ids, columns=targetID, keytype = keytype)
  if (species == 'mouse') id.list <- AnnotationDbi::select(org.Mm.eg.db, keys=ids, columns=targetID, keytype = keytype)
  if (species == 'rat') id.list <- AnnotationDbi::select(org.Rn.eg.db, keys=ids, columns=targetID, keytype = keytype)
  
  if (rm.na == T){
    data <- merge(id.list, data, by.x=keytype, by.y=col.id)
  } else {
    data <- merge(id.list, data, by.x=keytype, by.y=col.id, all.x=T)
  }
  
  if (delID==T){
    data <- data[, -which(colnames(data) %in% c(col.id))]
  }
  return(data)
}