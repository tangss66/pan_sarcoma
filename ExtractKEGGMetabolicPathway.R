#' @title ExtractKEGGMetabolicPathway
#' @description extract metabolic pathways from kegg database

library(KEGGREST)

#设置代理
library(httr)
set_config(
  use_proxy(url="127.0.0.1", port=19180)
)

org <- keggList('organism')
head(org) 
hsa_path <- keggLink("pathway","hsa")

meta= unique(hsa_path)[grepl('hsa00',unique(hsa_path))]
hsa_info <- lapply(meta, keggGet)
nm=unlist(lapply( hsa_info , function(x) x[[1]]$NAME))

library(stringr)
genes = unlist(lapply( hsa_info , function(x) {
  g = x[[1]]$GENE
  paste(str_split(g[seq(2,length(g),by=2)],';',simplify = T)[,1],collapse =';')
}))
df =data.frame(
  hsa= meta,
  nm=nm,
  genes =genes
)


library(tidyr)
df2 <- df %>% separate_rows('genes', sep = ";")

write.csv(df, 'D:/Projects/Pan-sarcoma/Analysis/R_function/result_ExtractKEGGMetabolicPathway.csv')
write.csv(df2, 'D:/Projects/Pan-sarcoma/Analysis/R_function/result_ExtractKEGGMetabolicPathway2.csv')



