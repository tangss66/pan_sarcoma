genesetH <- c(geneset_cys, geneset_cc, geneset_apo) %>% unique()

clinH <- clin %>% arrange(EBV_infection)
dataH <- ExpData[genesetH, clinH[['T']] ];dataH <- dataH[-c(2,6,10),]

dataH2 <- apply(dataH, 1, scale) %>% t(); colnames(dataH2) <- colnames(dataH)

ha_col <- c('#4dbbd5ff','#e64b35ff'); names(ha_col) <- c('negative','positive') 
ha <- columnAnnotation(EBV_infection = clinH$EBV_infection,
                       #ki_67 = class$ki67.pos.ratio, 
                       col =list(EBV_infection = ha_col),
                       na_col = 'white')
# col_fun <- circlize::colorRamp2(
#   c(-2,0,2),
#   c('blue','white','red')
# )
# col_fun <-circlize::colorRamp2(
#   c(-2, -0.1, 0.1, 2), 
#   c('#323896','#daf0f6', '#f5fbd2', '#a70226'))

col_fun <- circlize::colorRamp2(
  c(-4,0,4),
  c('#395F90','white','#FF0000')
)
col_fun <- circlize::colorRamp2(c(-2,0,2),c('#3993ff','white','#ec5058'))
Heatmap(dataH2, cluster_columns = F, cluster_rows = F,show_column_names = F,
        col = col_fun, top_annotation = ha)



##
library(ComplexHeatmap)
col_fun <- circlize::colorRamp2(c(-3,-0.5,0,0.5,3), c('#4575B4','#CCE6F0','#FAFDC7','#FEE598','#D73027'))

annotation_col = data.frame(
  CellType = factor(rep(c("CT1", "CT2"), 5)), 
  Time = 1:5
)
# dataH2 <- apply(dataH, 1, scale) %>% t()
# colnames(dataH2) <- colnames(dataH)
# Heatmap(dataH2, cluster_rows = F, cluster_columns = F, col = col_fun)

pheatmap(dataH, scale = 'row', cluster_rows = F, cluster_cols = F, color = col_fun, 
         cellwidth = 16, cellheight = 8, gaps_col = c(1,2,3,4,5,6),
         main = 'hierarchical cluster gene set enrichment heatmap', annotation_col = annotation_col
)