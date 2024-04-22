rm(list=ls())
setwd('D:/desk/remove_pollution')
library(Seurat)
library(celda)
rds='data_ob_v3.rds'

obj <- readRDS(rds)
obj = UpdateSeuratObject(obj)
obj_clean <- decontX(x=obj@assays$RNA@counts)

bl <- obj_clean$contamination > 0.1
cells = as.vector(as.character(obj@meta.data$rawbc[bl]))


new <- subset(x = obj,subset=rawbc %in% cells)

cells <- obj@meta.data$rawbc
wr <-obj_clean$contamination

df = data.frame(barcode=cells,contamination=wr)
write.table(x = df,file = 'contamination_stat.xls',sep = '\t',quote = F,row.names = F,col.names = T)
