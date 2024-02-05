rm(list=ls())
setwd('D:/desk/cellselect')
library(Seurat)

obj = readRDS('clusters13.rds')


plot <- DimPlot(obj, reduction = "umap") 
select.cells <- CellSelector(plot = plot) 
select.cells2 <- CellSelector(plot = plot) 
select.cells3 <- CellSelector(plot = plot)


sel1 = select.cells
sel2 = c(select.cells2,select.cells3)

bc1=c()
for(i in sel1){
  ret = obj@meta.data[i,'rawbc']
  bc1=c(bc1,ret)
}

bc2=c()
for(i in sel2){
  ret = obj@meta.data[i,'rawbc']
  bc2=c(bc2,ret)
}

write.table(bc1,file = '13_1.txt',quote = F,row.names = F,col.names = F)
write.table(bc2,file = '13_2.txt',quote = F,row.names = F,col.names = F)
