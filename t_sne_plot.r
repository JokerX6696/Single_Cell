rm(list=ls())
setwd('D:/desk/xtq_plot')
library('Seurat')
library('ggplot2')
library(patchwork)
obj <- readRDS('new_celltype_20230403.rds')
colors <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f", "#4575b4", "#91bfdb")
sx <- obj@meta.data

T_obj <- subset(obj,subset = (sampleid == 'T'))
N_obj <- subset(obj,subset = (sampleid == 'N'))

p1 <- TSNEPlot(object = T_obj,group.by='new_celltype') 
p1 <- p1 + 
  labs(title = 'Cancer tissue') +
  scale_color_manual(breaks = unique(sx$new_celltype),values = colors)+
  theme_bw() + 
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        ) +
  theme(legend.position = "none")


p2 <- TSNEPlot(object = N_obj,group.by='new_celltype') 
p2 <- p2 + 
  labs(title = 'Normal tissue') +
  scale_color_manual(breaks = unique(sx$new_celltype),values = colors)+
  theme_bw() + 
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
  ) 

p = p1 + p2
ggsave(plot = p,filename = 't-sne_plot.pdf',device = 'pdf',width = 12,height = 5)
ggsave(plot = p,filename = 't-sne_plot.png',device = 'png',width = 12,height = 5)

new_obj <- subset(obj,subset = (new_celltype == 'Endothelial cell'))
diffgene <- FindMarkers(object = new_obj,group.by = 'sampleid',ident.1 = 'T',ident.2 = 'N')

write.table(x = diffgene,file = 'Endothelial_cell_diff_T_N.xls',quote = F,row.names = T,col.names = T,sep = '\t')







