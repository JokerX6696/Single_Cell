library(Seurat)
library(ggplot2)
library(reshape2)
rds='/gpfs/oe-scrna/zhengfuxing/Project/scRNA/DZOE2023121146_human/20240310/rds/all_celltype.rds'
obj = readRDS(rds)
obj = subset(obj,subset=new_celltype %in% c('T_NK','Neutrophil') )

gene = 'SRGN'
cells = rownames(obj@meta.data)
# 构建矩阵
df = data.frame(barcode=cells,exp=obj@assays$RNA@data[gene,cells],group=obj@meta.data$group,celltype=obj@meta.data$new_celltype)
p <- ggplot(df,aes(x=celltype,y=exp,fill=group))+
            geom_violin()+
            geom_boxplot(width=0.05,position=position_dodge(0.9))+
            labs(x="",y="Expression", title=gene) +
            theme_classic()+
            theme(plot.title = element_text(hjust = 0.5))

ggsave(p,device='png',filename='violinplot.png')
ggsave(p,device='pdf',filename='violinplot.pdf')