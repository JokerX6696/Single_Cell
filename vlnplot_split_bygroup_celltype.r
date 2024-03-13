library(Seurat)
library(ggplot2)
library(ggpubr)
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
            theme(plot.title = element_text(hjust = 0.5))+
            stat_compare_means(aes(group=group),                       #按分组进行统计检验
                     method = "wilcox.test",
                     paired = F,                             #非配对t检验
                     symnum.args = list(cutpoint=c(0,0.001,0.01,0.05,1),
                                        symbols=c("***","**","*","ns")),
                     label = "p.signif",
                     #label.y = location$Proportion+0.02,      #添加显著性符号的位置
                     size=4.5)

ggsave(p,device='png',filename='violinplot.png')
ggsave(p,device='pdf',filename='violinplot.pdf')