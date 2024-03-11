rm(list=ls())
setwd('D:/desk/cmd')
library('Seurat')
library('pheatmap')
gene = as.character(read.table('top20_sampleid_M51-vs-S33_genes.xls',header=T,sep='\t')$gene)
obj = readRDS('data_ob_v3.rds')
################# 做别的参数改这里！
celltype='AT2'
sample1="M51"
sample2="S31"
###########################################
obj = subset(obj,subset=newcelltype==celltype & (sampleid == sample1 | sampleid == sample2))

cell_order = c(rownames(obj@meta.data)[obj@meta.data$sampleid==sample1] ,rownames(obj@meta.data)[obj@meta.data$sampleid==sample2])

mtx = data.frame(obj@assays$RNA@data[gene,cell_order])

tj = table(obj@meta.data$sampleid)
sample1_num = as.numeric(tj[sample1])
sample2_num = as.numeric(tj[sample2])
annotation_col = data.frame(
  identity = factor(rep(c(sample1, sample2), c(sample1_num, sample2_num)))
)
rownames(annotation_col) = colnames(mtx)
ann_colors = list(
  identity = c("#efbd25", "#cd4275")
)
names(ann_colors$identity) = c(sample1,sample2)
png(filename=paste0(celltype,'_',sample1,'_vs_',sample2,'.png'),width=900,height=600)
pheatmap(
    mat = mtx,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    cluster_col = FALSE,
    show_colnames=FALSE,
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    border=FALSE,
    scale = "row"
)
dev.off()

pdf(file=paste0(celltype,'_',sample1,'_vs_',sample2,'.pdf'),width=900,height=600)
pheatmap(
    mat = mtx,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    cluster_col = FALSE,
    show_colnames=FALSE,
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    border=FALSE,
    scale = "row"
)
dev.off()