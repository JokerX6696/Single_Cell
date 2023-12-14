rm(list=ls())
setwd('D:/desk/单细胞文件/pbmc3k/')
library('Seurat')

data_dir_1 <- 'filtered_gene_bc_matrices/hg19'
data_dir_2 <- 'filtered_gene_bc_matrices/hg38'
pre_vet <- c(data_dir_1,data_dir_2)
names(pre_vet) <- c('hg19','hg38')
### 读取 cellranger 生成的文件
counts <- Read10X(data.dir = pre_vet)
### 创建 seurat 对象   #################### 参数值得研究
scRNA1 <- CreateSeuratObject(
  counts = counts,
  min.cells=1,
  project = "run_test"
  )
head(counts[1:6,1:60])
scRNA1 = CreateSeuratObject(counts, min.cells=1)#依据表达矩阵创建seurat对象
scRNA1 = subset(scRNA1, downsample=1000,seed=42)#对每个样本随机采样1000个细胞，设置随机种子便于后续重复
dim(scRNA1)
table(scRNA1@meta.data$orig.ident)
#############################

data_dir_1 <- 'filtered_gene_bc_matrices/hg19'
data_dir_2 <- 'filtered_gene_bc_matrices/hg38'
pre_vet <- c(data_dir_1,data_dir_2)
names(pre_vet) <- c('hg19','hg38')
scRNA2list <- list()
for (i in 1:length(pre_vet)) {
  counts <- Read10X(data.dir = pre_vet[i])
  scRNA2list[[i]] <- CreateSeuratObject(counts, min.cells=1)
}
scRNA2 <- merge(scRNA2list[[1]], y=scRNA2list[[2]])
scRNA2 = subset(scRNA2, downsample=1000,seed=42)
table(scRNA2@meta.data$orig.ident)
head(scRNA2@meta.data)
sum(row.names(scRNA1@meta.data) == row.names(scRNA2@meta.data)) #检验两个object之间的细胞是否一致










