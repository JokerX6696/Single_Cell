rm(list=ls())
setwd('D:/desk/单细胞文件/pbmc3k/')
library('Seurat')

data_dir <- 'filtered_gene_bc_matrices/hg19'
### 读取 cellranger 生成的文件
counts <- Read10X(data.dir = data_dir)
### 创建 seurat 对象
obj <- CreateSeuratObject(
  counts = counts,
  min.cells=1,
  project = "run_test"
  )


















