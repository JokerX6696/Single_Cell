rm(list=ls())
setwd('D:/desk/单细胞文件/pbmc3k/')
library('Seurat')
library(homologene)
data_dir_1 <- 'filtered_gene_bc_matrices/hg19'
data_dir_2 <- 'filtered_gene_bc_matrices/hg38'
pre_vet <- c(data_dir_1,data_dir_2)
names(pre_vet) <- c('hg19','hg38')
library(dplyr)
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
scRNA1=NormalizeData(scRNA1)
scRNA1 <- ScaleData(scRNA1)
pbmc = scRNA1
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
Idents(pbmc)


rm(list=ls())
#############################
a=read.table('D:/desk/单细胞文件/pbmc3k/ratgenelist.txt')
d=c(a[,1][1:32623])
rat2mm=homologene(d, inTax = 10116, outTax = 10090)
setwd
write.table(x=rat2mm[,1:2],file = 'mmgenelist.txt',quote = F,row.names = F,col.names = T)
#################################################
g=c()
n=c()
counts=1
for (i in rat2mm[,1]) {
  
  if(! i %in% g){
    n=c(n,counts)
  }
  g=c(g,i)
  counts = counts + 1
}
newnames=rat2mm[n,2]
mj=rat2mm[n,1]
obj=readRDS('data_ob_v3.rds')

new=subset(obj,features = mj)
rownames(new) <- newnames
