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

rat=! duplicated(rat2mm[,1])
mm=! duplicated(rat2mm[,2])
chge <- rat2mm[rat & mm,1:2]




obj=readRDS('data_ob_v3.rds')

new=subset(obj,features = chge[,1])
rownames(new) <- newnames
rownames(new) <- chge[,2]
######################################
RenameGenesSeurat_v2 <- function(obj,newnames,gene.use=NULL,de.assay="RNA") {
  print("Run this before integration. It only changes obj@assays$*@counts, @data and @scale.data, @var.features,@reductions$pca@feature.loadings")
  lassays <- Assays(obj)
  assay.use <- obj@reductions$pca@assay.used
  DefaultAssay(obj) <- de.assay
  if (is.null(gene.use)) {
    all_genenames <- rownames(obj)
  }else{
    all_genenames <- gene.use
    obj <- subset(obj,features=gene.use)
  }
  
  order_name <- function(v1,v2,ref){
    v2 <- make.names(v2,unique=T)
    df1 <- data.frame(v1,v2)
    rownames(df1) <- df1$v1
    df1 <- df1[ref,]
    return(df1)
  }
  
  df1 <- order_name(v1=all_genenames,v2=newnames,ref=rownames(obj))
  all_genenames <- df1$v1
  newnames <- df1$v2
  
  if ('SCT' %in% lassays) {
    if ('SCTModel.list' %in%  slotNames(obj@assays$SCT)) {
      obj@assays$SCT@SCTModel.list$model1@feature.attributes <- obj@assays$SCT@SCTModel.list$model1@feature.attributes[all_genenames,]
      rownames(obj@assays$SCT@SCTModel.list$model1@feature.attributes) <- newnames
    }
  }
  change_assay <- function(a1=de.assay,obj,newnames=NULL,all_genenames=NULL){
    RNA <- obj@assays[a1][[1]]
    if (nrow(RNA) == length(newnames)) {
      if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
      if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
      if (length(RNA@var.features)) {
        df1 <- order_name(v1=all_genenames,v2=newnames,ref=RNA@var.features)
        all_genenames1 <- df1$v1
        newnames1 <- df1$v2
        RNA@var.features            <- newnames1
      }
      if (length(RNA@scale.data)){
        df1 <- order_name(v1=all_genenames,v2=newnames,ref=rownames(RNA@scale.data))
        all_genenames1 <- df1$v1
        newnames1 <- df1$v2
        rownames(RNA@scale.data)    <- newnames1
      }
      
    } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
    obj@assays[a1][[1]] <- RNA
    return(obj)
  }
  
  for (a in lassays) {
    DefaultAssay(obj) <- a
    df1 <- order_name(v1=all_genenames,v2=newnames,ref=rownames(obj))
    all_genenames1 <- df1$v1
    newnames1 <- df1$v2
    obj <- change_assay(obj=obj,a1=a,newnames=newnames1,all_genenames=all_genenames1)
  }
  
  hvg <- VariableFeatures(obj,assay=assay.use)
  if (length(obj@reductions$pca)){
    df1 <- order_name(v1=all_genenames,v2=newnames,ref=hvg)
    all_genenames1 <- df1$v1
    newnames1 <- df1$v2
    rownames(obj@reductions$pca@feature.loadings) <- newnames1
  }
  
  return(obj)
}
dim(new)
new_obj <- RenameGenesSeurat_v2(new,newnames = chge[,2],gene.use=NULL,de.assay="RNA")
