rm(list=ls())
setwd('D:/desk/wkdir/scRNA')
library(dplyr)
library(Seurat)
library(patchwork)
# 查看文件，
list.files('pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19')
# 数据导入
pbmc.counts <- Read10X(data.dir = "pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19")
# 创建Seurat对象
pbmc <- CreateSeuratObject(counts = pbmc.counts)
pbmc
str(pbmc)
# 计算线粒体read的百分比
pbmc <- CreateSeuratObject(counts = pbmc.counts, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# 显示前5个细胞的质控指标
head(pbmc@meta.data, 5)
# 查看特征与特征间的相关性
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
p1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
p2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(p1, p2))
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

#
# 识别前2000个特征
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# 识别前10的高异质性基因
top10 <- head(VariableFeatures(pbmc), 10)

# 绘图看看
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))


# 缩放数据
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
head(pbmc[["RNA"]]@scale.data,5)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
# 绘图
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)


# 5.确定数据集的维度
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
# 绘图看看
JackStrawPlot(pbmc, dims = 1:15)

# 聚类细胞
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
# 查看前5聚类
head(Idents(pbmc), 5)

# 使用UMAP聚类
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")
# 显示在聚类标签
DimPlot(pbmc, reduction = "umap", label = TRUE)

# 使用TSNE聚类
pbmc <- RunTSNE(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "tsne")
# 显示在聚类标签
DimPlot(pbmc, reduction = "tsne", label = TRUE)

