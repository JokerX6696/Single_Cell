obj = rd('clusters_all.h5seurat')
var_gene = obj@assays$RNA@var.features
f='/gpfs/oe-scrna/further_analysis/scRNA/10x/ZOE2023010555/20241010/decontX_DOE202213564/genelist.txt'
df = read.delim(f)
mk_lst = c()
for(i in df){mk_lst = c(mk_lst,i)}
mk_lst = mk_lst[mk_lst != ""]
mk_lst = as.character(CaseMatch(search = mk_lst,match = rownames(obj)))
other_var = var_gene[! var_gene %in% mk_lst]
len_mk_lst = length(mk_lst)
var = c(mk_lst,other_var[1:(200-len_mk_lst)])
obj@assays$RNA@var.features = var
obj <- ScaleData(
  obj,                   # Seurat 对象
  features = var,           # 要进行标准化的基因，默认为所有高变基因
  vars.to.regress = NULL,    # 需要回归掉的协变量（如细胞周期评分、批次效应）
  do.scale = TRUE,           # 是否缩放为均值为0，默认是
  do.center = TRUE           # 是否将数据居中，即减去均值，默认是
)
seurat_object = obj
# pca
seurat_object <- RunPCA(seurat_object, features = seurat_object@assays$RNA@var.features,npcs = 30, do.print = F, verbose = F)
# knn
seurat_object <- FindNeighbors( seurat_object, reduction = 'pca', dims = 1:30,features = seurat_object@assays$RNA@var.features,nn.eps = 0, force.recalc = T, verbose = F)
seurat_object = FindClusters(object = seurat_object, resolution = 0.4,algorithm = 1, verbose = F)
seurat_object <- RunUMAP(seurat_object, dims = 1:30)




sh(seurat_object,prefix='var_200')