obj = rd("/gpfs/oe-scrna/further_analysis/scRNA/Mobi/DOE202212736/20240703/3.Keratinocytes/3.1.newcelltype/newcelltype/seurat.h5seurat")



if(any(! rownames(obj@meta.data) == colnames(obj@assays$RNA@data))){
    stop('order error')
}
# 矩阵
exp = obj@assays$RNA@data
exp = normalize_mtx(exp)
names(exp)[1] = 'genes'
wq(exp,'exp')

# cls 文件
sy = 'rosacea'
dz = 'control'
cls = data.frame(col1="",col2="",col3=obj@meta.data$group)
cls$col1[1] = length(obj@meta.data$group)
cls$col1[2] = 2
cls$col1[3] = 1

cls$col2[1] = '#'
cls$col2[2] = sy
cls$col2[3] = dz

cls = t(cls)

write.table(cls,file='exp.cls',sep=' ',quote=F,col.names=F,row.names=F)



