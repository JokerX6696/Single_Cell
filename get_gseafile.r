obj = rd("/gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2023092903/20250107/1.featureplot/rds/data_ob_v3.rds")


get_gsea_config = function(object,case,control){
    if(any(! rownames(object@meta.data) == colnames(object@assays$RNA@data))){
        stop('order error')
    }
    # 矩阵
    exp = object@assays$RNA@data
    exp = normalize_mtx(exp)
    names(exp)[1] = 'genes'
    outname = paste(case,control,'exp',sep='_')
    wq(exp,prefix=outname)

    # cls 文件
    sy = case
    dz = control
    cls = data.frame(col1="",col2="",col3=object@meta.data$group)
    cls$col1[1] = length(object@meta.data$group)
    cls$col1[2] = 2
    cls$col1[3] = 1

    cls$col2[1] = '#'
    cls$col2[2] = sy
    cls$col2[3] = dz

    cls = t(cls)
    outname = paste(case,control,'exp.cls',sep='_')
    write.table(cls,file='exp.cls',sep=' ',quote=F,col.names=F,row.names=F)
}

