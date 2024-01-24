library('Seurat')
library(ggplot2)
rds="newcelltype.h5seurat"
data_ob = OESingleCell::ReadX(input = rds, informat = 'h5seurat', verbose = F)
out='/gpfs/oe-scrna/zhengfuxing/Project/scRNA/DZOE20230778683_Mouse/20240123/'
cell_type=unique(data_ob@meta.data$new_celltype)
genes=c('Il1b','Il2','Il4','Il6','Tnf','Ifng','Gzmb','Ccl2','Ccl3','Ccl5','Ccl7','Cxcl10','Cxcl12')
for(gene in genes){
    print(gene)
    if(! gene %in% rownames(data_ob)){print('no')} 
    for(cell in cell_type){
        print(cell)
    obj = subset(data_ob,subset=new_celltype == cell)
    p = VlnPlot(obj,features=gene,group.by='group') + scale_fill_manual(values=c('grey','red'))
    ggsave(filename=paste0(gene,'_',cell,'_vlnplot.png'),plot=p,width=9,height=6)
}
}

