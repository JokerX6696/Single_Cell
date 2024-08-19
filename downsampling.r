library('dplyr')
library('tibble')
metadata = obj@meta.data 
sub_data = metadata
sampled_cellmeta = sub_data %>% rownames_to_column() %>%
                                group_by( .dots= "new_celltype" ) %>%
                                sample_frac( size = 0.3,replace = F) %>% column_to_rownames()




cells = sampled_cellmeta$rawbc

df = table(obj@meta.data[,c('new_celltype','sampleid')]) <= 20
## 较少的交叉分组 全部保留！
add_cells = c()
for(i in 1:dim(df)[2]){
  for(k in 1:dim(df)[1]){
    info = df[k,i]
    if(info){
    new_celltype = rownames(df)[k]
    sampleid = colnames(df)[i]
    add_cell = rownames(obj@meta.data)[obj@meta.data$new_celltype == new_celltype & obj@meta.data$sampleid == sampleid]
    add_cells = c(add_cells,add_cell)
    }
  }
}
add_cells = obj@meta.data[cells,'rawbc']
cells = unique(c(cells,add_cells))
save_cells = cells
new_obj = subset(obj,subset=rawbc %in% save_cells)
saveRDS(object=new_obj,file='jcy.rds')

df1 = table(obj@meta.data[,c('new_celltype','sampleid')]) <= 20
df2 = table(obj@meta.data[,c('new_celltype','sampleid')]) <= 20 & table(new_obj@meta.data[,c('new_celltype','sampleid')]) <= 20