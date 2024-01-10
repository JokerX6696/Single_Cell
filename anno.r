get_celltype = function(cluster){
    cluster = as.character(cluster)
switch(cluster,
    '1' = "Naïve CD8 T",
    '2' = "CD8+T -DLG2 (Effector)",
    '3' = "CD8+T-effector",
    '4' = "4-Treg FOXP3",
    '5' = "CD4",
    '6' = "CD8+T-GZMD (Effector)",
    '7' = "CD8+T-Xcl1",
    '8' = "Proliferating T",
    '9' = "Th17",
    '10' = "Naïve CD4 T",
    '11' = "11-Treg FOXP3",
    '12' = "12-Treg FOXP3",
    "Unknown")
}

get_celltype2 = function(cluster){
    cluster = as.character(cluster)
switch(cluster,
    '1' = "CD8",
    '2' = "CD8",
    '3' = "CD8",
    '4' = "CD4",
    '5' = "CD4",
    '6' = "CD8",
    '7' = "CD8",
    '8' = "others",
    '9' = "CD4",
    '10' = "CD4",
    '11' = "CD4",
    '12' = "CD4",
    "Unknown")
}

library('Seurat')
obj = readRDS('data_ob_v3.rds')
yaqun = c()
yaqun2 = c()
clusters = obj@meta.data$clusters
for(ct in clusters){
    temp = get_celltype(ct)
    temp2 = get_celltype2(ct)
    if(temp=='Unknown' | temp2 == 'Unknow'){print('error 意料之外的 clusters');q()}
    yaqun=c(yaqun,temp)
    yaqun2=c(yaqun2,temp2)
}

obj@meta.data$yaqun_clusters = yaqun
obj@meta.data$yaqun_clusters2 = yaqun2