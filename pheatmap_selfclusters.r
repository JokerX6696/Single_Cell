set.seed = 1311
library(pheatmap)
obj = readRDS('heatdata.rds')
ret = kmeans(obj, centers = 3, nstart = 25)
df = data.frame(ret$cluster)
df = normalize_mtx(df)
names(df) = c('gene','cluster')
sce = readRDS('/gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2023121656/20240625/2.slingshot/2.2.bycelltype/slingshot_group1/sce.rds')
wq(df)
heatdata = obj[order(df$cluster),]
annotation_col = data.frame(sce$new_celltype)
rownames(annotation_col) = colnames(sce)
names(annotation_col) = 'new_celltype'
annotation_colors = list(
  new_celltype = c(B_cell = "#7fc97f", HSC = "#1b9e77", CLP = "#beaed4")
)
od = rev(rownames(annotation_col)[order(annotation_col$new_celltype)])
heatdata = heatdata[,od]
p = pheatmap( log1p(heatdata),
                                    scale = "row",
                                    cluster_cols = F,
                                    cluster_rows = F,
                                    show_colnames = F,
                                    show_rownames = F,
                                    color = colorRampPalette( c("#406AA8", "white", "#D91216") )(200),
                                    annotation_col = annotation_col,
                                    annotation_colors = annotation_colors,
                                    fontsize_row = 6, fontface="bold" )
                # ggsave(file.path(output_dir, paste0('Dynamic_genes_heatmap_',i,'.png')),plot=p, width = 8, height = 9)
                ggsave(paste0('Dynamic_genes_heatmap','.pdf'),plot=p, width = 8, height = 9,bg="white")
                ggsave(paste0('Dynamic_genes_heatmap','.png'),plot=p, width = 8, height = 9,bg="white",device='png')