rm(list=ls())
library(pheatmap)
library(ggplot2)
library(ggplotify)
file <- '1.2.centered_regulon_activity_groupby_design.xls'

df <- read.table(file, sep = '\t',header = T,row.names = 1)




# 定义一个函数来对数据框按行进行归一化
normalize_rows <- function(df) {
  normalized_df <- as.data.frame(t(apply(df, 1, function(row) {
    min_val <- min(row)
    max_val <- max(row)
    normalized_row <- -1 + 2 * (row - min_val) / (max_val - min_val)
    return(normalized_row)
  })))
  colnames(normalized_df) <- colnames(df)
  return(normalized_df)
}
df = normalize_rows(df)

group <- gsub("_.*$","",colnames(df))
annotation_col = data.frame(
  group = group
)
rownames(annotation_col) <- colnames(df)

p <- pheatmap(
  mat = df,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  cluster_col = FALSE,
  annotation_col = annotation_col,
  border=FALSE,
  angle_col=45,
  fontsize = 6
)

p <- as.ggplot(p)
p <- p + theme(plot.background = element_rect(fill = "white"),plot.margin = margin(l = 5, r = 10))

dev.off()
ggsave(device = 'png',plot = p,filename = "1.3.regulon_activity_heatmap.png",width = 9,height = 6)
ggsave(device = 'pdf',plot = p,filename = "1.3.regulon_activity_heatmap.pdf",width = 9,height = 6)


p2 <- pheatmap(
  mat = df,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  cluster_col = FALSE,
  cluster_row = FALSE,
  annotation_col = annotation_col,
  border=FALSE,
  angle_col=45,
  fontsize = 6
)
dev.off()
p2 <- as.ggplot(p2)
p2 <- p2 + theme(plot.background = element_rect(fill = "white"),plot.margin = unit(c(1,1,1,2),'cm'))
#png(filename = '1.3.regulon_activity_heatmap_no_cluster.png',width = 900,height = 600)
#p2

ggsave(device = 'png',plot = p2,filename = "1.3.regulon_activity_heatmap_no_clusters.png",width = 9,height = 6)
ggsave(device = 'pdf',plot = p2,filename = "1.3.regulon_activity_heatmap_no_clusters.pdf",width = 9,height = 6)






