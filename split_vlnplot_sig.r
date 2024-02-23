#!/home/dongjiaoyang/miniconda3/envs/OESingleCell/bin/R
#.libPaths(c("/home/dongjiaoyang/miniconda3/envs/OESingleCell/lib/R/library","/home/zhengfuxing/R/x86_64-conda-linux-gnu-library/4.0"))
rm(list=ls())
library(Seurat)
library(ggplot2)
library(gghalves)
library(tidyverse)
library(ggpubr)
data_ob = readRDS('data_ob_v3.rds')#OESingleCell::ReadX(input = '../change_sampleid.h5seurat', informat = 'h5seurat', verbose = F)
genelist=read.table('genelist.txt',header=T)[,1]
genelist=genelist[genelist %in% rownames(data_ob)]
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin,
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])

                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })
for(gene in genelist){
  cells = colnames(data_ob)
  value = as.numeric(data_ob@assays$RNA[gene,cells][1,])
  clusters = data_ob@meta.data[cells,'clusters']
  sample = data_ob@meta.data[cells,'sampleid']
  data_new = data.frame(clusters=clusters,cell=cells,data=value,sample=sample)
  
  geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                                draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                                show.legend = NA, inherit.aes = TRUE) {
    layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
          position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
          params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
  }
  
  ggplot(data_new, aes(x = clusters,y = data, fill = sample))+
    geom_split_violin(trim = T,colour="white")+
    geom_point(stat = 'summary',fun=mean,
               position = position_dodge(width = 0.2))+
    scale_fill_manual(values = c("#1ba7b3","#dfb424"))+
    stat_summary(fun.min = function(x){quantile(x)[2]},
                 fun.max = function(x){quantile(x)[4]},
                 geom = 'errorbar',color='black',
                 width=0.01,size=0.5,
                 position = position_dodge(width = 0.2))+
    stat_compare_means(data = data_new, aes(x = clusters,y = data),
                       # 修改显著性标注：
                       symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                        symbols = c("***", "**", "*", "-")),
                       label = "p.signif",
                       label.y = max(data_new$data),
                       hide.ns = F)+
    theme_bw()+
    xlab("clusters")+
    ylab(gene)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          legend.position = "top",
          #legend.key = element_rect(fill = c("#1ba7b3","#dfb424")),
          legend.justification = "right")
  
  ggsave(paste0(gene,"_violin_plot.pdf"), height = 5, width =8)
  ggsave(paste0(gene,"_violin_plot.png"), height = 5, width =8)
}
####################### 分半提琴图 ####################

# data <- read.table("draw.txt", sep = '\t',header = T)

# data_new <- data %>% 
#   pivot_longer(cols = !X, 
#                names_to = "Samples", 
#                values_to = "Values")

# colnames(data_new)[1] <- "Genes"

# data_new$group <- str_split(data_new$Samples, "_", simplify = T)[,2]

# head(data_new)
# data_new$Values <- log10(data_new$Values + 1)


