#!/gpfs/oe-scrna/zhengfuxing/conda/loupeR/bin/Rscript
#Sys.getenv("GITHUB_PAT")
#Sys.unsetenv("GITHUB_PAT")

#devtools::install_github("Japrin/STARTRAC")
#install.packages("tictoc")

library(Startrac)
library(ggplot2)
library(tictoc)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(tidyverse)
library(sscVis)
library(Seurat)
library(tidyverse)
library(readr)
library(qs)
library(BiocParallel)
library(ComplexHeatmap)
register(MulticoreParam(workers = 8, progressbar = TRUE)) 



rds = 'data_ob_v3.rds'


obj = rd(rds)

dat  <- obj@meta.data



Roe <- calTissueDist(dat,
         byPatient = F,
         colname.cluster = "new_celltype", # 不同细胞亚群
         colname.patient = "sampleid", # 不同样本
         colname.tissue = "group", # 不同组织
         method = "chisq", # "chisq", "fisher", and "freq" 
         min.rowSum = 0) 


library(reshape2)

df = melt(Roe)

colnames(df) = c('new_celltype','group','value')
df$col = ifelse(df$value>=1,'Enrichment','Depletion')

p = ggplot(data=df,mapping=aes(x=new_celltype,y=group,color=col,size=value)) + 
        geom_point() + 
        labs(color=' ',size='Ro/e') + 
        scale_color_manual(values=c('Enrichment'='red','Depletion'='blue')) + 
        theme_bw() + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave(plot=p,file=paste0(ct,'_roe_plot.png'),width=9,height=6)
ggsave(plot=p,file=paste0(ct,'_roe_plot.pdf'),width=9,height=6)


######################################################################

rds = '/gpfs/oe-scrna/further_analysis/scRNA/BD/DZOE2024052059/SH20241009-20241014/2.subcluster/Monocytes/Clustering/Monocytes.rds'

ct = 'Monocytes'
obj = rd(rds)

dat  <- obj@meta.data



Roe <- calTissueDist(dat,
         byPatient = F,
         colname.cluster = "clusters", # 不同细胞亚群
         colname.patient = "sampleid", # 不同样本
         colname.tissue = "group", # 不同组织
         method = "chisq", # "chisq", "fisher", and "freq" 
         min.rowSum = 0) 


library(reshape2)

df = melt(Roe)

colnames(df) = c('clusters','group','value')
df$clusters = factor(df$clusters)
df$col = ifelse(df$value>=1,'Enrichment','Depletion')

p = ggplot(data=df,mapping=aes(x=clusters,y=group,color=col,size=value)) + 
        geom_point() + 
        labs(color=' ',size='Ro/e') + 
        scale_color_manual(values=c('Enrichment'='red','Depletion'='blue')) + 
        theme_bw() + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave(plot=p,file=paste0(ct,'_roe_plot.png'),width=9,height=6)
ggsave(plot=p,file=paste0(ct,'_roe_plot.pdf'),width=9,height=6)