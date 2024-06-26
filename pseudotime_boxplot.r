library(monocle)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(combinat)
obj = readRDS('pseudotime_results.rds')
lev = c('WT','asxl1he_SRSF2','asxl1ho_SRSF2')
cells = rownames(obj@phenoData@data)
pseudotime = obj@phenoData@data$Pseudotime
ct = as.character(obj@phenoData@data$new_celltype)
sampleid = as.character(obj@phenoData@data$sampleid)
group = as.character(obj@phenoData@data$sampleid)
group[grepl('wt$',group)] = 'WT'
df = data.frame(cell=cells,pseudotime=pseudotime,celltype=ct,sample=sampleid,group=group)

df$group = factor(df$group,levels=lev)
rownames(df) = cells

for(ct in unique(df$celltype)){
plot_df = df[df$celltype == ct,]
write.table(x=plot_df,file=paste0(ct,'_plotdata.xls'),sep='\t',row.names=F,col.names=T,quote=F)

vs=combn(unique(plot_df$group), 2, simplify = FALSE)
wil_df = data.frame()
for(i in vs){
    ko = i[1]
    wt = i[2]
    x = plot_df$pseudotime[plot_df$group == ko]
    y = plot_df$pseudotime[plot_df$group == wt]
    ret = wilcox.test(x=x,y=y)$p.value
    temp_df = data.frame(group1=ko,group2=wt,wilcox_test=ret)
    wil_df = rbind(wil_df,temp_df)
}
write.table(x=wil_df,file=paste0(ct,'_wilcox_test.xls'),sep='\t',row.names=F,col.names=T,quote=F)
p = ggplot(data=plot_df,mapping=aes(x=group,y=pseudotime,fill=group)) + 
        geom_violin() + 
        geom_boxplot(width=0.2,outlier.shape = NA) + 
        scale_fill_manual(breaks=lev,values=c("grey","#FA250F","#00BED8")) + 
        theme_bw() +
        theme(
    panel.grid.major = element_blank() , # 删除主网格线
    panel.grid.minor = element_blank()   # 删除次网格线
  ) 
        #stat_compare_means(method = "wilcox.test", label = "p.format",comparisons = vs)    
ggsave(plot=p,filename=paste0(ct,'_vlnplot.png'),device='png',width=9,height=6)
ggsave(plot=p,filename=paste0(ct,'_vlnplot.pdf'),device='pdf',width=9,height=6)
}

