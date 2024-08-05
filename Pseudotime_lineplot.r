library(ggplot2)
library(ggpubr)

obj = readRDS('/gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2024032280/20240802/1.monocle/pseudotime_results.rds')
score_file = '/gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2024032280/20240802/6.addmodulescore/geneset_visualization/geneset.addmodeulescore.xls'
score_df = read.delim(score_file)
rownames(score_df) = score_df$Barcode
cells = rownames(obj@phenoData@data)
bk = unique(plot_df$new_celltype)
col = unique(plot_df$new_celltype_col)
col = col[match(levels(bk),bk)]
bk = bk[match(levels(bk),bk)]
s = data.frame(obj@reducedDimS[1,])
plot_df = data.frame(
                    cell = cells,
                    new_celltype = obj@phenoData@data$new_celltype,
                    new_celltype_col = obj@phenoData@data$new_celltype_col,
                    Component = s[cells,],
                    Pseudotime=obj@phenoData@data$Pseudotime,
                    Contraction = score_df[cells,'Contraction'],
                    Synthesis = score_df[cells,'Synthesis']
                    )
                    
labelx = (max(plot_df$Component) + min(plot_df$Component))/2
labely = max(plot_df$Contraction) * 1.1
p = ggplot(data=plot_df,mapping=aes(x=Component,y=Contraction,color=new_celltype)) + 
        geom_point(size=0.5) + 
        scale_color_manual(breaks = bk,values=col) + 
        geom_smooth(aes(group = 1),method = "loess",formula = y ~ x,se=F,color = "black") + 
        stat_cor(aes(group = 1),method = "pearson", label.x = labelx, label.y = labely) +
        theme_classic()

ggsave(filename='score2pse_lineplot_Contraction.png',plot=p,device='png',width=9,height=6)
ggsave(filename='score2pse_lineplot_Contraction.pdf',plot=p,device='pdf',width=9,height=6)




labelx = (max(plot_df$Component) + min(plot_df$Component))/2
labely = max(plot_df$Synthesis) * 1.1
p = ggplot(data=plot_df,mapping=aes(x=Component,y=Synthesis,color=new_celltype)) + 
        geom_point(size=0.5) + 
        scale_color_manual(breaks = bk,values=col) + 
        geom_smooth(aes(group = 1),method = "loess",formula = y ~ x,se=F,color = "black") + 
        stat_cor(aes(group = 1),method = "pearson", label.x = labelx, label.y = labely) +
        theme_classic()

ggsave(filename='score2pse_lineplot_Synthesis.png',plot=p,device='png',width=9,height=6)
ggsave(filename='score2pse_lineplot_Synthesis.pdf',plot=p,device='pdf',width=9,height=6)