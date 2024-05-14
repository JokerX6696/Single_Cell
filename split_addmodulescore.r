# library
library(Seurat)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(stringr)
## para
x_z = 'new_celltype'
rds='/gpfs/oe-scrna/further_analysis/scRNA/10x/DZQD2024010596/20240513/1.newcelltype_0411/seurat.h5seurat'
plot_data='/gpfs/oe-scrna/further_analysis/scRNA/10x/DZQD2024010596/20240513/4.addmodulescore/all/geneset_visualization/geneset.addmodeulescore_plot.xls'
obj = rd(rds)
info = read.delim(plot_data)
## 处理数据
cells = info$Barcode
info$x = obj@meta.data[match(cells,obj@meta.data$rawbc),x_z]
df = info
names(df)[3] = 'y'


df$group = factor(df$group,levels=c("control", "silica" , "H2"   ,   "tet"   ,  "H2_tet"))
## 计算显著性
my_comparisons <- list(c("control", "silica"), c("silica", "H2","tet","H2_tet"), c("H2", "tet","H2_tet"),c("tet", "H2_tet"))
sig = data.frame(x=unique(as.character(df$x)))
for(g in my_comparisons){
    t = g[1]
    con = g[2]
    add = c()
    for(ct in unique(as.character(df$x))){
        nums1 = df[df$x == ct & df$group == t,'y']
        nums2 = df[df$x == ct & df$group == con,'y']
        ret = wilcox.test(nums1,nums2)
        ret = ret$p.value
        add = c(add,ret)
    }
    sig = data.frame(sig,ret=add)
    names(sig)[length(sig)] = paste(t,con,sep="_")
}

get_sig = function(f,fh){
    all = c()
    for(i in f){
        if(i<0.001){
            ret = paste0(fh,fh,fh)
        }else if(i < 0.01){
            ret =paste0(fh,fh)
        }else if(i < 0.05){
            ret = fh
        }else{
            ret = "-"
        }
        all = c(all,ret)
    }
    return(all)
}
sig_df = sig
sig_df[,2] = get_sig(sig_df[,2],"*")
sig_df[,3] = get_sig(sig_df[,3],"#")
sig_df[,4] = get_sig(sig_df[,4],"e")
sig_df[,5] = get_sig(sig_df[,5],"&")
rownames(sig_df) = sig_df$x
sig_df = sig_df[,2:5]

# plot
p = ggplot(data=df,mapping=aes(x=x,y=y,fill=group)) + 
    geom_violin() + 
    geom_boxplot(width=0.2,position=position_dodge(0.9),outlier.size = 0) +
    labs(x=x_z,y='Score') + 
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5)) 

ymax = max(df$y)

for(ct in rownames(sig_df)){
    print(ct)
    ps = str_flatten(sig_df[ct,]," ")
    p = p + geom_text(x=ct, y=ymax, label=ps)
}

ggsave(plot=p,filename='addmodulescore_split_bygroup.png',width=16,height=6,device='png')
ggsave(plot=p,filename='addmodulescore_split_bygroup.pdf',width=16,height=6,device='pdf')