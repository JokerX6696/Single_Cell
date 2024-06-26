#!/home/dongjiaoyang/miniconda3/envs/OESingleCell/bin/Rscript
# Author: Zheng Fuxing
# Date: 2024-01-30
# Description: This script illustrates the visualization of expression levels of different genes across various cells using grouping.
# 装包
library(ggplot2)
library(ggsignif)
library(optparse)
library(stringr)
# 传参
option_list <- list(
  make_option(c('-i',"--input"), type="character", help="input seurat.h5seurat file"),
  make_option(c('-g',"--group"), type="character", help="group ,eg:group:A:B  or group:A:B,group:A:C,group:D:B"),
  make_option(c("--genelist", "-v"),  type="character", help="input genelist ")
)

opt <- parse_args(OptionParser(option_list=option_list))
# para
h5seurat_file = opt$input
genelist_file = opt$genelist
groups = opt$group
colors = c('#00468BFF','#ED0000FF')
# 绘图函数
sig= function(x,y) {
geom_signif(comparisons = x, #指定比较对象
            test = y, #指定检验方法
            size = 0.4, #指定标记中线条的尺寸
            textsize = 2.6, #指定标记中文字部分的大小
            vjust = -0.05, #指定标记中文字部分与横线之间的距离指定标记中文字部分与横线之间的距离
            step_increase = 0.1, #指定每根线条距离高低
            #tip_length = c(0.2, 0.2), #指定短竖线的长度
            map_signif_level =T 
            )
}
# 自定义 labeller 函数
custom_labeller <- function(variable, value) {
  str_wrap(as.character(value), width = 2)  # 根据需要调整宽度
}
# 读取 seurat 对象
obj = OESingleCell::ReadX(input = h5seurat_file, informat = 'h5seurat', verbose = F)
# 将基因列表保存至变量 并与基因簇取交集
genelist = read.table(genelist_file,header=F)
genelist = genelist[,'V1']
sel = grepl(genelist,rownames(obj),ignore.case = TRUE)
# 模糊匹配 忽略大小写！
nomatch=genelist[! sapply(genelist, function(x) any(grepl(x, rownames(obj), ignore.case = TRUE)))]  # save 没有匹配到的基因
write.table(file='gene_no_matched_list.txt',x=nomatch,row.names=F,col.names=F,quote=F) # 输出没有匹配到的基因
vec = c()
for(i in genelist){temp=rownames(obj)[grepl(paste0('^',i,'$'),rownames(obj),ignore.case=T)];vec = c(vec,temp)}
genelist = vec # 模糊匹配 

# 处理 group 输入数据
groups = strsplit(groups,'group:|,',perl=T)[[1]]
groups = groups[which("" != groups)]
##################################################################################################################
if (!file.exists('diff_plot')) {
  # 如果目录不存在，创建目录
  dir.create('diff_plot')
  print(paste("Directory", 'diff_plot', "created."))
} else {
  print(paste("Directory", 'diff_plot', "already exists."))
}
# 提取 基因 组别 以及归一化后的表达量
for(gp in groups){

    g1 = strsplit(gp,':')[[1]][1]
    g2 = strsplit(gp,':')[[1]][2]
    for(gene in genelist){
    # 提取不同组表达量
    group1=as.numeric(obj@assays$RNA@data[gene,which(obj@meta.data$group == g1)])
    ct1=obj@meta.data$new_celltype[which(obj@meta.data$group == g1)]
    df1=data.frame(group=rep(g1,length(group1)),counts=group1,new_celltype=ct1)
    group2=as.numeric(obj@assays$RNA@data[gene,which(obj@meta.data$group == g2)])
    ct2=obj@meta.data$new_celltype[which(obj@meta.data$group == g2)]
    df2=data.frame(group=rep(g2,length(group2)),counts=group2,new_celltype=ct2)
    # 根据需求拼接表格
    df=rbind(df1,df2)  
    df$new_celltype_space = gsub('_',' ',df$new_celltype)
    p=ggplot(df,aes(x=group,y=counts)) +                             
    geom_jitter(aes(color=group),size = 0.5,width=0.4,height = 0) +
    theme_bw() + 
    theme(panel.grid =element_blank()) +
    labs(x = "group", y = gene) + 
    facet_wrap(~ new_celltype_space, ncol = 4, scales = "free_y", labeller = labeller(new_celltype_space = label_wrap_gen(width = 17))) + 
    theme(legend.position = "none") +
    scale_color_manual(breaks=c(g1,g2),values=colors) + 
    sig(list(c(g1,g2)),"wilcox.test") 
    ggsave(p,filename=paste0('diff_plot/',gene,'_',g1,'_vs_',g2,'.png'),device='png',height=length(unique(obj@meta.data$new_celltype))) 
    ggsave(p,filename=paste0('diff_plot/',gene,'_',g1,'_vs_',g2,'.pdf'),device='pdf',height=length(unique(obj@meta.data$new_celltype))) 
}
}



