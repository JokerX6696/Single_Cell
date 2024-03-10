rm(list=ls())
library(Seurat)
library(ggplot2)
library(reshape2)
library(ggpubr)
#### para
# rds
# genelist



####
obj <- readRDS('data_ob_v3.rds')
gene <- read.table(file = 'genelist.txt',header=T)$gene
# Slc6a5 没有
gene <- paste(toupper(substring(gene, 1, 1)),
              tolower(substring(gene, 2)),
              sep = "")
gene <- gene[ gene %in% rownames(obj) ]

df <- data.frame(obj@assays$RNA@data[gene,])
df$gene=rownames(df)

df = melt(df,variable.name = 'gene')

names(df)[2] = 'barcode'

df$Group <- obj@meta.data[df$barcode,'group']

p <- ggplot(data = df, mapping = aes(x=gene,y=value,color=Group)) + 
  geom_boxplot() +
  scale_color_manual(values = c('blue','red')) + 
  theme_bw() +
  labs(x = "Gene", y = "Gene expression") +
  theme(legend.position = "top",
        axis.title = element_text(family = "serif"),
        axis.text.x = element_text(angle = 45, hjust = 1,size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
        ) + 
  stat_compare_means(aes(group=Group),                       #按分组进行统计检验
                     method = "wilcox.test",
                     paired = F,                             #非配对wilcox检验
                     symnum.args = list(cutpoint=c(0,0.001,0.01,0.05,1),
                                        symbols=c("***","**","*","")),
                     label = "p.signif",
                     size=4.5,
                     show.legend = FALSE
                     )
ggsave(filename = 'group_expression_boxplot.pdf',plot = p,device = 'pdf',width = 11,height = 6)
ggsave(filename = 'group_expression_boxplot.png',plot = p,device = 'png',width = 11,height = 6)  
