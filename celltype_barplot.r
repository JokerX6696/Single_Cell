library(Seurat)
library(ggplot2)
library(patchwork)
obj = readRDS('data_ob_v3.rds')


meta = obj@meta.data[,c('new_celltype','group',"new_celltype_col")]

df = as.data.frame(table(meta))

groups = unique(as.character(df$group))
move_num = 15
##########################################################################
df_plot = df[df$group == groups[2],]
df_plot$Freq = df_plot$Freq/sum(df_plot$Freq) * 100
df_plot$new_celltype = factor(df_plot$new_celltype,rev(levels(df_plot$new_celltype)))
fill_all = as.character(c(unique(df_plot$group),df$new_celltype))
fill_col = c('white',"#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02") 
df_plot$num = round(df_plot$Freq,2)
p1 = ggplot() + 
  geom_bar(data=df_plot,mapping=aes(x=new_celltype,y=100,color=group,fill=group),stat = "identity") + 
  scale_color_manual(values='#000000') + 
  scale_fill_manual(breaks = fill_all,values=fill_col) + 
  geom_bar(data=df_plot,mapping=aes(x=new_celltype,y=Freq,fill=new_celltype),stat = "identity") + 
  labs(x="",y="",title = fill_all[1]) + 
  geom_text(data=df_plot,mapping=aes(x=new_celltype,y=Freq+move_num,label=num)) + 
  theme(panel.background = element_blank(),       # 移除面板背景
        plot.background = element_blank(),        # 移除绘图背景
        panel.grid.major = element_blank(),       # 移除主网格线
        panel.grid.minor = element_blank(),       # 移除次网格线
        axis.title.x = element_blank(),           # 移除 x 轴标题
        #axis.title.y = element_blank(),           # 移除 y 轴标题
        axis.text.x = element_blank(),            # 移除 x 轴文本
        #axis.text.y = element_blank(),            # 移除 y 轴文本
        axis.ticks = element_blank(),             # 移除坐标轴刻度线
        legend.position = "none",                 # 移除图例
        plot.title = element_text(hjust = 0.5)
        ) + 
  coord_flip() 
##########################################################################
df_plot = df[df$group == groups[1],]
df_plot$Freq = df_plot$Freq/sum(df_plot$Freq) * 100
df_plot$new_celltype = factor(df_plot$new_celltype,rev(levels(df_plot$new_celltype)))
fill_all = as.character(c(unique(df_plot$group),df$new_celltype))
fill_col = c('white',"#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02") 
df_plot$num = round(df_plot$Freq,2)
p2 = ggplot() + 
  geom_bar(data=df_plot,mapping=aes(x=new_celltype,y=100,color=group,fill=group),stat = "identity") + 
  scale_color_manual(values='#000000') + 
  scale_fill_manual(breaks = fill_all,values=fill_col) + 
  geom_bar(data=df_plot,mapping=aes(x=new_celltype,y=Freq,fill=new_celltype),stat = "identity") + 
  labs(x="",y="",title = fill_all[1]) + 
  geom_text(data=df_plot,mapping=aes(x=new_celltype,y=Freq+move_num,label=num)) + 
  theme(panel.background = element_blank(),       # 移除面板背景
        plot.background = element_blank(),        # 移除绘图背景
        panel.grid.major = element_blank(),       # 移除主网格线
        panel.grid.minor = element_blank(),       # 移除次网格线
        axis.title.x = element_blank(),           # 移除 x 轴标题
        axis.title.y = element_blank(),           # 移除 y 轴标题
        axis.text.x = element_blank(),            # 移除 x 轴文本
        axis.text.y = element_blank(),            # 移除 y 轴文本
        axis.ticks = element_blank(),             # 移除坐标轴刻度线
        legend.position = "none",                 # 移除图例
        plot.title = element_text(hjust = 0.5)
  ) + 
  coord_flip() 
##########################################################################
df_plot = df[df$group == groups[3],]
df_plot$Freq = df_plot$Freq/sum(df_plot$Freq) * 100
df_plot$new_celltype = factor(df_plot$new_celltype,rev(levels(df_plot$new_celltype)))
fill_all = as.character(c(unique(df_plot$group),df$new_celltype))
fill_col = c('white',"#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02") 
df_plot$num = round(df_plot$Freq,2)
p3 = ggplot() + 
  geom_bar(data=df_plot,mapping=aes(x=new_celltype,y=100,color=group,fill=group),stat = "identity") + 
  scale_color_manual(values='#000000') + 
  scale_fill_manual(breaks = fill_all,values=fill_col) + 
  geom_bar(data=df_plot,mapping=aes(x=new_celltype,y=Freq,fill=new_celltype),stat = "identity") + 
  labs(x="",y="",title = fill_all[1]) + 
  geom_text(data=df_plot,mapping=aes(x=new_celltype,y=Freq+move_num,label=num)) + 
  theme(panel.background = element_blank(),       # 移除面板背景
        plot.background = element_blank(),        # 移除绘图背景
        panel.grid.major = element_blank(),       # 移除主网格线
        panel.grid.minor = element_blank(),       # 移除次网格线
        axis.title.x = element_blank(),           # 移除 x 轴标题
        axis.title.y = element_blank(),           # 移除 y 轴标题
        axis.text.x = element_blank(),            # 移除 x 轴文本
        axis.text.y = element_blank(),            # 移除 y 轴文本
        axis.ticks = element_blank(),             # 移除坐标轴刻度线
        legend.position = "none",                 # 移除图例
        plot.title = element_text(hjust = 0.5)
  ) + 
  coord_flip() 
##########################################################################
p = p1 + p2 + p3
ggsave(plot = p,filename = "celltype_barplot.png",width = 15,height = 9,device = 'png')
ggsave(plot = p,filename = "celltype_barplot.pdf",width = 15,height = 9,device = 'pdf')
