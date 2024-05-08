rm(list=ls())
setwd('D:/desk/diff_plot')
library(ggplot2)
library(tidyverse)
library(ggrepel)
CustomCol2 <- function(n){
  my_palette=c(
    "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
    "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
    "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
    "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17")
  return(my_palette[n])
} 
### 预处理
annos <- list.files(recursive = T,pattern = "*all_diffexp_genes_anno.xls")
groups <- c();for(i in strsplit(annos,'/')){groups=c(groups,i[2])};groups=unique(groups)
for (group in groups) {
  group1 <- annos[grepl(group,annos)]
  group_vs <- strsplit(group1[1],"/")[[1]][2]
  treat <- strsplit(group_vs,"-vs-")[[1]][1]
  con <- strsplit(group_vs,"-vs-")[[1]][2]
  df = data.frame()
  for(group in group1){
    cluster <- strsplit(group,"-")[[1]][1]
    temp_df <- read.delim(group)
    temp_df$cluster <- cluster
    df <- rbind(df,temp_df)
  }
  df$cluster = as.numeric(df$cluster)
  df <- df %>% arrange(cluster)
  df$cluster = as.character(df$cluster)
  #df <- read.delim('sampleid_Torpedo_embryo-vs-Globular_embryo-all_diffexp_genes_anno.xls')
  ## 适配脚本中的列名
  #names(df)[2] = 'p_val_adj'
  names(df)[8] = 'geneID'
  names(df)[3] = 'log2FC'
  # 添加显著性标签：
  df$label <- ifelse(df$q.value < 0.05, "Q_value < 0.05", "Q_value>=0.05")
  
  # 获取每个cluster中表达差异最显著的10个基因；
  
  all_list <- list()
  
  for(clr in unique(df$cluster)){
    k <- paste0('top10sig',clr)
    all_list[[k]] = filter(df, cluster == clr) %>%
      distinct(geneID, .keep_all = T) %>%
      top_n(5, abs(log2FC))
  }
  
  # 将提取所有cluster的Top10基因表格合并：
  top10sig <- data.frame()
  for(i in names(all_list)){
    top10sig <- rbind(top10sig,all_list[[i]])
  }
  
  
  # 新增一列，将Top10的差异基因标记为2，其他的标记为1；
  df$size <- case_when(
    !(df$geneID %in% top10sig$geneID) ~ 1,
    df$geneID %in% top10sig$geneID ~ 2
  )
  
  # 提取非Top10的基因表格；
  dt <- filter(df, size == 1)
  # 绘制每个Cluster Top10以外基因的散点火山图：
  p <- ggplot() +
    geom_jitter(
      data = dt,
      aes(x = cluster, y = log2FC, color = label),
      size = 0.85,
      width = 0.4
    )
  # 根据图p中log2FC区间确定背景柱长度：
  
  mx <- c();for(i in unique(df$cluster)){mx <- c(mx,max(df[df$cluster==i,'log2FC']))}
  mn <- c();for(i in unique(df$cluster)){mn <- c(mn,min(df[df$cluster==i,'log2FC']))}
  dfbar <- data.frame(
    x = unique(df$cluster),   ###
    y = mx 
  )
  dfbar1 <- data.frame(
    x = unique(df$cluster),
    y =mn
  )
  # # 绘制背景柱：
  # p1 <- ggplot() +
  #     geom_col(
  #         data = dfbar,
  #         mapping = aes(x = x, y = y),
  #         fill = "#dcdcdc", alpha = 0.6
  #     ) +
  #     geom_col(
  #         data = dfbar1,
  #         mapping = aes(x = x, y = y),
  #         fill = "#dcdcdc", alpha = 0.6
  #     )
  
  # 把散点火山图叠加到背景柱上：
  p2 <- ggplot() +
    geom_col(
      data = dfbar,
      mapping = aes(x = x, y = y),
      fill = "#dcdcdc", alpha = 0.6
    ) +
    geom_col(
      data = dfbar1,
      mapping = aes(x = x, y = y),
      fill = "#dcdcdc", alpha = 0.6
    ) +
    geom_jitter(
      data = dt,
      aes(x = cluster, y = log2FC, color = label),
      size = 0.85,
      position=position_jitter(seed=1)
    ) +
    geom_jitter(
      data = top10sig,
      aes(x = cluster, y = log2FC, color = label),
      size = 1,
      position=position_jitter(seed=1)
    )
  # 添加X轴的cluster色块标签：
  dfcol <- data.frame(
    x = c(1:length(unique(df$cluster))),
    y = 0,
    label = unique(df$cluster)
  )
  mycol <- CustomCol2(1:length(unique(df$cluster)))
  #mycol <- c("#E64B357F", "#4DBBD57F", "#00A0877F", "#3C54887F", "#F39B7F7F", "#8491B47F", "#91D1C27F")
  p3 <- p2 + geom_tile(
    data = dfcol,
    aes(x = x, y = y),
    height = 0.6,
    color = "black",
    fill = mycol,
    alpha = 0.6,
    show.legend = F
  )
  
  
  # 给每个Cluster差异表达前Top10基因加上标签：
  p4 <- p3 +
    geom_text_repel(
      data = top10sig,
      aes(x = cluster, y = log2FC, label = geneID),
      force = 1.2, position=position_jitter(seed=1) ,
      max.overlaps=50,
      arrow = arrow(
        length = unit(0.008, "npc"),
        type = "open", ends = "last"
      )
    )
  # 散点颜色调整：
  p5 <- p4 +
    scale_color_manual(
      name = NULL,
      values = c("red", "black")
    )
  # 修改X/Y轴标题和添加cluster数字：
  p6 <- p5 +
    labs(x = "Cluster", y = "average logFC") +
    geom_text(
      data = dfcol,
      aes(x = x, y = y, label = label),
      size = 6,
      color = "white"
    )
  # 自定义主题美化：
  p7 <- p6 +
    theme_minimal() +
    theme(
      axis.title = element_text(
        size = 13,
        color = "black",
        face = "bold"
      ),
      axis.line.y = element_line(
        color = "black",
        size = 1.2
      ),
      axis.line.x = element_blank(),
      axis.text.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = "top",
      legend.direction = "vertical",
      legend.justification = c(1, 0),
      legend.text = element_text(size = 15)
    )
  ggsave(plot = p7,filename = paste0(treat,'_vs_',con,'_plot.png'),bg = 'white',device = 'png',height = 6,width = 9)
  ggsave(plot = p7,filename = paste0(treat,'_vs_',con,'_plot.pdf'),bg = 'white',device = 'pdf',height = 6,width = 9)
}
