rm(list=ls())
library(ggplot2)
library(ggsignif)
library(gridExtra)
f = 'stat.tsv'

df = read.table(f,header = T,sep='\t',quote = "",comment.char = "")

df$group = gsub("_.*","",df$sampleid)


ct <- unique(df$new_celltype)
groups <- unique(df$group)
pcg <- c()
df$Percentage = 1
for (i in groups) {
  all <- sum(df$Freq[df$group==i])
  temp <- (df$Freq[df$group==i]/all)*100
  df$Percentage[df$group==i] = temp
}



new_celltype <- c()
group <- c()
Percentage <- c()
sd <- c()
for (j in ct) {
  for (k in groups) {
    cont <- df$new_celltype == j & df$group == k
    zb <- sum(df[cont,'Percentage'])
    bzc <- sd(df[cont,'Percentage'])
    new_celltype <- c(new_celltype,j)
    group <- c(group,k)
    Percentage <- c(Percentage,zb)
    sd <- c(sd,bzc)
  }
}

df_plot <- data.frame(new_celltype=new_celltype,group=group,Percentage=Percentage,sd=sd)




################################################################################################
p_list <- list()
for (celltype in ct) {
  df_sig <- df[df$new_celltype==celltype,]
  df_part <- df_plot[df_plot$new_celltype==celltype,]
  df_part$group = factor(df_part$group,levels = c('Control','A19','S2308'))
  ymax <- max(df_part$Percentage + df_part$sd) 
  size = ymax/20
  p_list[[celltype]] <- ggplot() + 
    geom_bar(data = df_part,mapping = aes(x=group,y=Percentage,fill=group),stat = "identity") + 
    theme_bw() +
    ggtitle(label = celltype) +
    labs(x="") +
    geom_signif(data = df_sig,mapping = aes(x=group,y=Percentage),comparisons = list(c("A19", "Control"),c("A19", "S2308"),c("Control", "S2308")),test = 't.test',y_position = c(ymax+1*size,ymax+2.5*size,ymax+4*size))+
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45, vjust = 0.5, size = 12, face = "bold"),
          axis.title.x = element_text(size = 14, face = "bold"))+
    geom_errorbar(data=df_part,aes(x=group,y=Percentage,ymin = Percentage - sd, ymax = Percentage + sd,color=group),width=0.,size=2)+
    scale_fill_manual(breaks = c('Control','A19','S2308'),values = c('#00468BFF','#ED0000FF','#42b540ff'))+
    scale_color_manual(breaks = c('Control','A19','S2308'),values = c('#00468BFF','#ED0000FF','#42b540ff'))  + 
    guides(fill = 'none',color='none') 
  ggsave(plot = p_list[[celltype]],filename = paste0(celltype,'.png'),width = 9,height = 6)
  ggsave(plot = p_list[[celltype]],filename = paste0(celltype,'.pdf'),width = 9,height = 6)
}


################################################################################################



p <- grid.arrange(p_list[[ct[1]]],
             p_list[[ct[2]]],
             p_list[[ct[3]]],
             p_list[[ct[4]]],
             p_list[[ct[5]]],
             p_list[[ct[6]]],
             p_list[[ct[7]]],
             p_list[[ct[8]]],
             ncol = 8,
             widths = c(8,8,8,8,8,8,8,8)
             )

ggsave(plot = p,filename = paste0('All_celltype','.png'),width = 15,height = 6)
ggsave(plot = p,filename = paste0('All_celltype','.pdf'),width = 15,height = 6)
