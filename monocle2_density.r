rm(list=ls())
library(ggplot2)
library(ggpointdensity)
library(viridis)
#obj = readRDS('pseudotime_results.rds')
df = readRDS('sample.rds')
df = data.frame(t(df))
names(df) = c('x','y')


ggplot() + 
  geom_pointdensity(data=df,mapping=aes(x=x,y=y)) +  # 密度图主函数 
  scale_color_viridis(option='magma') +
  theme_bw()


