CustomCol2 <- function(n){
  my_palette=c(
    "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
    "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
    "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
    "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17")
  return(my_palette[n])
} 
ggumap = DimPlot(object = data_ob, reduction = "umap", pt.size =0.5, group.by="new_celltype") +
    ggplot2::theme( plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::scale_colour_manual( values = CustomCol2(1:35)) + ggplot2::coord_fixed(ratio = 1)
ggplot2::ggsave(file.path("./",paste0("new_celltype",".pdf")), ggumap,width=8)
system(glue::glue("convert  -verbose -density 500 -trim  new_celltype.pdf  -quality 100  -flatten  new_celltype.png"))

library(dplyr)
simplified_meta = data_ob@meta.data %>%
                        dplyr::rename( "Barcode" = "rawbc") %>%
                        dplyr::select(Barcode,sampleid,clusters,group,new_celltype)
                                                
write.table(simplified_meta, quote = F,sep =",",row.names = F,
             "new_celltype.metadata.csv")          # 这个忘记了，还需要一个这个表格