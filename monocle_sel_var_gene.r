# 根据流程 拟时序分析 得到的 monocle 对象 提取top 10， 首要排序 qvalue  次要 vst.variance.standardized
pseudotime_results.rds 来源于/gpfs/oe-scrna/liuxuan/Project/MobiNova/DZOE2023092854_Mouse/houxu20231226/6.Monocle2/pseudotime_results.rds
Hypervariable_gene.txt 来源于 monocle 对象 老师要画 top 10 高变， 但是 pvalue 前几十个都是 0，因此排序方式 首选排序 qvalue，次要排序 vst.variance.standardized 得到该表
genelist_top10.txt 来源于 Hypervariable_gene.txt
library('Seurat')
library('dplyr')
obj=readRDS('pseudotime_results.rds')
data=obj@featureData@data
sorted_df=data %>% arrange(qval, desc(vst.variance.standardized))
out = data.frame(gene=s$gene_short_name,sorted_df[,-10])
write.table(x=out,file="Hypervariable_gene.txt",sep='\t',quote=F,row.names=F,col.names=T)
