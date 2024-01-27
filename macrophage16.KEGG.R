#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")


#引用包
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

pvalueFilter=0.05      #p值过滤条件
qvalueFilter=0.05      #矫正后的p值过滤条件

#定义图形的颜色
colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}

setwd("D:\\生信视频\\159.macrophage资料\\159.macrophage资料\\16.KEGG - 10846")      #设置工作目录
rt=read.table("corSigResult.txt", header=T, sep="\t", check.names=F)     #读取输入文件

#提取共表达基因, 将基因名字转换为基因id
genes=unique(as.vector(rt[,2]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=data.frame(genes, entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        #去除基因id为NA的基因
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#kegg富集分析
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$genes[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
#保存显著富集的结果
write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)

#定义显示通路的数目
showNum=30
if(nrow(KEGG)<showNum){
	showNum=nrow(KEGG)
}

#柱状图
pdf(file="barplot.pdf", width=6, height=3)
barplot(kk, drop=TRUE, showCategory=showNum, label_format=50, color=colorSel)
dev.off()

#气泡图
pdf(file="bubble.pdf", width=10, height=10)
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=50, color=colorSel)
dev.off()

#基因和通路的关系图
pdf(file="cnetplot.pdf", width=8, height=4)
cnet=cnetplot(kk, circular=TRUE, showCategory=5, colorEdge=TRUE)
print(cnet)
dev.off()



