#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")


#���ð�
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

pvalueFilter=0.05      #pֵ��������
qvalueFilter=0.05      #�������pֵ��������

#����ͼ�ε���ɫ
colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}

setwd("D:\\������Ƶ\\159.macrophage����\\159.macrophage����\\16.KEGG - 10846")      #���ù���Ŀ¼
rt=read.table("corSigResult.txt", header=T, sep="\t", check.names=F)     #��ȡ�����ļ�

#��ȡ���������, ����������ת��Ϊ����id
genes=unique(as.vector(rt[,2]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=data.frame(genes, entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        #ȥ������idΪNA�Ļ���
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#kegg��������
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$genes[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
#�������������Ľ��
write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)

#������ʾͨ·����Ŀ
showNum=30
if(nrow(KEGG)<showNum){
	showNum=nrow(KEGG)
}

#��״ͼ
pdf(file="barplot.pdf", width=6, height=3)
barplot(kk, drop=TRUE, showCategory=showNum, label_format=50, color=colorSel)
dev.off()

#����ͼ
pdf(file="bubble.pdf", width=10, height=10)
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=50, color=colorSel)
dev.off()

#�����ͨ·�Ĺ�ϵͼ
pdf(file="cnetplot.pdf", width=8, height=4)
cnet=cnetplot(kk, circular=TRUE, showCategory=5, colorEdge=TRUE)
print(cnet)
dev.off()


