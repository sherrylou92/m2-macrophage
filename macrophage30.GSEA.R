#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")


#引用包
library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

expFile="symbol 10846.txt"         #表达数据文件
riskFile="risk.TCGA.txt"     #风险文件
gmtFile="c2.cp.kegg.symbols.gmt"     #基因集文件
setwd("D:\\生信视频\\159.macrophage资料\\159.macrophage资料\\30.GSEA -扩大p值")     #设置工作目录

#读取文件,并对输入文件进行整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]



#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(risk), colnames(data))
risk=risk[sameSample,]
data=data[,sameSample]

#高低风险比较，得到logFC
dataL=data[,row.names(risk[risk[,"Risk"]=="low",])]
dataH=data[,row.names(risk[risk[,"Risk"]=="high",])]
meanL=rowMeans(dataL)
meanH=rowMeans(dataH)
meanL[meanL<0.00001]=0.00001
meanH[meanH<0.00001]=0.00001
logFC=log2(meanH)-log2(meanL)
logFC=sort(logFC,decreasing=T)
genes=names(logFC)

#读取基因集文件
gmt=read.gmt(gmtFile)

#进行GSEA富集分析
kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff=1, minGSSize=15, maxGSSize=500)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$pvalue<0.1,]
write.table(kkTab,file="GSEA.result.txt",sep="\t",quote=F,row.names = F)


#输出低风险富集的图形
termNum=5     #设置展示通路的数目，展示前5个富集最显著的通路
kkDown=kkTab[kkTab$NES<0,]
if(nrow(kkDown)>=termNum){
	showTerm=row.names(kkDown)[1:termNum]      #获取展示通路的名称
	gseaplot=gseaplot2(kk, showTerm, base_size=8, title="Enriched in low risk group")
	pdf(file="GSEA.lowrisk.pdf", width=7, height=5.5)
	print(gseaplot)
	dev.off()
}


