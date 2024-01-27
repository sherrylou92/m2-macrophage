#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("tidyverse")
#install.packages("ggplot2")
#install.packages("ggpubr")
#install.packages("ggExtra")


#引用包
library(limma)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggExtra)

cellName="Macrophages M2"     #巨噬细胞名称
corFilter=0.3       #相关系数的过滤条件
pFilter=0.01      #相关性检验pvalue的过滤条件
expFile="symbol 10846.txt"      #表达数据文件
immFile="CIBERSORT-Results.txt"     #免疫细胞浸润的结果文件
setwd("D:\\生信视频\\159.macrophage资料\\159.macrophage资料\\12.cor - 10846")      #设置工作目录

#读取免疫细胞浸润结果文件，并对数据进行整理
immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
immune=immune[immune[,"P-value"]<0.05,]
immune=as.matrix(immune[,1:(ncol(immune)-3)])

#去除正常样品
group=sapply(strsplit(row.names(immune),"\\-"), "[", 2)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
immune=t(immune[group==0,])

#读取表达数据文件，并对输入文件进行整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]
data=log2(data+1)

#去除正常样品
group=sapply(strsplit(colnames(data),"\\-"), "[", 2)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
exp=data[,group==0]

#数据合并
sameSample=intersect(colnames(immune), colnames(exp))
immune=immune[cellName,sameSample,drop=F]
exp=exp[,sameSample,drop=F]
data=rbind(immune, exp)

#提取巨噬细胞的含量
x=as.numeric(immune[cellName,])
outTab=data.frame()
#对基因进行循环，进行相关性检验
for(j in rownames(data)){
	if(cellName==j){next}
    #提取基因的表达量
    y=as.numeric(data[j,])
	corT=cor.test(x, y, method = 'spearman')
	cor=corT$estimate
	pvalue=corT$p.value
	outTab=rbind(outTab, cbind(Cell=cellName, Gene=j, cor, pvalue))
	#保存满足条件的基因
	if((abs(cor)>corFilter) & (pvalue<pFilter)){
		#可视化
		df1=as.data.frame(cbind(x,y))
		p1=ggplot(df1, aes(x, y)) + 
			xlab(cellName)+ ylab(paste0(j, " expression"))+
			geom_point()+ geom_smooth(method="lm", formula=y~x) + theme_bw()+
			stat_cor(method = 'spearman', aes(x =x, y =y))
		p2=ggMarginal(p1, type="histogram", xparams=list(fill = "orange"), yparams=list(fill = "blue"))
		pdf(file=paste0("cor.", j, ".pdf"), width=5, height=4.6)
		print(p2)
		dev.off()
	}
}

#对相关性的结果进行过滤, 得到相关性显著的结果
outTab=outTab[order(as.numeric(as.vector(outTab$pvalue))),]
outTab=outTab[abs(as.numeric(outTab$cor))>corFilter & as.numeric(outTab$pvalue)<pFilter,]
write.table(file="corSigResult.txt", outTab, sep="\t", quote=F, row.names=F)

#输出共表达基因的表达量
expOut=data[c(cellName,as.vector(outTab[,2])),]
expOut=cbind(id=row.names(expOut), expOut)
write.table(file="corGeneExp.txt", expOut, sep="\t", quote=F, row.names=F)


