#install.packages('e1071')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("preprocessCore")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


library("limma")           #引用包
inputFile="symbol 10846.txt"     #表达数据文件
setwd("D:\\生信视频\\159.macrophage资料\\159.macrophage资料\\07.CIBERSORT - 副本")      #设置工作目录

#读取输入文件，并对输入文件整理
rt=read.table(inputFile, header=T, sep="\t", check.names=F)      #读取输入文件
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]

#输出整理好的数据
out=rbind(ID=colnames(data),data)
write.table(out,file="uniq.symbol.txt",sep="\t",quote=F,col.names=F)

#运行CIBERSORT，得到免疫细胞浸润的结果
source("CIBERSORT.R")
results=CIBERSORT("ref.txt", "uniq.symbol.txt", perm=1000)


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱: seqbio@foxmail.com
######光俊老师微信: eduBio

