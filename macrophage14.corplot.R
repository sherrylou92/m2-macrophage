#install.packages("corrplot")


library(corrplot)      #引用包
inputFile="corGeneExp.txt"     #输入文件
setwd("D:\\生信视频\\159.macrophage资料\\159.macrophage资料\\14.corrplot - 10846")       #设置工作目录

#读取文件
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
rt=rt[1:31,]     #选取相关性最显著的30个基因进行可视化(其中的1行是巨噬细胞)
rt=t(rt)      #数据转置
M=cor(rt)     #相关性矩阵

#绘制相关性图形
pdf(file="corpot.pdf", width=7, height=7)
corrplot(M,
         method = "circle",      #图形展示的形式
         order = "original",     #基因的排序方式
         type = "upper",         #图形展示的位置
         col=colorRampPalette(c("green", "white", "red"))(50)     #设置颜色,正相关用红色,负相关用绿色表示
         )
dev.off()


