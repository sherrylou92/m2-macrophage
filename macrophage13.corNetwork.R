#install.packages("reshape2")
#install.packages("igraph")


#引用包
library(igraph)
library(reshape2)

inputFile="corSigResult.txt"     #输入文件
setwd("D:\\生信视频\\159.macrophage资料\\159.macrophage资料\\13.network - 10846")     #设置工作目录

#读取输入文件
data=read.table(inputFile, header=T, sep="\t", check.names=F)

#定义网络图的节点和边
corweight = data$cor
weight = corweight+abs(min(corweight))+5
d = data.frame(p1=data$Cell, p2=data$Gene, weight=data$cor)
g = graph.data.frame(data, directed=FALSE)

#设置共表达网络的参数
E(g)$color = ifelse(corweight>0,rgb(254/255,67/255,101/255,abs(corweight)),rgb(0/255,0/255,255/255,abs(corweight)))
E(g)$weight = weight       #边的权重
V(g)$size = 8              #节点的大小
V(g)$shape = "circle"      #节点的形状
V(g)$lable.cex = 1.5         #节点字体的大小
V(g)$color = "white"       #节点的颜色

#可视化
pdf(file="network.pdf", width=7, height=6)
layout(matrix(c(1,1,1,0,2,0),byrow=T,nc=3), height=c(7,1), width=c(3,4,3))
par(mar=c(1.5,2,2,2))
vertex.frame.color = NA
plot(g,layout=layout_nicely,vertex.label.cex=V(g)$lable.cex,edge.width = E(g)$weight,edge.arrow.size=0,vertex.label.color="black",vertex.frame.color=vertex.frame.color,edge.color=E(g)$color,vertex.label.cex=V(g)$lable.cex,vertex.label.font=2,vertex.size=V(g)$size,edge.curved=0.4)

#绘制图例
color_legend = c(rgb(254/255,67/255,101/255,seq(1,0,by=-0.01)),rgb(0/255,0/255,255/255,seq(0,1,by=0.01)))
par(mar=c(2,2,1.5,2), xpd=T, cex.axis=1.6, las=1)
barplot(rep(1,length(color_legend)),border = NA, space = 0,ylab="",xlab="",xlim=c(1,length(color_legend)),horiz=FALSE,
        axes = F, col=color_legend,main="")
axis(3,at=seq(1,length(color_legend),length=5),c(1,0.5,0,-0.5,-1),tick=FALSE)
dev.off()


