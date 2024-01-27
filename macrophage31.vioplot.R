#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("vioplot")


#���ð�
library(limma)
library(vioplot)

immFile="CIBERSORT-Results.txt"      #����ϸ���������ļ�
riskFile="risk.TCGA.txt"             #�����ļ�
pFilter=0.05                         #CIBERSORT�����������
setwd("D:\\������Ƶ\\159.macrophage����\\159.macrophage����\\31.vioplot")     #���ù���Ŀ¼

#��ȡ���߽���ļ����������ݽ�������
immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
immune=immune[immune[,"P-value"]<pFilter,]
immune=as.matrix(immune[,1:(ncol(immune)-0)])



#��ȡ���˵ķ����ļ�
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#���շ��ն���Ʒ����
lowName=row.names(risk[risk[,"Risk"]=="low",])
highName=row.names(risk[risk[,"Risk"]=="high",])
lowImm=intersect(row.names(immune), lowName)
highImm=intersect(row.names(immune), highName)
rt=rbind(immune[lowImm,], immune[highImm,])
lowNum=length(lowImm)
highNum=length(highImm)

#���С����ͼ
outTab=data.frame()
pdf("vioplot.pdf", height=8, width=12)
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(rt))
y=c(1:ncol(rt))
plot(x,y,
     xlim=c(0,63), ylim=c(min(rt), max(rt)+0.02),
     main="", xlab="", ylab="Fraction",
     pch=21,
     col="white",
     xaxt="n")

#��ÿ������ϸ��ѭ��������С����ͼ���ͷ�������ɫ��ʾ���߷����ú�ɫ��ʾ
for(i in 1:ncol(rt)){
	  if(sd(rt[1:lowNum,i])==0){
	    rt[1,i]=0.00001
	  }
	  if(sd(rt[(lowNum+1):(lowNum+highNum),i])==0){
	    rt[(lowNum+1),i]=0.00001
	  }
	  lowData=rt[1:lowNum,i]
	  highData=rt[(lowNum+1):(lowNum+highNum),i]
	  vioplot(lowData,at=3*(i-1),lty=1,add = T,col = 'green')
	  vioplot(highData,at=3*(i-1)+1,lty=1,add = T,col = 'red')
	  wilcoxTest=wilcox.test(lowData, highData)
	  p=wilcoxTest$p.value
	  if(p<pFilter){
	      cellPvalue=cbind(Cell=colnames(rt)[i], pvalue=p)
		  outTab=rbind(outTab,cellPvalue)
	  }
	  mx=max(c(lowData,highData))
	  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
	  text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",sprintf("%.03f",p))), cex=0.8)
}
legend("topright", 
       c("Low risk", "High risk"),
       lwd=3.5, bty="n", cex=1.2,
       col=c("green","red"))
text(seq(1,64,3),-0.05,xpd = NA,labels=colnames(rt),cex = 0.9,srt = 45,pos=2)
dev.off()

#�������ϸ�������ƺͲ����pvalue
write.table(outTab,file="diff.result.txt",sep="\t",row.names=F,quote=F)

