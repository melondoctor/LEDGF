######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")

library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

setwd("C:\\Users\\melon doctor\\Desktop\\98panCancer\\25.GSEA")      #设置工作目录
gene="PSIP1"                               #基因名称
gmtFile="c2.cp.kegg.v2022.1.Hs.symbols.gmt"          #基因集文件

#读入gmt文件
gmt=read.gmt(gmtFile)

#获取目录下的所有肿瘤数据文件
files=dir()
files=grep("^symbol.",files,value=T)

for(i in files){
	#读取肿瘤数据文件
	rt=read.table(i,sep="\t",header=T,check.names=F)
	CancerType=gsub("^symbol\\.|\\.txt$","",i)
	
	#如果一个基因占了多行，取均值
	rt=as.matrix(rt)
	rownames(rt)=rt[,1]
	exp=rt[,2:ncol(rt)]
	dimnames=list(rownames(exp),colnames(exp))
	data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
	data=avereps(data)
	#删除正常样品
	group=sapply(strsplit(colnames(data),"\\-"),"[",4)
	group=sapply(strsplit(group,""),"[",1)
	group=gsub("2","1",group)
	data=data[,group==0]
	
	#按目标基因将样品分成高低表达两组
	dataL=data[,(data[gene,]<=median(data[gene,]))]
	dataH=data[,(data[gene,]>median(data[gene,]))]
	meanL=rowMeans(dataL)
	meanH=rowMeans(dataH)
	meanL[meanL<0.00001]=0.00001
	meanH[meanH<0.00001]=0.00001
	logFC=log2(meanH/meanL)
	logFC=sort(logFC,decreasing=T)

    #富集分析
    kk=GSEA(logFC,TERM2GENE=gmt, nPerm=100,pvalueCutoff = 1)
	kkTab=as.data.frame(kk)
	kkTab=kkTab[kkTab$pvalue<0.05,]
	write.table(kkTab,file=paste0("KEGG.",CancerType,".txt"),sep="\t",quote=F,row.names = F)
	
	#输出富集的图形
	termNum=5
	if(nrow(kkTab)>=termNum){
		gseaplot=gseaplot2(kk, row.names(kkTab)[1:termNum],base_size =8)
		pdf(file=paste0("Term.",CancerType,".pdf"),width=8,height=6)
		print(gseaplot)
		dev.off()
	}
}


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056
