######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")

library(limma)
library(ggpubr)

setwd("C:\\Users\\melon doctor\\Desktop\\13.cliCor")                     #修改工作目录
file="singleGeneExp.txt"                                                      #输入文件
rt=read.table(file,sep="\t",header=T,check.names=F,row.names=1)              #读取表达数据文件
rt$CancerType <- factor(rt$CancerType)
cli=read.table("clinical.txt",sep="\t",header=T,check.names=F,row.names=1)    #读取临床数据文件
gene=colnames(rt)[1]
clinical=colnames(cli)[1]

#对肿瘤类型循环，进行临床相关性分析
outTab=data.frame()
for(i in levels(rt[,"CancerType"])){
  rt1=rt[(rt[,"CancerType"]==i),]
  #交集样品
  data=cbind(rt1,gene=rt1[,gene])
  data=as.matrix(data[,c(gene,"gene")])
  if(nchar(row.names(data)[1])!=nchar(row.names(cli)[1])){
    row.names(data)=gsub(".$","",row.names(data))}
  data=avereps(data)
  sameSample=intersect(row.names(data),row.names(cli))
  sameData=data[sameSample,]
  sameClinical=cli[sameSample,]
  cliExpData=cbind(as.data.frame(sameClinical),sameData)
  if(nrow(cliExpData)==0){next}
  #设置比较组
  
  group=levels(factor(cliExpData$sameClinical))
  comp=combn(group,2)
  my_comparisons=list()
  for(j in 1:ncol(comp)){my_comparisons[[j]]<-comp[,j]}
  #绘制boxplot
  cliExpData$sameClinical <- factor(cliExpData$sameClinical)
  boxplot=ggboxplot(cliExpData, x="sameClinical", y="gene", color="sameClinical",
                    xlab=clinical,
                    ylab=paste(gene,"expression"),
                    legend.title=clinical,
                    title=paste0("Cancer: ",i),
                    add = "jitter")+ 
    stat_compare_means(comparisons = my_comparisons)
  cliExpData <- cbind(rownames(cliExpData),cliExpData)
  colnames(cliExpData)[1] <- "sample"
  write.table(cliExpData,file = paste0(clinical,".",i,".txt"),quote = F,row.names = F,sep="\t",col.names = T)
  pdf(file=paste0(clinical,".",i,".pdf"),width=5.5,height=5)
  print(boxplot)
  dev.off()
}



######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056
