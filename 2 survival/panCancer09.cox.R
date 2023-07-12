######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("survival")
#install.packages("forestplot")

library(survival)
library(forestplot)

setwd("C:\\Users\\melon doctor\\Desktop\\98panCancer\\09.cox")                       #设置工作目录
rt=read.table("expTime.txt",header=T,sep="\t",check.names=F,row.names=1)    #读取输入文件
rt$futime=rt$futime/365
gene=colnames(rt)[3]
rt$CancerType <- factor(rt$CancerType)

#对肿瘤类型进行循环
outTab=data.frame()
for(i in levels(rt[,"CancerType"])){
  rt1=rt[(rt[,"CancerType"]==i),]
  cox=coxph(Surv(futime, fustat) ~ rt1[,gene], data = rt1)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(cancer=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxP) )
}
write.table(outTab,file="cox.result.txt",sep="\t",row.names=F,quote=F)    #输出基因和p值表格文件


############绘制森林图函数############
#读取输入文件
coxFile="cox.result.txt"
forestFile="forest.pdf"
forestCol="red"
rt=read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
data=as.matrix(rt)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
clrs <- fpColors(box=forestCol,line="blue", summary="royalblue")      #定义颜色
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )   #定义图片文字
pdf(file=forestFile,width = 9,height = 10,onefile = FALSE)
forestplot(tabletext, 
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(50, "mm"),
           xlog=T,
           lwd.ci=4,
           boxsize=0.6,
           xlab="Hazard ratio",
           txt_gp=fpTxtGp(ticks=gpar(cex=1.1),xlab=gpar(cex = 1.25))
)
dev.off()
############绘制森林图函数############



######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056