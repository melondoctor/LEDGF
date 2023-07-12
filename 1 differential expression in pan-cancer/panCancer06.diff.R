
#install.packages("ggpubr")

library(ggpubr)

setwd("C:\\Users\\melon doctor\\Desktop\\98panCancer\\06.diff")                 #设置工作目录
data=read.table("singleGeneExp.txt",sep="\t",header=T,check.names=F)   #读取箱线图输入文件
gene=colnames(data)[2]
colnames(data)[2]="expression"
#绘制箱型图
p=ggboxplot(data, x="CancerType", y="expression", color = "Type", 
     ylab=paste0(gene," expression"),
     xlab="",
     palette = c("blue","red") )
p=p+rotate_x_text(60)
pdf(file="boxplot.pdf",width=8,height=5)    #输出图片文件
p+stat_compare_means(aes(group=Type),
      method="wilcox.test",
      symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
      label = "p.signif")
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056
