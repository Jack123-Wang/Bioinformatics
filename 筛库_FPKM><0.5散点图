setwd("C:\\Users\\Administrator\\Desktop\\生信分析记录\\FPKM计算\\output（ggplot2）")

#合并基因表达矩阵dat和基因长度gen_len为dat2
dat<-read.table("GSE78144 gene expression matrix.csv",header = T,sep=",")
dat<-dat[,c(1:5)]
dat$aver<-rowMeans(dat[,c(2:5)])
dat<-dat[,c(-2:-5)]
names(dat)<-(c("Gene_Symbol","LIF_Aver"))
gen_len<-read.table("C:\\Users\\Administrator\\Desktop\\生信分析记录\\课题：mESC筛库（YM）\\output\\质量控制\\FPKM图\\output（基因长度）\\All_mmu10gene_len.txt"
                    ,header=T)
names(gen_len)<-(c("Gene_Symbol","Length"))
dat2<-merge(dat,gen_len)
μFLD=1
dat2$effLength <- dat2$Length - μFLD + 1

#计算FPKM\TPM
countToTpm <- function(counts, effLen)
{
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}
countToFpkm <- function(counts, effLen)
{
  N <- sum(counts)
  exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}
dat2$tpm <- with(dat2, countToTpm(LIF_Aver, effLength))
dat2$fpkm <- with(dat2, countToFpkm(LIF_Aver, effLength))

#筛选出前10%和后10%的基因
T50=16830*0.5
D50=16830*0.5
dat2<-dat2[order(dat2$fpkm, decreasing= T), ]
dat2_Top50<-dat2[c(1:T50),]
dat2_Low50<-dat2[c(D50:16830),]
dat_Screen<-read.table("C:\\Users\\Administrator\\Desktop\\生信分析记录\\课题：mESC筛库（YM）\\原始数据\\Input_vs_total.sgrna_summary.xls"
                       ,header=T)
library(dplyr)
dat_Screen1<-select(dat_Screen,Gene,control_mean,treat_mean)
names(dat_Screen1)[1]<-"Gene_Symbol"
dat_screen1T<-merge(dat_Screen1,dat2_Top10,by="Gene_Symbol")
dat_screen1T<-select(dat_screen1T,Gene_Symbol,control_mean,treat_mean)
dat_screen1D<-merge(dat_Screen1,dat2_Low10,by="Gene_Symbol")
dat_screen1D<-select(dat_screen1D,Gene_Symbol,control_mean,treat_mean)

#对top前10%后10%的基因做散点图

library(ggplot2)
p <- ggplot()
#散点图+X/Y轴的名称修改；
p + geom_point(data = test1, mapping = aes(x = control_mean, y = treat_mean),size=1,color="blue")+ 
  labs(title="FPKM less than 0.5", x="control_mean", y="treat_mean") +
  ylim(0,1000)+xlim(0,1000)

p <- ggplot()
#散点图+X/Y轴的名称修改；
p + geom_point(data = test2, mapping = aes(x = control_mean, y = treat_mean),size=1,color="red")+ 
  labs(title="FPKM more than 0.5", x="control_mean", y="treat_mean") +
  ylim(0,1000)+xlim(0,1000)
