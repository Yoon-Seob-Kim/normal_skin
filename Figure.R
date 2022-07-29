#load the libraries 
library(copynumber)
library(dplyr)
library(plyr)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(data.table)
library(MutationalPatterns)
library(deconstructSigs)
library(BSgenome.Hsapiens.UCSC.hg19)
library(RColorBrewer)
library(SigProfilerExtractorR)
library(reticulate)
library("dndscv")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(pheatmap)
library(plyr)
library(scales)

##Count
setwd("/data/Delta_data4/kysbbubbu/normal_skin/08_Counts")
dat.n=read.delim("Skin.txt",sep="\t",header=T)
head(dat.n)
write.csv(as.data.frame(table(dat.n$Sample)),"count.csv")
write.csv(as.data.frame(table(dat.n$Sample,dat.n$UV.signature)),"UV_count.csv")
write.csv(as.data.frame(table(dat.n$Sample,dat.n$Type)),"function_count.csv")
write.csv(as.data.frame(table(dat.n$Sample,dat.n$Type1)),"type.csv")

###Figure 1. Mutation count and VAFs
###Figure 1. Mutation count and VAFs
###Figure 1. Mutation count and VAFs
##Mutation counts
setwd("/data/Delta_data4/kysbbubbu/normal_skin/08_Counts")
dat.n=read.delim("Skin.txt",sep="\t",header=T)
table(dat.n$Type.1)
head(dat.n)
dat.n$Type=plyr::mapvalues(x=dat.n$Type, from =c("Noncoding","Frameshift","Inframe","Splicing","Nonsense","Missense","Silent"), to =c("Noncoding","Indels","Indels","Splicing","Nonsense","Missense","Silent"))
dat.n2=as.data.frame(table(dat.n$Sample,dat.n$Type))
dat.n3=as.data.frame(table(dat.n$Sample))
dat.n2$Var1=factor(dat.n2$Var1,levels=dat.n2[order(dat.n3$Freq),]$Var1)
dat.n2$Var2=factor(dat.n2$Var2,levels=c("Noncoding","Indels","Splicing","Nonsense","Missense","Silent"))
#Figure1b bar plot 
ggplot(dat.n2, aes(x=Var1, y=Freq, fill=Var2)) +geom_bar(stat = "identity", width=0.7) + theme_classic() + theme(axis.title.x = element_blank(),axis.text.x= element_text(size = 13,face="bold",color="black",angle=90,vjust=0.5),axis.text.y= element_text(size = 13,face="bold",color="black"), axis.title.y = element_text(size=16, face="bold"))+theme(legend.title =element_blank(),legend.text = element_text(size=13), legend.position="right")+scale_fill_manual(values=brewer.pal(6,"Set1"))+ylab("Mutation counts")+ scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+ theme(legend.position = "none")

?scale_y_continuous
##age and sun-exposure state 
setwd("/data/Delta_data4/kysbbubbu/normal_skin/08_Counts")
dat.n=read.delim("Skin_annotation.txt",sep="\t",header=T)
head(dat.n)
dat.n2=as.data.frame(t(dat.n$Age))
rownames(dat.n2)=c("Age")
colnames(dat.n2)=dat.n$Sample
head(dat.n2)
dat.n2=as.data.frame(dat.n2[,order(dat.n$Count)])
dat.n=dat.n[order(dat.n$Count),]
rownames(dat.n)=dat.n$Sample
my_sample_col=as.data.frame(dat.n$Sun.exposure)
my_sample_col
rownames(my_sample_col)=dat.n$Sample
colnames(my_sample_col)="Site"
my_sample_col$Site=factor(my_sample_col$Site,levels=c("Exposed", "Non-exposed"))
Var2 = c("#FF0000","#0000FF")
names(Var2) = c("Exposed", "Non-exposed")
ann_colors = list(Site=Var2)
#Figure 1b heatmap
pheatmap(dat.n2,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,annotation_col =my_sample_col,annotation_colors=ann_colors,color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(100),annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_legend = T,legend	=T)
pheatmap(dat.n2,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,annotation_col =my_sample_col,annotation_colors=ann_colors,color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(100),annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_legend = F,legend	=F,cellheight = 10)
?pheatmap

##age and sun-exposure state 
setwd("/data/Delta_data4/kysbbubbu/normal_skin/08_Counts")
dat.n=read.delim("Skin_annotation.txt",sep="\t",header=T)
head(dat.n)
dat.n2=as.data.frame(t(dat.n$UV_proportion))
rownames(dat.n2)=c("Age")
colnames(dat.n2)=dat.n$Sample
head(dat.n2)
dat.n2=as.data.frame(dat.n2[,order(dat.n$Count)])
dat.n=dat.n[order(dat.n$Count),]
rownames(dat.n)=dat.n$Sample
my_sample_col=as.data.frame(dat.n$Sun.exposure)
my_sample_col
rownames(my_sample_col)=dat.n$Sample
colnames(my_sample_col)="Site"
my_sample_col$Site=factor(my_sample_col$Site,levels=c("Exposed", "Non-exposed"))
Var2 = c("#FF0000","#0000FF")
names(Var2) = c("Exposed", "Non-exposed")
ann_colors = list(Site=Var2)
max(dat.n2)
#Figure 1b heatmap
pheatmap(dat.n2,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,annotation_col =my_sample_col,annotation_colors=ann_colors,color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(100),annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_legend = T,legend	=T)
pheatmap(dat.n2,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,annotation_col =my_sample_col,annotation_colors=ann_colors,color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(100),annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_legend = F,legend	=F,cellheight = 10)
?pheatmap

##Boxplot for mutation count 
setwd("/data/Delta_data4/kysbbubbu/normal_skin/08_Counts")
dat.n=read.delim("Skin_annotation.txt",sep="\t",header=T)
head(dat.n)
my_comparisons <- list( c("Exposed", "Non-exposed"))
dat.n$Sun.exposure=factor(dat.n$Sun.exposure,levels=c("Non-exposed","Exposed"))
#Figure 1c
ggboxplot(dat.n, x = "Sun.exposure", y = "Count",color = "Sun.exposure", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "Sun.exposure",outlier.shape = NA)+ theme_classic() +theme(axis.title.x = element_blank(),axis.text.x = element_text(size = 18,face="bold",color="black"),axis.text.y = element_text(size = 16), axis.title.y = element_text( size=19,face="bold"),legend.position = "none")+ scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+ylab("Mutation counts")+scale_y_log10()
wilcox.test(Count ~ Sun.exposure, data = dat.n)


#scatter plot for mutation count and age
setwd("/data/Delta_data4/kysbbubbu/normal_skin/08_Counts")
dat.n=read.delim("Skin_annotation.txt",sep="\t",header=T)
head(dat.n)
dat.n$Sun.exposure=factor(dat.n$Sun.exposure,levels=c("Non-exposed","Exposed"))
##Figure 1d
ggplot(dat.n, aes(x=Age, y=Count, color=Sun.exposure, shape=Sun.exposure)) + geom_point(size=2)+ labs(x = "Age (years)", y="Mutation counts")+ theme_classic()+theme(axis.text = element_text(size = 13), axis.title = element_text( size=16, face="bold"))+  scale_colour_manual(name="Skin_site", values= c("blue","red"))+scale_y_log10()+theme(legend.title=element_blank(),legend.text=element_text(size=14))+ theme(legend.position = "none")
cor.test(dat.n$Age,dat.n$Count)
cor.test(dat.n$Age,log10(dat.n$Count))


median(dat.n[dat.n$V11==c("Exposed"),]$V6)
max(dat.n[dat.n$V11==c("Exposed"),]$V6)
min(dat.n[dat.n$V11==c("Exposed"),]$V6)

median(dat.n[dat.n$V11==c("Non-exposed"),]$V6)
max(dat.n[dat.n$V11==c("Non-exposed"),]$V6)
min(dat.n[dat.n$V11==c("Non-exposed"),]$V6)

head(dat.n)


##Calculate VAFs 
setwd("/data/Delta_data4/kysbbubbu/normal_skin/08_Counts")
dat.n=read.delim("Skin.txt",sep="\t",header=T)
dat.n=dat.n[dat.n$Type1=="SNP",]
dat.n2=cbind(as.data.frame(aggregate(T_VAF ~ Sample, dat.n, FUN = max)),as.data.frame(aggregate(T_VAF ~ Sample, dat.n, FUN = median)),as.data.frame(aggregate(T_VAF ~ Sample, dat.n, FUN = mean)))
colnames(dat.n2)=c("max","median","mean")
fwrite(x=dat.n2,file="VAF.csv",row.names=T,col.names=T)


##VAF comparison
setwd("/data/Delta_data4/kysbbubbu/normal_skin/08_Counts")
dat.n=read.delim("Skin_annotation.txt",sep="\t",header=T)
dat.n$Sun.exposure=factor(dat.n$Sun.exposure,levels=c("Non-exposed","Exposed"))
dat.n2=dat.n[!dat.n$Sample %in% c("NS13","NS12","NS04","NS08","NS16","NS05","NS39","NS19"),]
###Figure 1e
#Comparison of median VAFs of SNVs between exposed and non-exposed group 
#400/600
ggboxplot(dat.n2, x = "Sun.exposure", y = "SNV_Median_VAF",color = "Sun.exposure", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "Sun.exposure",outlier.shape = NA)+ theme_classic() +theme(axis.title.x = element_blank(),axis.text.x = element_text(size = 18,face="bold",color="black"),axis.text.y = element_text(size = 16), axis.title.y = element_text( size=19, face="bold"))+ scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+ylab("Meidan VAF of SNVs (%)")+theme(legend.position="none")
wilcox.test(SNV_Median_VAF ~ Sun.exposure, data = dat.n2)
#600/500
#Correlation between Age and median VAF of SNVs 
setwd("/data/Delta_data4/kysbbubbu/normal_skin/08_Counts")
dat.n=read.delim("Skin_annotation.txt",sep="\t",header=T)
dat.n$Sun.exposure=factor(dat.n$Sun.exposure,levels=c("Non-exposed","Exposed"))
dat.n2=dat.n[!dat.n$Sample %in% c("NS13","NS12","NS04","NS08","NS16","NS05","NS39","NS19"),]
###Figure 1f
ggplot(dat.n2, aes(x=Age, y=SNV_Median_VAF, color=Sun.exposure, shape=Sun.exposure)) + geom_point(size=2)+ labs(x = "Age (years)", y="Meidan VAF of SNVs (%)")+ theme_classic()+theme(axis.text = element_text(size = 13), axis.title = element_text( size=16, face="bold"))+  scale_colour_manual(name="Sun.exposure", values= c("blue","red"))+theme(legend.title=element_blank(),legend.text=element_text(size=14))+ theme(legend.position = "none")
cor.test(dat.n2$Age,dat.n2$SNV_Median_VAF)



####VAF distribution 
setwd("/data/Delta_data4/kysbbubbu/normal_skin/08_Counts")
dat.n=read.delim("Skin.txt",sep="\t",header=T)
head(dat.n)
dat.n$Sun.exposure=factor(dat.n$Sun.exposure,levels=c("Non-exposed","Exposed"))
dia_me <- ddply(dat.n, .(Sun.exposure), numcolwise(median))
head(dia_mea)
#Figure S2
ggplot(dat.n, aes(x=T_VAF,color=Sun.exposure)) +geom_density(alpha=0.2)+geom_vline(data=dia_me, aes(xintercept=T_VAF, colour=Sun.exposure),linetype="dashed", size=0.2) +theme_classic()+theme(text = element_text(size=12))+scale_color_manual(values=c("#0000FF","#FF0000"))+theme(axis.title=element_text(color="black", size=16, face="bold"),axis.text.x=element_text(color="black", size=14, face="bold"),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=14, face="bold"), strip.text = element_text(size = 14, face="bold"))+ylab("Denstiy")+ scale_x_continuous(name = "VAF (%)")+ theme(legend.title=element_blank(),legend.text= element_text(size = 14))+theme(legend.position = "none")
#Figure S2b
ggboxplot(dat.n, x = "Sun.exposure", y = "T_VAF",color = "Sun.exposure", palette =c("#0000FF","#FF0000"),shape = "Sun.exposure",outlier.shape = NA)+ theme_classic() +theme(axis.title.x = element_blank(),axis.text.x = element_text(size = 18,face="bold",color="black"),axis.text.y = element_text(size = 16), axis.title.y = element_text( size=19, face="bold"))+ scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+ylab("VAF (%)")+theme(legend.position="none")+ylim(0,10)
wilcox.test(T_VAF ~ Sun.exposure, data = dat.n)

##UV-proportion
setwd("/data/Delta_data4/kysbbubbu/normal_skin/08_Counts")
dat.n=read.delim("Skin_annotation.txt",sep="\t",header=T)
dat.n$Sun.exposure=factor(dat.n$Sun.exposure,levels=c("Non-exposed","Exposed"))
dat.n2=dat.n[!dat.n$Sample %in% c("NS13","NS12","NS04","NS08","NS16","NS05","NS39","NS19"),]
###Figure S3
#400/600
#Comparison of UV proportion  between exposed and non-exposed group 
ggboxplot(dat.n2, x = "Sun.exposure", y = "UV_proportion",color = "Sun.exposure", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "Sun.exposure",outlier.shape = NA)+ theme_classic() +theme(axis.title.x = element_blank(),axis.text.x = element_text(size = 18,face="bold",color="black"),axis.text.y = element_text(size = 16), axis.title.y = element_text( size=19, face="bold"))+ scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+ylab("UV-signature (%)")+theme(legend.position="none")
wilcox.test(UV_proportion ~ Sun.exposure, data = dat.n2)
#600/500
###Figure S3
setwd("/data/Delta_data4/kysbbubbu/normal_skin/08_Counts")
dat.n=read.delim("Skin_annotation.txt",sep="\t",header=T)
dat.n$Sun.exposure=factor(dat.n$Sun.exposure,levels=c("Non-exposed","Exposed"))
dat.n2=dat.n[!dat.n$Sample %in% c("NS13","NS12","NS04","NS08","NS16","NS05","NS39","NS19"),]
#Correlation between Age and UV proportion
ggplot(dat.n2, aes(x=Age, y=UV_proportion, color=Sun.exposure, shape=Sun.exposure)) + geom_point(size=2)+ labs(x = "Age (years)", y="UV-signature (%)")+ theme_classic()+theme(axis.text = element_text(size = 13), axis.title = element_text( size=16, face="bold"))+  scale_colour_manual(name="Sun.exposure", values= c("blue","red"))+theme(legend.title=element_blank(),legend.text=element_text(size=14))+ theme(legend.position = "none")
cor.test(dat.n2$Age,dat.n2$UV_proportion)

###Figure 2. Significantly mutated genes
###Figure 2. Significantly mutated genes 
###Figure 2. Significantly mutated genes 
##Run dndscv
setwd("/data/Delta_data4/kysbbubbu/normal_skin/10_dndscv")
mutations=read.delim("dndscv_50.txt",header=TRUE,row.names=NULL)
target_genes=read.table("tcga_list.txt",header=FALSE)
target_genes=target_genes$V1
dndsskin = dndscv(mutations, gene_list=target_genes,outmats=T)
write.csv(dndsskin$sel_cv,"sel_cv_sensus.csv")


##Prepare dndscv input (duplicate removal)
setwd("/data/Delta_data4/kysbbubbu/normal_skin/10_dndscv")
mutations=read.delim("dndscv_50.txt",header=TRUE,row.names=NULL)
mutations = mutations[order(mutations$sampleID,mutations$chr,mutations$pos),]
ind = which(diff(mutations$pos)==1)
fwrite(x=as.data.frame(mutations[unique(sort(c(ind,ind+1))),]),file="duplicates.csv",row.names=T,col.names=T)
?read.delim


##Run dndscv
setwd("/data/Delta_data4/kysbbubbu/normal_skin/10_dndscv")
mutations=read.delim("dndscv_50_modify.txt",header=TRUE,row.names=NULL)
#mutations=mutations[duplicated(mutations[-1])==F,]
target_genes=read.table("TCGA_drivers.txt",header=FALSE)
target_genes=target_genes$V4
dndsskin = dndscv(mutations, gene_list=target_genes,outmats=T)
write.csv(dndsskin$sel_cv,"sel_cv_sensus.csv")

###SMG dndscv result visualization 
df=read.csv("sel_cv_sensus.csv",header=TRUE)
df=df[df$gene_name %in% c("NOTCH1","FAT1","TP53","PPM1D","KMT2D","ASXL1"),]
df2=df[,c(2,8,10,11)]
df2<- rename(df2,wmis_cv="Missense", wspl_cv="Nonsense\nSplicing", wind_cv="Indel") 
df3 <- reshape2::melt(df2, id.vars = "gene_name")
df3$gene_name=factor(df3$gene_name,levels=c("NOTCH1","FAT1","TP53","PPM1D","KMT2D","ASXL1"))
#1000/400
#Figure2a
ggplot(data=df3, aes(x=gene_name, y=value, fill=variable)) + geom_bar(stat="identity", color="black", position=position_dodge(),width=0.7)+theme_classic()+scale_fill_brewer(palette="Set1")+ylab("dN/dS")+theme(axis.title.x=element_blank(),axis.title.y = element_text( size=16, face="bold"),axis.text = element_text( size=13,,face="bold"))+xlab("Significantly mutated genes")+theme(legend.text=element_text(size=13),legend.title=element_blank())+theme(legend.position = "none")

##SMG count per genes 
setwd("/data/Delta_data4/kysbbubbu/normal_skin/10_dndscv")
dat.n=read.delim("SMG_list.txt",header=T,row.names=NULL)
dat.n$Type=plyr::mapvalues(x=dat.n$Type, from =c("Silent","Missense","Nonsense","Splicing","Frameshift","Inframe","Noncoding"), to =c("Silent","Missense","Nonsense","Splicing","Indels","Indels","Noncoding"))
dat.n=dat.n[dat.n$Type %in% c("Missense","Nonsense","Splicing","Indels"),]
table(dat.n$Gene)
dat.n$Type=factor(dat.n$Type,levels=rev(c("Missense","Nonsense","Splicing","Indels")))
dat.n2=as.data.frame(table(dat.n$Type,dat.n$Gene))
head(dat.n2)
dat.n2$Var2=factor(dat.n2$Var2,levels=c("NOTCH1","FAT1","TP53","PPM1D","KMT2D","ASXL1"))

#Figure2b
ggplot(data=dat.n2, aes(x=Var2, y=Freq, fill=Var1)) +  geom_bar(stat="identity", width=0.7)+  scale_fill_brewer(palette="Set1")+theme_classic() +theme(legend.text=element_text(size=13),legend.title=element_blank())+ylab("Mutation count")+theme(axis.title.x=element_blank(),axis.title.y = element_text( size=16, face="bold"),axis.text = element_text(size=13,face="bold"))+xlab("Significantly mutated genes")+theme(legend.position = "none")


##Preparation: Proportion of samples (%) according to the nonsilent mutation count 
setwd("/data/Delta_data4/kysbbubbu/normal_skin/10_dndscv")
dat.n=read.delim("SMG_skin.txt",header=T,row.names=NULL)
head(dat.n)
table(dat.n$V8)
dat.n=dat.n[dat.n$Type!="Noncoding",]
dat.n=dat.n[dat.n$Type!="Silent",]
head(dat.n)
df=as.data.frame(table(dat.n$Sample,dat.n$Gene))
head(df)
table(df$Freq)
df$Freq=factor(df$Freq,levels=c(0,1,2,3,4,5,6))
df=as.data.frame(df[df$Freq!=0,])
df$Freq=plyr::mapvalues(x=df$Freq, from =c(1,2,3,4,5,6), to =c(1,2,3,4,">4",">4"))
write.csv(df,"SMG_count.csv")

##Proportion of samples (%) according to the nonsilent mutation count 
setwd("/data/Delta_data4/kysbbubbu/normal_skin/10_dndscv")
dat.n=read.delim("SMG_count.txt",header=T,row.names=NULL)
dat.n$Mutation=factor(dat.n$Mutation,levels=c("NOTCH1","FAT1","TP53","PPM1D","KMT2D","ASXL1"))
dat.n$Number=factor(dat.n$Number,levels=c(">3","3","2","1","0"))
#Figure 2c
ggplot(data=dat.n, aes(x=Mutation, y=Count, fill=Number)) +  geom_bar(stat="identity", width=0.7,position = "fill")+theme_classic() +theme(legend.text=element_text(size=13),legend.title=element_blank())+theme(axis.title.x=element_blank(),axis.title.y = element_text(size=16, face="bold"),axis.text = element_text(size=13,face="bold"))+ylab("Proportion of samples (%)")+  scale_fill_brewer(palette="Paired",direction=-1)+ scale_y_continuous(labels = function(x) paste0(x*100))+theme(legend.position="none")


##VAF of SMGs
setwd("/data/Delta_data4/kysbbubbu/normal_skin/10_dndscv")
dat.n=read.delim("SMG_skin.txt",header=T,row.names=NULL)
head(dat.n)
table(dat.n$Sample)
dat.n=dat.n[dat.n$Type!="Noncoding",]
dat.n=dat.n[dat.n$Type!="Silent",]
dat.n$Gene=factor(dat.n$Gene,levels=c("NOTCH1","FAT1","TP53","PPM1D","KMT2D","ASXL1"))
head(dat.n)
dat.n2=dat.n[,c("Gene","T_VAF")]
#Figure 2d
ggboxplot(dat.n2, x = "Gene", y = "T_VAF", add = "jitter", outlier.shape = NA)+  theme_classic()+theme(legend.text=element_text(size=13),legend.title=element_blank())+theme(axis.title.x=element_blank(), axis.title.y= element_text( size=16, face="bold"),axis.text = element_text(size=13,face="bold"))+xlab("Significantly mutated genes")+ylab("VAF (%)")+theme(legend.position = "none")


median(dat.n2$T_VAF)
min(dat.n2$T_VAF)
max(dat.n2$T_VAF)

median(dat.n2[dat.n2$Gene=="NOTCH1",]$T_VAF)
min(dat.n2[dat.n2$Gene=="NOTCH1",]$T_VAF)
max(dat.n2[dat.n2$Gene=="NOTCH1",]$T_VAF)

median(dat.n2[dat.n2$Gene=="FAT1",]$T_VAF)
min(dat.n2[dat.n2$Gene=="FAT1",]$T_VAF)
max(dat.n2[dat.n2$Gene=="FAT1",]$T_VAF)

median(dat.n2[dat.n2$Gene=="TP53",]$T_VAF)
min(dat.n2[dat.n2$Gene=="TP53",]$T_VAF)
max(dat.n2[dat.n2$Gene=="TP53",]$T_VAF)

median(dat.n2[dat.n2$Gene=="PPM1D",]$T_VAF)
min(dat.n2[dat.n2$Gene=="PPM1D",]$T_VAF)
max(dat.n2[dat.n2$Gene=="PPM1D",]$T_VAF)

median(dat.n2[dat.n2$Gene=="KMT2D",]$T_VAF)
min(dat.n2[dat.n2$Gene=="KMT2D",]$T_VAF)
max(dat.n2[dat.n2$Gene=="KMT2D",]$T_VAF)

median(dat.n2[dat.n2$Gene=="ASXL1",]$T_VAF)
min(dat.n2[dat.n2$Gene=="ASXL1",]$T_VAF)
max(dat.n2[dat.n2$Gene=="ASXL1",]$T_VAF)

#Correlation between ages (years) and non-silent mutations in drvier genes
setwd("/data/Delta_data4/kysbbubbu/normal_skin/08_Counts")
dat.n=read.delim("Skin_annotation.txt",header=T,row.names=NULL)
dat.n2=dat.n[!dat.n$Sample %in% c("NS13","NS12","NS04","NS08","NS16","NS05","NS39","NS19"),]
colnames(dat.n2)
dat.n$Sun.exposure=factor(dat.n$Sun.exposure,levels=c("Non-exposed","Exposed"))
#Figure S5a
ggplot(dat.n, aes(x=Age, y= Driver, color=Sun.exposure, shape=Sun.exposure)) + geom_point(size=2)+ labs(x = "Age (years)", y="Non-silent mutations\nin driver genes (n)")+ theme_classic()+theme(text = element_text(size=14))+theme(axis.text.x = element_text( size = 14,face="bold"),axis.text.y = element_text(size = 14,face="bold"), axis.title.x = element_text( size=14, face="bold"),axis.title.y = element_text( size=14, face="bold"))+  scale_colour_manual(name="Sun.exposure", values= c("blue", "red"))+theme(legend.title=element_blank(),legend.text=element_text(size=14))+ theme(legend.position = "none")
cor.test(dat.n$Age,dat.n$Driver)

#Comparison of VAFs between samples with and without driver mutations 
setwd("/data/Delta_data4/kysbbubbu/normal_skin/08_Counts")
dat.n=read.delim("Skin.txt",header=T,row.names=NULL)
dat.n=dat.n[!dat.n$Sample %in% c("NS13","NS12","NS04","NS08","NS16","NS05","NS39","NS19"),]
table(dat.n$Type1)
dat.n2=dat.n[dat.n$Type1=="SNP",]
table(dat.n2$Type1)
median(dat.n2[dat.n2$Driver=="With_driver",]$T_VAF)
median(dat.n2[dat.n2$Driver=="Without_driver",]$T_VAF)
dat.n2$Driver=plyr::mapvalues(x=dat.n2$Driver, from =c("With_driver","Without_driver"), to =c("Samples with\nnon-silent mutations\nin driver genes","Samples without\nnon-silent mutations\nin driver genes"))
#Figure S5b
ggboxplot(dat.n2, x = "Driver", y = "T_VAF", color = "Driver", palette =c("#FF0000","#0000FF"), shape = "Driver", outlier.shape = NA, bxp.errorbar = T)+  theme_classic()+theme(legend.text=element_text(size=13),legend.title=element_blank())+theme(axis.title.x=element_blank(), axis.title.y= element_text( size=16, face="bold"),axis.text = element_text(size=13,face="bold"))+xlab("Genes")+ylab("VAFs of SNVs (%)")+ylim(0,8)+theme(legend.position = "none")
wilcox.test(T_VAF ~ Driver, data = dat.n2)

###Figure 3. Protein loci of non-silent mutations in driver gene 
###Figure 3. Protein loci of non-silent mutations in driver gene 
###Figure 3. Protein loci of non-silent mutations in driver gene 

#Preparation for figure 3 (normal skin)
setwd("/data/Delta_data4/kysbbubbu/normal_skin/08_Counts")
dat.n=read.delim("Skin.txt",sep="\t",header=T)
dat.n2=dat.n[dat.n$Gene %in% c("NOTCH1","FAT1","TP53","PPM1D","KMT2D","ASXL1"),]
write.csv(dat.n2,"SMG_list.csv")

#Preparation for figure 3 (cSCC)
setwd("/data/Delta_data4/kysbbubbu/normal_skin/11_SMG")
dat.n=read.delim("41525_2021_226_MOESM5_ESM.txt",sep="\t",header=T)
head(dat.n)
dat.n2=dat.n[dat.n$Hugo_Symbol %in% c("NOTCH1","FAT1","TP53","PPM1D","KMT2D","ASXL1"),]
write.csv(dat.n2,"SMG_cSCC.csv",row.names=F,col.names=T)


###Figure 4. Mutational signature 
###Figure 4. Mutational signature 
###Figure 4. Mutational signature 

##prepare mutational patterns input 
setwd("/data/Delta_data4/kysbbubbu/normal_skin/09_Patterns")
dat.n=read.delim("dndscv_50.txt",sep="\t",header=F)
head(dat.n)
sample=unique(dat.n$V8)
for (list1 in sample) {
  mutation1=dat.n[dat.n$V8==list1,]
  colnames(mutation1)=c("#CHROM","POS",	"ID",	"REF",	"ALT",	"QUAL",	"FILTER",	"INFO",	"FORMAT",	"NORMAL",	"TUMOR")
  write.table(x=mutation1,  quote = FALSE, col.names	=TRUE, row.names=FALSE, file=list1, sep="\t")}

#96-context input preparation for Sigprofiler
setwd("/data/Delta_data4/kysbbubbu/normal_skin/09_Patterns")
dat.n=read.delim("dndscv_50.txt",sep="\t",header=F)
table(dat.n$V1)
head(dat.n)
sigs.input <- mut.to.sigs.input(mut.ref =dat.n, sample.id = "V8", chr = "V1", pos = "V2", ref = "V4", alt = "V5", bsg = BSgenome.Hsapiens.UCSC.hg19)
dim(dat.n)

##prepare Sigprofiler input
fwrite(x=as.data.frame(t(sigs.input)),file="Sigprofiler.csv",row.names=T,col.names=T)


#96-context input preparation for Sigprofiler
setwd("/data/Delta_data4/kysbbubbu/normal_skin/09_Patterns")
dat.n=read.delim("dndscv_50.txt",sep="\t",header=F)
sigs.input <- mut.to.sigs.input(mut.ref =dat.n, sample.id = "V8", chr = "V1", pos = "V2", ref = "V4", alt = "V5", bsg = BSgenome.Hsapiens.UCSC.hg19)
fwrite(x=as.data.frame(t(sigs.input)),file="Sigprofiler.csv",row.names=T,col.names=T)

#96-context input preparation for Sigprofiler
setwd("/data/Delta_data4/kysbbubbu/normal_skin/09_Patterns")
dat.n=read.delim("tumor_mutation_list.txt",sep="\t",header=F)
sigs.input <- mut.to.sigs.input(mut.ref =dat.n, sample.id = "V8", chr = "V1", pos = "V2", ref = "V4", alt = "V5", bsg = BSgenome.Hsapiens.UCSC.hg19)
fwrite(x=as.data.frame(t(sigs.input)),file="Sigprofiler_tumor.csv",row.names=T,col.names=T)

##Run sigprofiler
sigprofilerextractor("matrix", "/data/Delta_data4/kysbbubbu/normal_skin/07_Sigprofiler/Epidermis", "/data/Delta_data4/kysbbubbu/normal_skin/07_Sigprofiler/Epidermis/Sigprofiler.txt", exome=T,context_type="96",cosmic_version=3.2)
sigprofilerextractor("matrix", "/data/Delta_data4/kysbbubbu/normal_skin/07_Sigprofiler/Tumor", "/data/Delta_data4/kysbbubbu/normal_skin/07_Sigprofiler/Tumor/Sigprofiler_tumor.txt", exome=T,context_type="96",cosmic_version=3.2)

##visualization 96-mutational context of each sample 
setwd("/data/Delta_data4/kysbbubbu/normal_skin/09_Patterns")
sample=read.delim("vcf_list.txt",sep="\t",header=F)
vcf_files =sample$V4
sample_names=sample$V3
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
tissues=rep(c("skin"),24)
grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome,predefined_dbs_mbs=T)
mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)
mut_mat2=mut_mat[,order(colSums(mut_mat))]
colSums(mut_mat2)
#1000/1600
plot_96_profile(mut_mat2[,1:12], ymax = 0.2, condensed =F)+ theme(axis.text.x= element_text(size = 8,color="black"), axis.text.y= element_text(size = 10,color="black"), axis.title.x = element_text(size=17, face="bold"),axis.title.y = element_text(size=17, face="bold"))+theme(legend.title =element_blank(),legend.text = element_text(size=16), legend.position="right") + theme(strip.text.x = element_text(size = 14,face="bold"),strip.text.y = element_text(size = 14,face="bold"))+xlab("Mutational context") 
plot_96_profile(mut_mat2[,13:24], ymax = 0.2, condensed =F)+ theme(axis.text.x= element_text(size = 8,color="black"), axis.text.y= element_text(size = 10,color="black"), axis.title.x = element_text(size=17, face="bold"),axis.title.y = element_text(size=17, face="bold"))+theme(legend.title =element_blank(),legend.text = element_text(size=16), legend.position="right") + theme(strip.text.x = element_text(size = 14,face="bold"),strip.text.y = element_text(size = 14,face="bold"))+xlab("Mutational context") 


#De novo extracted signature visualization 
setwd("/data/Delta_data4/kysbbubbu/normal_skin/09_Patterns")
dat.n=read.delim("SBS96_De-Novo_Signatures.txt",sep="\t",row.names="Type",header=T)
head(dat.n)
mut_mat
##mutational context of de novo extracted signatures
dat.n2=as.data.frame(dat.n[,"SBS96A"])
rownames(dat.n2)=rownames(dat.n)
#Figure 4a
#1000/300
plot_96_profile(dat.n2, ymax = 0.2, condensed =F)+ theme(axis.text.x= element_text(size = 8,color="black"), axis.text.y= element_text(size = 12,color="black"), axis.title.x = element_text(size=16, face="bold"),axis.title.y = element_text(size=16, face="bold"))+theme(legend.title =element_blank(),legend.text = element_text(size=16), legend.position="right") + theme(strip.text.x = element_text(size = 14,face="bold"),strip.text.y = element_blank())+xlab("Mutational context")  + scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2"))

##Strand bias visualization
setwd("/data/Delta_data4/kysbbubbu/normal_skin/09_Patterns")
sample=read.delim("vcf_list.txt",sep="\t",header=F)
vcf_files =sample$V4
sample_names=sample$V3
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
tissues=rep(c("skin"),24)
grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome,predefined_dbs_mbs=T)
genes_hg19 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
mut_mat_s <- mut_matrix_stranded(grl, ref_genome, genes_hg19)
strand_counts <- strand_occurrences(mut_mat_s, by = tissues)
strand_bias <- strand_bias_test(strand_counts)
?strand_bias_test
strand <- mut_strand(grl[[1]], genes_hg19)
#Figure 4b
plot_strand(strand_counts, mode = "relative")+ theme(axis.text.x= element_text(size = 12,color="black"), axis.text.y= element_text(size = 12,color="black"), axis.title.x = element_text(size=16, face="bold"),axis.title.y = element_text(size=16, face="bold"))+theme(legend.title =element_blank(),legend.text = element_text(size=16), legend.position="right") + theme(strip.text = element_blank())+xlab("Mutational type")+theme(legend.position="none")
strand_bias <- strand_bias_test(strand_counts)
?p.adjust


##mutational context of the average mutational signature
dat.n3=as.data.frame(dat.n[,"Average"])
rownames(dat.n3)=rownames(dat.n)
#Figure S7a
#1000/300
plot_96_profile(dat.n3, ymax = 0.2, condensed =F)+ theme(axis.text.x= element_text(size = 8,color="black"), axis.text.y= element_text(size = 12,color="black"), axis.title.x = element_text(size=16, face="bold"),axis.title.y = element_text(size=16, face="bold"))+theme(legend.title =element_blank(),legend.text = element_text(size=16), legend.position="right") + theme(strip.text.x = element_text(size = 14,face="bold"),strip.text.y = element_blank())+xlab("Mutational context")  + scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2"))
+
cos.sim=function(ma, mb){
  mat=tcrossprod(ma, mb)
  t1=sqrt(apply(ma, 1, crossprod))
  t2=sqrt(apply(mb, 1, crossprod))
  mat / outer(t1,t2)}
cos.sim(t(dat.n2),t(dat.n3))


#Non-exposed skin visualization
setwd("/data/Delta_data4/kysbbubbu/normal_skin/09_Patterns")
sample=read.delim("vcf_list.txt",sep="\t",header=F)
vcf_files =sample$V4
sample_names=sample$V1
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
tissues=rep(c("skin"),24)
grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome,predefined_dbs_mbs=T)
mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)
sigs2=mut_mat[,c("NS41","NS49")]
#FigureS9
plot_96_profile(sigs2, ymax = 0.2, condensed =F)+ theme(axis.text.x= element_text(size = 8,color="black"), axis.text.y= element_text(size = 13,color="black"), axis.title.x = element_text(size=16, face="bold"),axis.title.y = element_text(size=16, face="bold"))+theme(legend.title =element_blank(),legend.text = element_text(size=16), legend.position="right") + theme(strip.text.x = element_text(size = 14,face="bold"),strip.text.y = element_text(size = 14,face="bold"))+xlab("Mutational context") 

#De novo extracted signature visualization 
setwd("/data/Delta_data4/kysbbubbu/normal_skin/09_Patterns")
dat.n=read.delim("SBS96_De-Novo_Signatures.txt",sep="\t",row.names="Type",header=T)
head(dat.n)
mut_mat
##mutational context of de novo extracted signatures
dat.n2=as.data.frame(dat.n[,"SBS96A"])
rownames(dat.n2)=rownames(dat.n)
dim(dat.n2)
dim(sigs2)
cos.sim(as.matrix(t(dat.n2)),as.matrix(t(sigs2)))


###piechart visualization
barplot = as.data.frame(c("SBS1", "SBS5", "SBS7a", "SBS7b"))
colnames(barplot)="SBS"
# Signature SBS1 (3.98%) & Signature SBS5 (17.84%) & Signature SBS7a (21.62%) & Signature SBS7b (56.56%)
barplot$proportion=c("0.039","0.179","0.216","0.566")
barplot$proportion=as.numeric(as.character(barplot$proportion))
library(scales)
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )
barplot$SBS=factor(barplot$SBS, levels= c("SBS1", "SBS5", "SBS7a", "SBS7b"))
#Figure 4b
ggplot(data = barplot, aes(x = "", y = proportion, fill = SBS )) +   geom_bar(stat = "identity", position = position_fill()) +  blank_theme + coord_polar(theta = "y")  +    theme(axis.text=element_blank(),legend.position='right',legend.title=element_blank(),legend.text=element_text(size=14))+ scale_fill_manual(values=c(brewer.pal(n = 4, name ="Set2")))


#Mutation count
##800/150
dat.n=t(mut_matrix(vcf_list = grl, ref_genome = ref_genome))
dat.n=as.matrix(dat.n[order(rowSums(dat.n)),])
Mut=as.data.frame(rowSums(dat.n))
colnames(Mut)="SNV"
min(Mut)
max(Mut)
#Figure 4c
#800/150
ggplot(Mut, aes(x=factor(1:length(SNV)), y=SNV)) +theme_classic()+geom_bar(stat = "identity", width=1,colour="black")+theme(text = element_text(size=12)) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=12, face="bold"))+ theme(legend.position = "none",panel.grid.major.y=element_line(linetype="dotted",colour="black"))+ scale_y_continuous(trans="log10",limits=c(100,3000),oob = rescale_none)

##Heatmap visualization of 96-mutaiotnal context 
dat.n <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)
my_sample_col=as.data.frame(colnames(dat.n))
colnames(my_sample_col)="Site"
rownames(my_sample_col)=my_sample_col$Site
my_sample_col$Site=sample$V5
my_gene_col <- data.frame(sample = rep(c("C>A", "C>G","C>T","T>A","T>C","T>G"), c(16,16,16,16,16,16)))
rownames(my_gene_col)=rownames(dat.n)
colnames(my_gene_col)="Type"
dat.n = dat.n[,order(colSums(dat.n))]
Var1 = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")
names(Var1) = c("C>A", "C>G","C>T","T>A","T>C","T>G")
Var2 = c("#FF0000","#0000FF")
names(Var2) = c("Exposed", "Non-exposed")
ann_colors = list(Type = Var1, Site=Var2)

##900/600
dat.n2=as.data.frame(prop.table(dat.n,margin=2))
colSums(dat.n2)
#Figure 4c
pheatmap(dat.n2,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 0.2, by = 0.002),annotation_col =my_sample_col,annotation_row =my_gene_col,colorRampPalette(brewer.pal(n = 9, name ="YlOrRd"))(100),annotation_names_row=FALSE,annotation_names_col=FALSE,annotation_colors=ann_colors, annotation_legend = FALSE,legend	=FALSE)

pheatmap(dat.n2,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 0.2, by = 0.002),annotation_col =my_sample_col,annotation_row =my_gene_col,colorRampPalette(brewer.pal(n = 9, name ="YlOrRd"))(100),annotation_names_row=FALSE,annotation_names_col=FALSE,annotation_colors=ann_colors, annotation_legend = TRUE,legend	=T)

pheatmap(dat.n,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 0.2, by = 0.002),annotation_col =my_sample_col,annotation_row =my_gene_col,colorRampPalette(brewer.pal(n = 9, name ="YlOrRd"))(100),annotation_names_row=FALSE,annotation_names_col=FALSE,annotation_colors=ann_colors, annotation_legend = TRUE,legend	=T,cellwidth = 0.1)


#Cosine similarity
##800/150
sample=read.delim("vcf_list.txt",sep="\t",header=F)
vcf_files =sample$V4
sample_names=sample$V1
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
tissues=rep(c("skin"),24)
grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome,predefined_dbs_mbs=T)
dat.n <- t(mut_matrix(vcf_list = grl, ref_genome = ref_genome))
dat.n=as.matrix(dat.n[order(rowSums(dat.n)),])
sig=as.data.frame(prop.table(dat.n,margin=1))
dat.n=read.delim("SBS96_De-Novo_Signatures.txt",sep="\t",row.names="Type",header=T)
dat.n2=as.data.frame(t(dat.n[,"SBS96A"]))
colnames(dat.n2)=rownames(dat.n)
dat.n2=as.matrix(dat.n2)
sig=as.matrix(sig)
cos.matrix=as.data.frame(t(cos.sim(dat.n2,sig)))
colnames(cos.matrix)="cos.sim"
#Figure 4d
#800/100
ggplot(cos.matrix, aes(x=factor(1:length(cos.sim)), y=cos.sim)) +theme_classic()+geom_bar(stat = "identity", width=0.7,colour="steelblue",fill="steelblue")+theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=15, face="bold")) +ylab("Cosine similarity")+ theme(legend.position = "none")+scale_y_continuous(limits=c(0,1),breaks = c(0,0.9))+ geom_hline(yintercept=c(0.9), linetype="dotted")

##Heatmap/Barplot visualization of proportions of mutational types 
type_occurrences <- mut_type_occurrences(grl, ref_genome)
dat.n=type_occurrences[,1:6]
sample=read.delim("vcf_list.txt",sep="\t",header=F)
dat.n=as.matrix(dat.n[order(rowSums(dat.n)),])
dat.n2=as.data.frame(t(prop.table(dat.n,margin=1)))
my_sample_col=as.data.frame(sample$V5)
colnames(my_sample_col)="Site"
rownames(my_sample_col)=sample$V3
my_gene_col <- data.frame(Type =c("C>A", "C>G","C>T","T>A","T>C","T>G"))
rownames(my_gene_col)=my_gene_col$Type
Var1 = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")
names(Var1) =c("C>A", "C>G","C>T","T>A","T>C","T>G")
Var2 = c("#FF0000","#0000FF")
names(Var2) = c("Exposed", "Non-exposed")
ann_colors = list(Type = Var1, Site=Var2)
pheatmap(dat.n2,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 1, by = 0.01),annotation_col =my_sample_col,annotation_row =my_gene_col,colorRampPalette(brewer.pal(n = 9, name ="YlOrRd"))(100),annotation_names_row=FALSE,annotation_names_col=FALSE,annotation_colors=ann_colors, annotation_legend = TRUE,legend	=T)
dat.n2$Sig=rownames(dat.n2)
dat.n2 <- reshape2::melt(dat.n2, id.vars = "Sig")
#Figure 4d
ggplot(dat.n2, aes(x=variable, y=value, fill=Sig)) +geom_bar(stat="identity",position = "fill", width=0.7) + theme_classic() + theme(axis.title.x = element_blank(),axis.text.x= element_text(size = 15,face="bold",color="black",angle=90,vjust=0.5),axis.text.y= element_text(size = 13,face="bold",color="black"), axis.title.y = element_text(size=16, face="bold"))+theme(legend.title =element_blank(),legend.text = element_text(size=16), legend.position="right")+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2"))+ scale_y_continuous(breaks = scales::pretty_breaks(n = 2))+ylab("Proportion of SNVs (%)")+scale_y_continuous(labels = function(x) paste0(x*100))+ theme(legend.position = "none")

##Heatmap visualization of proportions of mutational signatures
setwd("/data/Delta_data4/kysbbubbu/normal_skin/09_Patterns")
dat.n=as.data.frame(read.delim("COSMIC_SBS96_Activities.txt",sep="\t",header=T,row.names="Samples"))
head(dat.n)
dat.n=as.matrix(dat.n[order(rowSums(dat.n)),])
dat.n2=as.data.frame(t(prop.table(dat.n,margin=1)))
my_sample_col=as.data.frame(sample$V5)
sample=read.delim("vcf_list.txt",sep="\t",header=F)
my_sample_col=as.data.frame(sample$V5)
colnames(my_sample_col)="Site"
rownames(my_sample_col)=sample$V3
my_gene_col <- data.frame(Sig =c("SBS1", "SBS5","SBS7a","SBS7b"))
rownames(my_gene_col)=my_gene_col$Sig
Var1 = brewer.pal(n = 4, name ="Set2")
names(Var1) = c("SBS1", "SBS5","SBS7a","SBS7b")
Var2 = c("#FF0000","#0000FF")
names(Var2) = c("Exposed", "Non-exposed")
ann_colors = list(Sig = Var1, Site=Var2)
pheatmap(dat.n2,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(0, 1, by = 0.01),annotation_col =my_sample_col,annotation_row =my_gene_col,colorRampPalette(brewer.pal(n = 9, name ="YlOrRd"))(100),annotation_names_row=FALSE,annotation_names_col=FALSE,annotation_colors=ann_colors, annotation_legend = TRUE,legend	=T)
dat.n2$Sig=rownames(dat.n2)
dat.n2 <- reshape2::melt(dat.n2, id.vars = "Sig")
#Figure 4d
ggplot(dat.n2, aes(x=variable, y=value, fill=Sig)) +geom_bar(stat="identity",position = "fill", width=0.7) + theme_classic() + theme(axis.title.x = element_blank(),axis.text.x= element_text(size = 15,face="bold",color="black",angle=90,vjust=0.5),axis.text.y= element_text(size = 13,face="bold",color="black"), axis.title.y = element_text(size=16, face="bold"))+theme(legend.title =element_blank(),legend.text = element_text(size=16), legend.position="right")+scale_fill_manual(values=brewer.pal(4,"Set2"))+ scale_y_continuous(breaks = scales::pretty_breaks(n = 2))+ylab("Proportion of SNVs (%)")+scale_y_continuous(labels = function(x) paste0(x*100))+theme(legend.position = "none")

##DBS context
setwd("/data/Delta_data4/kysbbubbu/normal_skin/09_Patterns")
sample=read.delim("vcf_list.txt",sep="\t",header=F)
vcf_files =sample$V4
sample_names=sample$V1
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
tissues=rep(c("skin"),24)
grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome,predefined_dbs_mbs=F)
?read_vcfs_as_granges
head(grl[[1]])
dbs_grl <- get_mut_type(grl, type = "dbs",predefined_dbs_mbs = T)
head(dbs_grl[[1]])
dbs_grl <- get_dbs_context(dbs_grl)
head(dbs_grl[[1]])
dbs_counts <- count_dbs_contexts(dbs_grl)
plot_dbs_contexts(dbs_counts, same_y = TRUE)
plot_main_dbs_contexts(dbs_counts, same_y = TRUE)
dbs_counts2=as.data.frame(rowSums(dbs_counts))
colnames(dbs_counts2)="Proportion"
colSums(dbs_counts2)
dbs_counts2$Proportion=dbs_counts2$Proportion/1287
#1200/300
plot_dbs_contexts(dbs_counts2, same_y = TRUE, condensed = F)+ theme(axis.text.x= element_text(size = 8,color="black"), axis.text.y= element_text(size = 10,color="black"), axis.title.x = element_text(size=17, face="bold"),axis.title.y = element_text(size=17, face="bold"))+theme(legend.title =element_blank(),legend.text = element_text(size=16), legend.position="right") + theme(strip.text.y = element_blank(),strip.text.x = element_text(size = 14,face="bold"))+xlab("Mutational context")+ylab("Proportion of DBS")

#600/400
plot_main_dbs_contexts(dbs_counts2, same_y = TRUE)+ theme(axis.text.x= element_text(size = 14,color="black"), axis.text.y= element_text(size = 14,color="black"), axis.title.x = element_text(size=17, face="bold"),axis.title.y = element_text(size=17, face="bold"))+theme(legend.title =element_blank(),legend.text = element_text(size=16), legend.position="right") + theme(strip.text.y = element_blank(),strip.text.x = element_text(size = 14,face="bold"))+xlab("Mutational type")+ylab("Proportion of DBS")

###Figure 5 CNA visualization
###Figure 5 CNA visualization
###Figure 5 CNA visualization

setwd("/data/Delta_data4/kysbbubbu/normal_skin/06_FACET")
dat.n=read.delim("copy_number.txt",sep="\t",header=T)
dat.n=dat.n[dat.n$sample2=="Epidermis",]
dat.n=as.data.frame(dat.n[,1:8])
#1000/600, N=39
#Figure S10a
plotAberration(dat.n,thres.gain=0.25,sample.labels=TRUE,sep.samples=0)


setwd("/data/Delta_data4/kysbbubbu/normal_skin/06_FACET")
dat.n=read.delim("copy_number.txt",sep="\t",header=T)
dat.n=dat.n[dat.n$sample2=="Tumor",]
dat.n=as.data.frame(dat.n[,1:8])
head(dat.n)
#1000/220
#Figure S10a
plotAberration(dat.n,thres.gain=0.25,sample.labels=TRUE,sep.samples=0)

##Calculate legnth 
setwd("/data/Delta_data4/kysbbubbu/normal_skin/06_FACET")
dat.n=read.delim("copy_number.txt",sep="\t",header=T)
dat.n2=dat.n[dat.n$mean!=0,]
dat.n3=as.data.frame(aggregate(length ~ sampleID, dat.n2, FUN = sum))
fwrite(x=as.data.frame(dat.n3),file="CNA.csv",row.names=T,col.names=T)

setwd("/data/Delta_data4/kysbbubbu/normal_skin/08_Counts")
dat.n=read.delim("Skin_annotation.txt",sep="\t",header=T)
dat.n2=dat.n[dat.n$Tumor..analyzed=="O",c("Sample","CNA_skin ","CNA_tumor")]
dat.n2=dat.n[,c("Sample","CNA_skin","CNA_tumor")]
colnames(dat.n2)=c("Case","Skin","Cancer")
dat.n3 <- reshape2::melt(dat.n2, id.vars = "Case")
dat.n3$value=dat.n3$value/3036303846*100
#Figure S10b
ggboxplot(dat.n3, x = "variable", y = "value", color = "variable", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "variable", outlier.shape = NA)+  theme_classic()+theme(legend.text=element_text(size=13),legend.title=element_blank())+theme(axis.title.x=element_blank(), axis.title.y= element_text( size=18, face="bold"),axis.text.x= element_text( size=18, face="bold"),axis.text.y = element_text(size=15,face="bold"))+ylab("Genome with CNA (%)")+theme(legend.position = "none")
wilcox.test(value ~ variable, data = dat.n3)

###Figure 6 comparison of normal skin and paired cancer 
###Figure 6 comparison of normal skin and paired cancer 
###Figure 6 comparison of normal skin and paired cancer 

#Confirm mutational counts for each sample 
setwd("/data/Delta_data4/kysbbubbu/normal_skin/08_Counts")
dat.n=read.delim("Skin.txt",sep="\t",header=T)
fwrite(x=as.data.frame(table(dat.n$Sample)),file="sample_count.csv",row.names=T,col.names=T)
fwrite(x=as.data.frame(table(dat.n$Sample,dat.n$UV.signature)),file="UV_count.csv",row.names=T,col.names=T)

dat.n=read.delim("Tumor.txt",sep="\t",header=T)
head(dat.n)
fwrite(x=as.data.frame(table(dat.n$Sample)),file="tumor_count.csv",row.names=T,col.names=T)
fwrite(x=as.data.frame(table(dat.n$Sample,dat.n$UV.signature)),file="tumor_UV_count.csv",row.names=T,col.names=T)

dat.n=dat.n[dat.n$Type1=="SNP",]
dat.n2=cbind(as.data.frame(aggregate(T_VAF ~ Sample, dat.n, FUN = max)),as.data.frame(aggregate(T_VAF ~ Sample, dat.n, FUN = median)),as.data.frame(aggregate(T_VAF ~ Sample, dat.n, FUN = mean)))
fwrite(x=as.data.frame(dat.n2),file="tumor_VAF.csv",row.names=T,col.names=T)

##mutation count comparision
setwd("/data/Delta_data4/kysbbubbu/normal_skin/08_Counts")
dat.n=read.delim("Skin_annotation.txt",sep="\t",header=T)
dat.n2=dat.n[dat.n$Tumor..analyzed=="O",c("Sample","Count","Tumor_count")]
colnames(dat.n2)=c("Case","Skin", "Cancer")
dat.n3 <- reshape2::melt(dat.n2, id.vars = "Case")
dat.n3$Case=factor(dat.n3$Case,levels=c("NS39","NS32","NS41","NS28","NS47","NS42","NS43","NS18","NS45","NS09"))
dat.n3$value=as.numeric(dat.n3$value)
#Figure 6a
ggplot(data=dat.n3, aes(x=Case, y=value, fill=variable))  +  geom_bar(stat="identity", position=position_dodge(),width=0.7)+ scale_fill_brewer(palette="Set1",direction=-1)+theme_classic() +theme(axis.text.x = element_text(size = 14, face="bold",angle=0,vjust=0.5),axis.text.y = element_text(size = 12), axis.title.x = element_blank(),axis.title.y = element_text( size=16, face="bold"), legend.text=element_text(size=12),legend.title=element_blank())+ylab("Mutation counts")+theme(legend.position = "none")
#Figure S11a
ggboxplot(dat.n3, x = "variable", y = "value", color = "variable", palette =c("#0000FF","#FF0000"),  add = "jitter", shape = "variable", outlier.shape = NA)+  theme_classic()+theme(legend.text=element_text(size=13),legend.title=element_blank())+theme(axis.title.x=element_blank(), axis.title.y= element_text( size=18, face="bold"),axis.text.x= element_text( size=18, face="bold"),axis.text.y = element_text(size=15,face="bold"))+xlab("Samples")+ylab("Mutation counts")+theme(legend.position = "none")+scale_y_log10(limits=c(10,10000))
?scale_y_log10
wilcox.test(value ~ variable, data = dat.n3)

?wilcox.test

setwd("/data/Delta_data4/kysbbubbu/normal_skin/08_Counts")
dat.n=read.delim("Skin_annotation.txt",sep="\t",header=T)
head(dat.n)
dat.n2=dat.n[dat.n$Tumor..analyzed=="O",c("Sample","Count","Tumor_count")]
colnames(dat.n2)=c("Case","Skin","Cancer")
dat.n2$Case=factor(dat.n2$Case,levels=c("NS39","NS32","NS41","NS28","NS47","NS42","NS43","NS18","NS45","NS09"))
dat.n2$Disease2=c(rep("BCC",4),rep("cSCC in situ",2),rep("BCC",4))
dat.n2$Skin=as.numeric(dat.n2$Skin)
dat.n2$Cancer=as.numeric(dat.n2$Cancer)
#Figure S6b
ggplot(dat.n2, aes(x=Skin, y=Cancer, color=Disease2, shape=Disease2)) + geom_point(size=3)+ labs(x = "Mutation counts (skin)", y="Mutation counts (cancer)")+ theme_classic()+theme(axis.title = element_text( size=16, face="bold"),axis.text = element_text( size=13,face="bold"),legend.text=element_text(size=14))+  scale_colour_manual(name="Disease2", values= c("blue", "red"))+scale_x_log10()+scale_y_log10()+theme(legend.title=element_blank())+ theme(legend.position = "none")
cor.test(dat.n2$Epidermis,dat.n2$Tumor, method = "pearson")
cor.test(log10(dat.n2$Epidermis),log10(dat.n2$Tumor), method = "pearson")



##UV proportion comparision
setwd("/data/Delta_data4/kysbbubbu/normal_skin/08_Counts")
dat.n=read.delim("Skin_annotation.txt",sep="\t",header=T)
head(dat.n)
dat.n[,18]
dat.n2=as.data.frame(dat.n[dat.n$Tumor..analyzed=="O",c(1,18,39)])
colnames(dat.n2)=c("Case","Epidermis","Tumor")
dat.n3 <- reshape2::melt(dat.n2, id.vars = "Case")
dat.n3$Case=factor(dat.n3$Case,levels=c("NS39","NS32","NS41","NS28","NS47","NS42","NS43","NS18","NS45","NS09"))
dat.n3$value=as.numeric(dat.n3$value)
#Figure S11b
ggplot(data=dat.n3, aes(x=Case, y=value, fill=variable))  +  geom_bar(stat="identity", position=position_dodge(),width=0.7)+ scale_fill_brewer(palette="Set1",direction=-1)+theme_classic() +theme(axis.text.x = element_text(size = 14, face="bold",angle=0,vjust=0.5),axis.text.y = element_text(size = 12), axis.title.x = element_blank(),axis.title.y = element_text( size=16, face="bold"), legend.text=element_text(size=12),legend.title=element_blank())+ylab("Proportin of UV (%)")+ylim(0,100)+theme(legend.position = "none")

dat.n2$Disease2=c(rep("BCC",4),rep("cSCC in situ",2),rep("BCC",4))
dat.n2$Epidermis=as.numeric(dat.n2$Epidermis)
dat.n2$Tumor=as.numeric(dat.n2$Tumor)
#Figure S11b
ggplot(dat.n2, aes(x=Epidermis, y=Tumor, color=Disease2, shape=Disease2)) + geom_point(size=3)+ labs(x = "Proportion of UV (%) (skin)", y="Proportion of UV (%) (cancer)")+ theme_classic()+theme(axis.title = element_text( size=16, face="bold"),axis.text = element_text( size=13,face="bold"),legend.text=element_text(size=14))+  scale_colour_manual(name="Disease2", values= c("blue", "red"))+theme(legend.title=element_blank())+ theme(legend.position = "none")
cor.test(dat.n2$Epidermis,dat.n2$Tumor, method = "pearson")

##Comparision of Median VAFs of SNVs  
setwd("/data/Delta_data4/kysbbubbu/normal_skin/08_Counts")
dat.n=read.delim("Skin_annotation.txt",sep="\t",header=T)
dat.n2=as.data.frame(dat.n[dat.n$Tumor..analyzed=="O",c("Sample","SNV_Median_VAF","Tumor_SNV_Median_VAF")])
colnames(dat.n2)=c("Case","Epidermis","Tumor")
dat.n3 <- reshape2::melt(dat.n2, id.vars = "Case")
dat.n3$Case=factor(dat.n3$Case,levels=c("NS39","NS32","NS41","NS28","NS47","NS42","NS43","NS18","NS45","NS09"))
dat.n3$value=as.numeric(dat.n3$value)
#Figure S11c
ggplot(data=dat.n3, aes(x=Case, y=value, fill=variable))  +  geom_bar(stat="identity", position=position_dodge(),width=0.7)+ scale_fill_brewer(palette="Set1",direction=-1)+theme_classic() +theme(axis.text.x = element_text(size = 14, face="bold",angle=0,vjust=0.5),axis.text.y = element_text(size = 12), axis.title.x = element_blank(),axis.title.y = element_text( size=16, face="bold"), legend.text=element_text(size=12),legend.title=element_blank())+ylab("Median VAFs of SNVs (%)")+theme(legend.position = "none")+ylim(0,35)

dat.n2$Disease2=c(rep("BCC",4),rep("cSCC in situ",2),rep("BCC",4))
dat.n2$Epidermis=as.numeric(dat.n2$Epidermis)
dat.n2$Tumor=as.numeric(dat.n2$Tumor)
#Figure S11c
ggplot(dat.n2, aes(x=Epidermis, y=Tumor, color=Disease2, shape=Disease2)) + geom_point(size=3)+ labs(x = "Median VAFs of SNVs (%) (skin)", y="Median VAFs of SNVs (%) (cancer)")+ theme_classic()+theme(axis.title = element_text( size=16, face="bold"),axis.text = element_text( size=13,face="bold"),legend.text=element_text(size=14))+  scale_colour_manual(name="Disease2", values= c("blue", "red"))+theme(legend.title=element_blank())+ theme(legend.position = "none")
cor.test(dat.n2$Epidermis,dat.n2$Tumor, method = "pearson")


