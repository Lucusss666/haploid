setwd()
count<-count[count$class=="Gene",c(1:11,14)]
rownames(count)=gsub("[.].*","",gsub("Gohir.","Gohir_",rownames(count)))

# OG
ogQ<-read.table("./orthohomoeolog052421.txt", sep="\t", header=FALSE,stringsAsFactor = FALSE)
og<-data.frame(gorai=c(paste0(ogQ$V5,".A"),paste0(ogQ$V5,".D")),gohir=c(ogQ$V3,ogQ$V4))
og$gohir<-gsub("[.].*","",gsub("Gohir.","Gohir_",og$gohir))
og<-og[order(og$gorai),]
head(og)

# prep
combat_Expr.new <- read.delim("combat_Expr.new.txt")
count<-cbind(combat_Expr.new[og$gohir,],og)
rownames(count)=count$gorai
tpm<-cbind(hAD1merged.TPM[og$gohir,],og)
rownames(tpm)=tpm$gorai

#correlation
correlation_matrix <- cor(combat_Expr.new, method = "pearson")

# make genome list
AD1 <-list(count = count[,15:20])
A2  <-list(count = count[,4:6])
F1  <-list(count = count[,10:14])
D5  <-list(count = count[,1:3])
AD2  <-list(count = count[,21:26])
hAD1  <-list(count = count[,27:31])
hAD2  <-list(count = count[,32:37])

head(A2$count)

getTotal=function(x){y=x[grep(".A$",rownames(x)),]+x[grep(".D$",rownames(x)),];rownames(y)=gsub(".A$","",rownames(y));return(y)}
A2$count = getTotal(A2$count)
A2$tpm = getTotal(A2$tpm)
D5$count = getTotal(D5$count)
D5$tpm = getTotal(D5$tpm)

#################
## DE analysis ##
#################
library(DESeq2)
#install.packages("gplots")
library(gplots)

# count sig gene number
getSig<-function(res,fc.threshold=0,direction=NULL){
  sig<- res[res$padj<0.05 & !is.na(res$padj) & abs(res$log2FoldChange)>=fc.threshold,]
  if(is.null(direction)){
    n<-nrow(sig)
  }else if(direction=="up"){
    n<-nrow(sig[sig$log2FoldChange>0,])
  }else if(direction=="down"){
    n<-nrow(sig[sig$log2FoldChange<0,])
  }
  return(n)
}
# volcano plot
plotVolcano<-function(res, title)
{
  plot(res[,"log2FoldChange"], -log2(res[,"padj"]), main=title, xlab="log2FoldChange", ylab="-log2padj",pch=".",ylim=c(0,200))
  abline(h=-log2(0.05))
}

any(is.na(A2$count))
which(is.na(A2$count))
head(A2$count)

# A2 vs D5
count = cbind(A2$count,D5$count)
info = data.frame(sample=names(count), genome = c("A2","A2","A2","D5","D5","D5"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res <- results(DESeq(dds), contrast=c("genome","A2","D5"))
print( summary(res,alpha=.05) ) # Higher: A 3403 D 3914;
write.table(res, file="DE.A2vsD5.txt",  sep="\t")
A=res

# F1: At vs Dt
count = cbind(F1$count[grep("A$",rownames(F1$count)),], F1$count[grep("D$",rownames(F1$count)),])
rownames(count) =gsub(".A$","",rownames(count))
names(count)=paste(names(count),c("A","A","A","A","A","D","D","D","D","D"),sep=".")
info = data.frame(sample=names(count), genome = c("A","A","A","A","A","D","D","D","D","D"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
colSums(count) 
res = results(DESeq(dds),contrast=c("genome","A","D"))
print( summary(res,alpha=.05) ) # Higher: A 1793 D 1818; 
write.table(res, file="DE.F1.AtvsF1.Dt.txt",  sep="\t")
B=res

# AD1: At vs Dt
count = cbind(AD1$count[grep("A$",rownames(AD1$count)),], AD1$count[grep("D$",rownames(AD1$count)),])
rownames(count) =gsub(".A$","",rownames(count))
names(count)=paste(names(count),c("A","A","A","A","A","A","D","D","D","D","D","D"),sep=".")
info = data.frame(sample=names(count), genome = c("A","A","A","A","A","A","D","D","D","D","D","D"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
colSums(count)
res = results(DESeq(dds),contrast=c("genome","A","D"))
print( summary(res,alpha=.05) ) # Higher: A 3874 D 3982;
write.table(res, file="DE.AD1.AtvsAD1.Dt.txt",  sep="\t")
# homoeolog expression divergence in polyploid, test for Bp=0
Bp1=res

# hAD1: At vs Dt
count = cbind(hAD1$count[grep("A$",rownames(hAD1$count)),], hAD1$count[grep("D$",rownames(hAD1$count)),])
rownames(count) =gsub(".A$","",rownames(count))
names(count)=paste(names(count),c("A","A","A","A","A","D","D","D","D","D"),sep=".")
info = data.frame(sample=names(count), genome = c("A","A","A","A","A","D","D","D","D","D"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
colSums(count)
res = results(DESeq(dds),contrast=c("genome","A","D"))
print( summary(res,alpha=.05) ) # Higher: A 2310 D 2326;
write.table(res, file="DE.hAD1.AtvshAD1.Dt.txt",  sep="\t")
Hp1=res

# AD2: At vs Dt
count = cbind(AD2$count[grep("A$",rownames(AD2$count)),], AD2$count[grep("D$",rownames(AD2$count)),])
rownames(count) =gsub(".A$","",rownames(count))
names(count)=paste(names(count),c("A","A","A","A","A","A","D","D","D","D","D","D"),sep=".")
info = data.frame(sample=names(count), genome = c("A","A","A","A","A","A","D","D","D","D","D","D"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
colSums(count)
res = results(DESeq(dds),contrast=c("genome","A","D"))
print( summary(res,alpha=.05) ) # Higher: A 3779 D 3761;
write.table(res, file="DE.AD2.AtvsAD2.Dt.txt",  sep="\t")
Bp2=res

# hAD2: At vs Dt
count = cbind(hAD2$count[grep("A$",rownames(hAD2$count)),], hAD2$count[grep("D$",rownames(hAD2$count)),])
rownames(count) =gsub(".A$","",rownames(count))
names(count)=paste(names(count),c("A","A","A","A","A","A","D","D","D","D","D","D"),sep=".")
info = data.frame(sample=names(count), genome = c("A","A","A","A","A","A","D","D","D","D","D","D"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
colSums(count)
res = results(DESeq(dds),contrast=c("genome","A","D"))
print( summary(res,alpha=.05) ) # Higher: A 2969 D 3010;
write.table(res, file="DE.hAD2.AtvshAD2.Dt.txt",  sep="\t")
Hp2=res

# F1 vs AD1
count = cbind(F1$count,AD1$count)
info = data.frame(sample=names(count), genome = c("F1","F1","F1","F1","F1","AD1","AD1","AD1","AD1","AD1","AD1"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","F1","AD1"))
print( summary(res,alpha=.05) ) # higher:  F1 6593, AD1 5735;
write.table(res, file="DE.F1vsAD1.txt",  sep="\t")
W1 =res

# F1 vs AD2
count = cbind(F1$count,AD2$count)
info = data.frame(sample=names(count), genome = c("F1","F1","F1","F1","F1","AD2","AD2","AD2","AD2","AD2","AD2"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","F1","AD2"))
print( summary(res,alpha=.05) ) # higher: F1 6640, AD2 6520;
write.table(res, file="DE.AD2vsF1.txt",  sep="\t")
W2=res

# F1 vs hAD2
count = cbind(F1$count,hAD2$count)
info = data.frame(sample=names(count), genome = c("F1","F1","F1","F1","F1","hAD2","hAD2","hAD2","hAD2","hAD2","hAD2"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","F1","hAD2"))
print( summary(res,alpha=.05) ) # higher: hAD2 5576, F1 5845;
write.table(res, file="DE.F1vshAD2.txt",  sep="\t")
hAD2vsF1=res

# F1 vs hAD1
count = cbind(F1$count,hAD1$count)
info = data.frame(sample=names(count), genome = c("F1","F1","F1","F1","F1","hAD1","hAD1","hAD1","hAD1","hAD1"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","F1","hAD1"))
print( summary(res,alpha=.05) ) # higher: F1 4699, hAD1 3697
write.table(res, file="DE.F1vshAD1.txt",  sep="\t")
hAD1vsF1=res

# hAD1 vs hAD2
count = cbind(hAD1$count,hAD2$count)
info = data.frame(sample=names(count), genome = c("hAD1","hAD1","hAD1","hAD1","hAD1","hAD2","hAD2","hAD2","hAD2","hAD2","hAD2"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","hAD1","hAD2"))
print( summary(res,alpha=.05) ) # higher: hAD1 1627, hAD2 1980; 
write.table(res, file="DE.hAD1vshAD2.txt",  sep="\t")
hAD1vshAD2=res

#AD1 vs AD2
count = cbind(AD1$count,AD2$count)
info = data.frame(sample=names(count), genome = c("AD1","AD1","AD1","AD1","AD1","AD1","AD2","AD2","AD2","AD2","AD2","AD2"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","AD1","AD2"))
print( summary(res,alpha=.05) ) # higher: AD1 3266, AD2 3864; 
write.table(res, file="DE.AD1vsAD2.txt",  sep="\t")
AD1vsAD2=res

#AD1 vs hAD1
count = cbind(AD1$count,hAD1$count)
info = data.frame(sample=names(count), genome = c("AD1","AD1","AD1","AD1","AD1","AD1","hAD1","hAD1","hAD1","hAD1","hAD1"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","AD1","hAD1"))
print( summary(res,alpha=.05) ) # higher: AD1 2008, hAD1 1734; 
write.table(res, file="DE.AD1vshAD1.txt",  sep="\t")
AD1vshAD1=res

#AD2 vs hAD2
count = cbind(AD2$count,hAD2$count)
info = data.frame(sample=names(count), genome = c("AD2","AD2","AD2","AD2","AD2","AD2","hAD2","hAD2","hAD2","hAD2","hAD2","hAD2"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","AD2","hAD2"))
print( summary(res,alpha=.05) ) # higher: AD2 269, hAD2 33; 
write.table(res, file="DE.AD2vshAD2.txt",  sep="\t")
AD2vshAD2=res

# Total comparison
############################
getTotal=function(x){y=x[grep(".A$",rownames(x)),]+x[grep(".D$",rownames(x)),];rownames(y)=gsub(".A$","",rownames(y));return(y)}
F1.t = getTotal(F1$count)
AD1.t = getTotal(AD1$count)
AD2.t = getTotal(AD2$count)
hAD1.t = getTotal(hAD1$count)
hAD2.t = getTotal(hAD2$count)
libTotal =c(colSums(A2$count[,1:3]),colSums(D5$count[,1:3]))
libSize = libTotal/mean(libTotal)
mid = as.data.frame(cbind(A2$count[,1]/libSize[1]+D5$count[,1]/libSize[4], A2$count[,1]/libSize[1]+D5$count[,2]/libSize[5],A2$count[,1]/libSize[1]+D5$count[,3]/libSize[6], A2$count[,2]/libSize[2]+D5$count[,1]/libSize[4], A2$count[,2]/libSize[2]+D5$count[,2]/libSize[5],A2$count[,2]/libSize[2]+D5$count[,3]/libSize[6], A2$count[,3]/libSize[3]+D5$count[,1]/libSize[4], A2$count[,3]/libSize[3]+D5$count[,2]/libSize[5],A2$count[,3]/libSize[3]+D5$count[,3]/libSize[6]))/2
names(mid)=paste0("mid-",c("S14","S15","S16","S24","S25","S26","S34","S35","S36"))
total = cbind(A2$count,D5$count,F1.t,AD1.t,AD2.t,hAD1.t,hAD2.t,mid)
colSums(total)
write.table(total, file="total.txt",  sep="\t")
# check total grouping
library(ggplot2)
library(scales)
library(ape)
plotGrouping <- function(norm_log, color, shape, text, tip, save = "plotGrouping.pdf"){
  # norm<-sweep(total,2,info$lib_size,"/")*10^6
  # norm_log <- log2(norm+1)
  pca=prcomp(t(norm_log))
  dat = as.data.frame(pca$x)
  proportion<-summary(pca)$importance[2,1:2]*100
  proportion<-paste0(names(proportion)," (", proportion, "%)")
  p<-ggplot(aes(PC1, PC2, color=color,    shape=shape),data=dat) + geom_point() +xlab(proportion[1]) + ylab(proportion[2])
  pdf(save)
  print( p + geom_text(aes_string(x = "PC1", y = "PC2", label = text), color="grey", hjust = 0, nudge_x = 0.09) )
  
  hc<-hclust( dist(t(norm_log)) )
  tre <- as.phylo(hc)
  tre$tip.label <- as.character(tip)
  # try to match ggplot color: library(scale);
  # show_col(col4<-hue_pal()(4))
  tipCol <- color
  levels(tipCol) <- hue_pal()(nlevels(tipCol))
  plot(tre,tip.col =as.character(tipCol),  type="unrooted",cex=0.6, no.margin=TRUE)
  plot(tre,tip.col =as.character(tipCol), type="fan",cex=0.6, no.margin=TRUE)
  dev.off()   }
# plots
info = data.frame(sample=names(total), genome = c("A2","A2","A2","D5","D5","D5","F1","F1","F1","F1","F1","AD1","AD1","AD1","AD1","AD1","AD1","AD2","AD2","AD2","AD2","AD2","AD2","hAD1","hAD1","hAD1","hAD1","hAD1","hAD2","hAD2","hAD2","hAD2","hAD2","hAD2","mid","mid","mid","mid","mid","mid","mid","mid","mid"))
libTotal=colSums(total)
libSize = libTotal/mean(libTotal)
rpm = sweep(total,2,libSize,FUN="/")
info$genome <- as.character(info$genome)
info$sample <- as.character(info$sample)

if (length(unique(info$sample)) < 2) {
  stop("info$genome are less than 2")
}


color_palette <- hue_pal()(length(unique(info$genome)))
color_vector <- color_palette[as.factor(info$genome)]
plotGrouping(log2(rpm+1), color=color_vector, tip=info$sample, text=1:43, save = "plotGrouping.total.log2rpm.pdf")

# F1.total vs Mid
col_indices <- c(7:11, 35:43)
count = total[, col_indices]
info = data.frame(sample=names(count), genome = c("F1","F1","F1","F1","F1","mid","mid","mid","mid","mid","mid","mid","mid","mid"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","F1","mid"))
print( summary(res,alpha=.05) ) # higher: F1 4586 Mid 4361
write.table(res, file="DE.F1.tvsMid.txt",  sep="\t")
F1vsMid =res

# AD1.total vs Mid
col_indices <- c(12:17, 35:43)
count = total[, col_indices]
info = data.frame(sample=names(count), genome = c("AD1","AD1","AD1","AD1","AD1","AD1","mid","mid","mid","mid","mid","mid","mid","mid","mid"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","AD1","mid"))
print( summary(res,alpha=.05) ) # higher: AD1 5561, Mid 5761;
write.table(res, file="DE.AD1.tvsMid.txt",  sep="\t")
AD1vsMid =res

# AD2.total vs Mid
col_indices <- c(18:23, 35:43)
count = total[, col_indices]
info = data.frame(sample=names(count), genome = c("AD2","AD2","AD2","AD2","AD2","AD2","mid","mid","mid","mid","mid","mid","mid","mid","mid"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","AD2","mid"))
print( summary(res,alpha=.05) ) # higher: AD2 6303, Mid 6258;
write.table(res, file="DE.AD2.tvsMid.txt",  sep="\t")
AD2vsMid =res

# hAD1.total vs Mid
col_indices <- c(24:28, 35:43)
count = total[, col_indices]
info = data.frame(sample=names(count), genome = c("hAD1","hAD1","hAD1","hAD1","hAD1","mid","mid","mid","mid","mid","mid","mid","mid","mid"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","hAD1","mid"))
print( summary(res,alpha=.05) ) # higher: hAD1 5936, Mid 5790;
write.table(res, file="DE.hAD1.tvsMid.txt",  sep="\t")
hAD1vsMid =res

# hAD2.total vs Mid
col_indices <- c(29:34, 35:43)
count = total[, col_indices]
info = data.frame(sample=names(count), genome = c("hAD2","hAD2","hAD2","hAD2","hAD2","hAD2","mid","mid","mid","mid","mid","mid","mid","mid","mid"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","hAD2","mid"))
print( summary(res,alpha=.05) ) # higher: hAD2 6145, Mid 5640;
write.table(res, file="DE.hAD2.tvsMid.txt",  sep="\t")
hAD2vsMid =res

# A2 vs F1.total
col_indices <- c(7:11, 1:3)
count = total[, col_indices]
info = data.frame(sample=names(count), genome = c("F1","F1","F1","F1","F1","A2","A2","A2"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","F1","A2"))
print( summary(res,alpha=.05) ) # higher: F1 3504, A2 2673;
write.table(res, file="DE.F1.tvsA2.txt",  sep="\t")
F1vsA2 =res

#  A2 vs AD1.total
col_indices <- c(12:17, 1:3)
count = total[, col_indices]
info = data.frame(sample=names(count), genome = c("AD1","AD1","AD1","AD1","AD1","AD1","A2","A2","A2"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","AD1","A2"))
print( summary(res,alpha=.05) ) # higher AD1 3780, A2 3517;
write.table(res, file="DE.AD1.tvsA2.txt",  sep="\t")
AD1vsA2 =res

#  A2 vs AD2.total
col_indices <- c(18:23, 1:3)
count = total[, col_indices]
info = data.frame(sample=names(count), genome = c("AD2","AD2","AD2","AD2","AD2","AD2","A2","A2","A2"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","AD2","A2"))
print( summary(res,alpha=.05) ) # higher AD1 4032, A2 3538
write.table(res, file="DE.AD2.tvsA2.txt",  sep="\t")
AD2vsA2=res

#  A2 vs hAD1.total
col_indices <- c(24:28, 1:3)
count = total[, col_indices]
info = data.frame(sample=names(count), genome = c("hAD1","hAD1","hAD1","hAD1","hAD1","A2","A2","A2"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","hAD1","A2"))
print( summary(res,alpha=.05) ) # higher hAD1 3401, A2 3362
write.table(res, file="DE.hAD1vsA2.txt",  sep="\t")
hAD1vsA2 =res

#  A2 vs hAD2.total
col_indices <- c(1:3, 29:34)
count = total[, col_indices]
info = data.frame(sample=names(count), genome = c("A2","A2","A2","hAD2","hAD2","hAD2","hAD2","hAD2","hAD2"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","hAD2","A2"))
print( summary(res,alpha=.05) ) # higher A2 3061, hAD2 3650
write.table(res, file="DE.A2vshAD2.txt",  sep="\t")
hAD2vsA2=res

#  D5 vs F1.total
col_indices <- c(4:6, 7:11)
count = total[, col_indices]
info = data.frame(sample=names(count), genome = c("D5","D5","D5","F1","F1","F1","F1","F1"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","F1","D5"))
print( summary(res,alpha=.05) ) # higher D5 1813, F1 1890
write.table(res, file="DE.D5vsF1.t.txt",  sep="\t")
F1vsD5 =res

#  D5 vs AD1.total
col_indices <- c(4:6, 12:17)
count = total[, col_indices]
info = data.frame(sample=names(count), genome = c("D5","D5","D5","AD1","AD1","AD1","AD1","AD1","AD1"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","AD1","D5"))
print( summary(res,alpha=.05) ) # higher AD1 3732, D5 3764
write.table(res, file="DE.D5vsAD1.t.txt",  sep="\t")
AD1vsD5 =res

#  D5 vs AD2.total
col_indices <- c(4:6, 18:23)
count = total[, col_indices]
info = data.frame(sample=names(count), genome = c("D5","D5","D5","AD2","AD2","AD2","AD2","AD2","AD2"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","AD2","D5"))
print( summary(res,alpha=.05) ) # higher D5 4441, AD2 4569
write.table(res, file="DE.D5vsAD2.t.txt",  sep="\t")
AD2vsD5 =res

#  D5 vs hAD1.total
col_indices <- c(24:28, 4:6)
count = total[, col_indices]
info = data.frame(sample=names(count), genome = c("hAD1","hAD1","hAD1","hAD1","hAD1","D5","D5","D5"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","hAD1","D5"))
print( summary(res,alpha=.05) ) # higher  hAD1 3230, D5 3613
write.table(res, file="DE.hAD1.tvsD5.txt",  sep="\t")
hAD1vsD5 =res

#  D5 vs hAD2.total
col_indices <- c(4:6, 29:34)
count = total[, col_indices]
info = data.frame(sample=names(count), genome = c("D5","D5","D5","hAD2","hAD2","hAD2","hAD2","hAD2","hAD2"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","hAD2","D5"))
print( summary(res,alpha=.05) ) # higher D5 3588, hAD2 3867
write.table(res, file="DE.D5vshAD2.t.txt",  sep="\t")
hAD2vsD5 =res

# Mid.total vs A2
col_indices <- c(1:3, 35:43)
count = total[, col_indices]
info = data.frame(sample=names(count), genome = gsub("SD5-|-S.*","",names(count)))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","mid","A2"))
print( summary(res,alpha=.05) ) # higher Mid 4438, A2 2383
write.table(res, file="DE.MidvsA2.txt",  sep="\t")
MidvsA2 =res


# Mid.total vs D5
count = total[,grep("mid|-D5-",names(total))]
info = data.frame(sample=names(count), genome = gsub("SD5-|-S.*","",names(count)))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","mid","D5"))
print( summary(res,alpha=.05) ) # higher Mid 2755, D5 1546
write.table(res, file="DE.MidvsD5.txt",  sep="\t")
MidvsD5 =res

# F1.t vs AD1.t
col_indices <- c(7:11, 12:17)
count = total[, col_indices]
info = data.frame(sample=names(count), genome = c("F1","F1","F1","F1","F1","AD1","AD1","AD1","AD1","AD1","AD1"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","F1","AD1"))
print( summary(res,alpha=.05) ) # higher F1 3540, AD1 3047
write.table(res, file="DE.F1.tvsAD1.t.txt",  sep="\t")
F1.vsAD1.t =res

# F1.t vs AD2.t
col_indices <- c(7:11, 18:23)
count = total[, col_indices]
info = data.frame(sample=names(count), genome = c("F1","F1","F1","F1","F1","AD2","AD2","AD2","AD2","AD2","AD2"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","F1","AD2"))
print( summary(res,alpha=.05) ) # higher F1 3524, AD2 3599
write.table(res, file="DE.F1.tvsAD2.t.txt",  sep="\t")
F1.vsAD2.t =res

# F1.t vs hAD1.t
col_indices <- c(7:11, 24:28)
count = total[, col_indices]
info = data.frame(sample=names(count), genome = c("F1","F1","F1","F1","F1","hAD1","hAD1","hAD1","hAD1","hAD1"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","F1","hAD1"))
print( summary(res,alpha=.05) ) # higher F1 2468, hAD1 1973
write.table(res, file="DE.F1.tvshAD1.t.txt",  sep="\t")
F1.vshAD1.t =res

# F1.t vs hAD2.t
col_indices <- c(7:11, 29:34)
count = total[, col_indices]
info = data.frame(sample=names(count), genome = c("F1","F1","F1","F1","F1","hAD2","hAD2","hAD2","hAD2","hAD2","hAD2"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","F1","hAD2"))
print( summary(res,alpha=.05) ) # higher F1 2942, hAD2 3244
write.table(res, file="DE.F1.tvshAD2.t.txt",  sep="\t")
F1.vshAD2.t =res

# AD1.t vs hAD1.t
col_indices <- c(12:17, 24:28)
count = total[, col_indices]
info = data.frame(sample=names(count), genome = c("AD1","AD1","AD1","AD1","AD1","AD1","hAD1","hAD1","hAD1","hAD1","hAD1"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","AD1","hAD1"))
print( summary(res,alpha=.05) ) # higher AD1 1304, hAD1 1168
write.table(res, file="DE.AD1.tvshAD1.t.txt",  sep="\t")
AD1.vshAD1.t =res

# AD2.t vs hAD2.t
col_indices <- c(18:23, 29:34)
count = total[, col_indices]
info = data.frame(sample=names(count), genome = c("AD2","AD2","AD2","AD2","AD2","AD2","hAD2","hAD2","hAD2","hAD2","hAD2","hAD2"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","AD2","hAD2"))
print( summary(res,alpha=.05) ) # higher AD2 150, hAD2 18
write.table(res, file="DE.AD2.tvshAD2.t.txt",  sep="\t")
AD2.vshAD2.t =res

# save
save(list=c("A","B","Bp1","Bp2","W","W2","Hp1","Hp2", "total","mid", ls(pattern="vs")), file="diffExprPattern.rdata")

# summary
unique(rownames(A)==rownames(B)) # TRUE
unique(rownames(A)==rownames(Bp1)) # TRUE
unique(rownames(A)==rownames(Bp2)) # TRUE
unique(rownames(A)==rownames(W)) #TRUE
unique(rownames(A)==rownames(W2)) #TRUE
unique(rownames(A)==rownames(Hp1)) # TRUE
unique(rownames(A)==rownames(Hp2)) # TRUE

# make summary table
sigT<-c("sample 1","sample 2","Total","DE (q<0.05)","1>2","2>1")
for(file in list.files(pattern="^DE"))
{
  res<-read.table(file,sep="\t",header=TRUE)
  sigRes <- c(unlist(strsplit(gsub("DE.|.txt","",file),split="vs") ), nrow(res), getSig(res),getSig(res,direction="up"),getSig(res,direction="down"))
  sigT<-rbind(sigT,sigRes)
}
T<-as.data.frame(sigT[-1,], row.names=FALSE)
names(T)<-sigT[1,]
T
#pdf plot
pdf("checkDE.pdf")
# table DE
textplot(T)
mtext("DE analysis result summary")
# volcano plot
plotVolcano(A, "Diploid A2 vs D5")
plotVolcano(B, "Hybrid alleic At vs Dt")
plotVolcano(Bp1, "Polyploid homoeologous At vs Dt")
plotVolcano(Bp2, "AD2_Polyploid homoeologous At vs Dt")
plotVolcano(W, "Genome doubling F1 vs AD1")
plotVolcano(W2, "Genome doubling F1 vs AD2")
plotVolcano(Hp1, "Genome haploidization F1 vs AD1")
plotVolcano(Hp2, "Genome haploidization F1 vs AD2")
for(i in ls(pattern="vs")){plotVolcano(get(i),i)}
# compare log2 Fold change
plot( A[,"log2FoldChange"], B[,"log2FoldChange"],pch=".",main="log2 Fold Change",xlab="A - diploids",ylab="B - F1")
lines(c(-6,6),c(-6,6),col="blue")
plot( A[,"log2FoldChange"], Bp1[,"log2FoldChange"],pch=".",main="log2 Fold Change",xlab="A - diploids",ylab="Bp1 - Polyploids")
lines(c(-6,6),c(-6,6),col="blue")
plot( A[,"log2FoldChange"], Bp2[,"log2FoldChange"],pch=".",main="log2 Fold Change",xlab="A - diploids",ylab="Bp2 - Polyploids")
lines(c(-6,6),c(-6,6),col="blue")
plot( B[,"log2FoldChange"], Bp[,"log2FoldChange"],pch=".",main="log2 Fold Change",xlab="B - F1",ylab="Bp - Polyploids")
lines(c(-6,6),c(-6,6),col="blue")
plot( W[grep("A$",rownames(W)),"log2FoldChange"], W[grep("D$",rownames(W)),"log2FoldChange"],pch=".",main="log2 Fold Change",xlab="W.At",ylab="W.Dt")
lines(c(-6,6),c(-6,6),col="blue")
plot( W[grep("A$",rownames(W)),"log2FoldChange"], W[grep("D$",rownames(W)),"log2FoldChange"],pch=".",main="log2 Fold Change",xlab="W.At",ylab="W.Dt")
lines(c(-6,6),c(-6,6),col="blue")
plot( F1vsA2[,"log2FoldChange"], F1vsD5[,"log2FoldChange"],pch=".",main="log2 Fold Change")
lines(c(-6,6),c(-6,6),col="blue")
plot( AD1vsA2[,"log2FoldChange"], AD1vsD5[,"log2FoldChange"],pch=".",main="log2 Fold Change")
lines(c(-6,6),c(-6,6),col="blue")
plot( MidvsA2[,"log2FoldChange"], MidvsD5[,"log2FoldChange"],pch=".",main="log2 Fold Change")
lines(c(-6,6),c(-6,6),col="blue")

dev.off()

########################
## cis-trans analysis##
########################
load("diffExprPattern.rdata")->l
l
### comparing A-B= 0 is tricky, both are log2FoldChange and its standard error lfcse
# maybe I can compare with t test
#### T test from means and standard errors ####
# m1, m2: the sample means
# s1, s2: the sample standard deviations
# n1, n2: the same sizes
# se1, se2: the sample standard errors
# se1 <- s1/sqrt(n)
# m0: the null value for the difference in means to be tested for. Default is 0.
# equal.variance: whether or not to assume equal variance. Default is FALSE.
t.test2 <- function(m1,m2,se1,se2,n1,n2,m0=0,equal.variance=FALSE)
{
  if( equal.variance==FALSE )
  {
    # se <- sqrt( (s1^2/n1) + (s2^2/n2) )
    se <- sqrt( (se1^2) + (se2^2) )
    # welch-satterthwaite df
    # df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
    df <- ( (se1^2 + se2^2)^2 )/( (se1^2)^2/(n1-1) + (se2^2)^2/(n2-1) )
  } else
  {
    # pooled standard deviation, scaled by the sample sizes
    # se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) )
    se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) )
    df <- n1+n2-2
  }
  t <- (m1-m2-m0)/se
  dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))
  names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
  return(dat)
}

x1 = rnorm(3)
x2 = rnorm(3)
# you'll find this output agrees with that of t.test when you input x1,x2
t.test(x1,x2)
t.test2( mean(x1),  mean(x2), sd(x1)/sqrt(3), sd(x2)/sqrt(3), 3,3)

# function to make a categorization table
criteria<- as.data.frame(rbind(c("A!=0;B!=0;A=B", "1.Cis only"),
                               c("A!=0;B=0;A!=B", "2.Trans only"),
                               c("A!=0;B!=0;A!=B", "Cis+Trans"),
                               c("A=0;B!=0;A!=B", "5.Compensatory"),
                               c("A=0;B=0;A=B", "6.Conserved") ))
names(criteria) <- c("class","category")
classCisTrans<-function(A.res, B.res, A.n, B.n, log2fc.threshold=0,plotTitle=NULL)
{
  # A = log2(A2/D5), cis + trans
  A <- A.res[,c("log2FoldChange", "lfcSE", "padj")]
  names(A) <- c("A", "A.SE", "A.padj")
  # B = log2(F1_t/F1_m), cis
  B <- B.res[,c("log2FoldChange", "lfcSE", "padj")]
  names(B) <- c("B", "B.SE", "B.padj")
  
  A<-A[rownames(B),]
  table <- cbind(A,B)
  table$AminusB <- table$A - table$B
  table$AminusB.pvalue <- apply(table,1,function(x) t.test2(m1=x[1],m2=x[4],se1=x[2], se2=x[5], n1=A.n, n2=B.n)["p-value"])
  
  table$cisNtrans <- ifelse(abs(table$A)>=log2fc.threshold & table$A.padj<0.05 & !is.na(table$A.padj), "A!=0", "A=0")
  table$cis <- ifelse(abs(table$B)>=log2fc.threshold & table$B.padj<0.05 & !is.na(table$B.padj), "B!=0", "B=0")
  table$trans <- ifelse(abs(table$AminusB)>=log2fc.threshold & table$AminusB.pvalue<0.05, "A!=B", "A=B")
  table$class <- paste(table$cisNtrans,table$cis,table$trans,sep=";")
  table$category <- as.character(criteria$category[ match(table$class,criteria$class)])
  table$category[is.na(table$category)] <- "7.Ambiguous"
  table$category[ table$category=="Cis+Trans" & table$B*table$AminusB >0 ] <- "3.Cis+Trans: enhancing"
  table$category[ table$category=="Cis+Trans" & table$B*table$AminusB <0 ] <- "4.Cis+Trans: compensating"
  
  colors <- c("red","blue","purple","brown","green","black","grey")
  if(!is.null(plotTitle)){
    p<- ggplot( table, aes(x=A, y=B, color=category)) + geom_point(alpha=0.8) +  xlab("Cis + Trans") + ylab("Cis") + ggtitle(plotTitle) + scale_color_manual(values=colors) + geom_vline(xintercept=0) + geom_hline(yintercept=0) +  geom_abline(intercept = 0, slope = 1) + theme_bw()
    # p2<-ggplot( table, aes( factor(category), fill=category)) + geom_bar(stat="count") + scale_fill_manual(values=colors)
    print(p)
  }
  return(table)
}

# function to plot cis trans results
plotCisTrans<-function(table, plotTitle="")
{
  colors <- c("red","blue","purple","brown","green","black","grey")
  p<- ggplot( table, aes(x=A, y=B, color=category)) + geom_point(alpha=0.8) +  xlab("Cis + Trans") + ylab("Cis") + ggtitle(plotTitle) + scale_color_manual(values=colors) + geom_vline(xintercept=0) + geom_hline(yintercept=0) +  geom_abline(intercept = 0, slope = 1) + theme_bw()
  # p2<-ggplot( table, aes( factor(category), fill=category)) + geom_bar(stat="count") + scale_fill_manual(values=colors)
  print(p)
}

# categorization using A and B
res.F1 <- classCisTrans(A.res = A, B.res = B, A.n = 3, B.n=5, log2fc.threshold=0)
res.AD1 <- classCisTrans(A.res = A, B.res = Bp1, A.n = 3, B.n=6, log2fc.threshold=0)
res.AD2 <- classCisTrans(A.res = A, B.res = Bp2, A.n = 3, B.n=6, log2fc.threshold=0)
res.hAD1 <- classCisTrans(A.res = A, B.res = Hp1, A.n = 3, B.n=5, log2fc.threshold=0)
res.hAD2 <- classCisTrans(A.res = A, B.res = Hp2, A.n = 3, B.n=6, log2fc.threshold=0)

# make plot
library("ggplot2")
pdf("plotCistrans.pdf")
textplot(data.frame(table(res.F1$category)),cex=0.6)
textplot(data.frame(table(res.AD1$category)),cex=0.6)
textplot(data.frame(table(res.AD2$category)),cex=0.6)
textplot(data.frame(table(res.hAD1$category)),cex=0.6)
textplot(data.frame(table(res.hAD2$category)),cex=0.6)
plotCisTrans(as.data.frame(res.F1), "F1.Cis/Trans")
plotCisTrans(as.data.frame(res.AD1), "AD1.Cis/Trans")
plotCisTrans(as.data.frame(res.AD2), "AD2.Cis/Trans")
plotCisTrans(as.data.frame(res.hAD1), "hAD1.Cis/Trans")
plotCisTrans(as.data.frame(res.hAD2), "hAD2.Cis/Trans")

dev.off()

# ouput
write.table(res.F1, file ="F1.CistransregPattern.txt",sep="\t")
write.table(res.AD1, file ="AD1.CistransregPattern.txt",sep="\t")
save(res.F1, res.AD1, res.AD2, res.hAD1, res.hAD2, file="CistransregPattern.rdata")

#############################
## Genome Evolution Impact-AD1##
#############################
load("diffExprPattern.rdata")
load("regPattern.rdata")
# Categorization of additional effects according to extended analytic framework
# Hr = B - A
# Pr = Bp - A
# Wr = Bp - B

# Sr = Bp1 - Bp2
# Dr1 = Bp1 - Hp1
# Dr2 = Bp2 - Hp2
classEffects<-function(A.res, B.res, Bp.res, Hp.res, A.n, B.n, Bp.n, Hp.n, log2fc.threshold=0)
{
  # A = log2(A2/D5)
  A <- A.res[,c("log2FoldChange", "lfcSE", "padj")]
  names(A) <- c("A", "A.SE", "A.padj")
  # B = log2(F1_At/F1_Dt)
  B <- B.res[,c("log2FoldChange", "lfcSE", "padj")]
  names(B) <- c("B", "B.SE", "B.padj")
  # Bp = log2(AD1_At/AD1_Dt)
  Bp <- Bp.res[,c("log2FoldChange", "lfcSE", "padj")]
  names(Bp) <- c("Bp", "Bp.SE", "Bp.padj")
  # Hp = log2(hAD1_At/hAD1_Dt)
  Hp <- Hp.res[,c("log2FoldChange", "lfcSE", "padj")]
  names(Hp) <- c("Hp", "Hp.SE", "Hp.padj")
  
  # make sure all rownames identical
  table <- cbind(A,B,Bp,Hp)
  
  # Hr
  table$Hr <- table$B - table$A
  table$Hr.pvalue <- apply(table,1,function(x) t.test2(m1=x[4],m2=x[1],se1=x[5], se2=x[2], n1=B.n, n2=A.n)["p-value"])
  
  # Pr
  table$Pr <- table$Bp - table$A
  table$Pr.pvalue <- apply(table,1,function(x) t.test2(m1=x[7],m2=x[1],se1=x[8], se2=x[2], n1=Bp.n, n2=A.n)["p-value"])
  
  # Wr
  table$Wr <- table$Bp - table$B
  table$Wr.pvalue <- apply(table,1,function(x) t.test2(m1=x[7],m2=x[4],se1=x[8], se2=x[5], n1=Bp.n, n2=B.n)["p-value"])
  
  #Dr
  table$Dr <- table$Hp - table$Bp
  table$Dr.pvalue <- apply(table,1,function(x) t.test2(m1=x[10],m2=x[7],se1=x[11], se2=x[5], n1=Hp.n, n2=Bp.n)["p-value"])
  
  # interpretation
  table$Hr.reg <- ifelse(abs(table$Hr)>=log2fc.threshold & table$Hr.pvalue<0.05 & !is.na(table$Hr.pvalue), "Hr!=0", "Hr=0")
  table$Hr.reg[table$Hr>=log2fc.threshold & table$Hr.reg=="Hr!=0"] <- "Hr>0, stronger At"
  table$Hr.reg[table$Hr<=-log2fc.threshold & table$Hr.reg=="Hr!=0"] <- "Hr<0, stronger Dt"
  
  table$Pr.reg <- ifelse(abs(table$Pr)>=log2fc.threshold & table$Pr.pvalue<0.05 & !is.na(table$Pr.pvalue), "Pr!=0", "Pr=0")
  table$Pr.reg[table$Pr>=log2fc.threshold & table$Pr.reg=="Pr!=0"] <- "Pr>0, stronger At"
  table$Pr.reg[table$Pr<=-log2fc.threshold & table$Pr.reg=="Pr!=0"] <- "Pr<0, stronger Dt"
  
  table$Wr.reg <- ifelse(abs(table$Wr)>=log2fc.threshold & table$Wr.pvalue<0.05 & !is.na(table$Wr.pvalue), "Wr!=0", "Wr=0")
  table$Wr.reg[table$Wr>=log2fc.threshold & table$Wr.reg=="Wr!=0"] <- "Wr>0, stronger At"
  table$Wr.reg[table$Wr<=-log2fc.threshold & table$Wr.reg=="Wr!=0"] <- "Wr<0, stronger Dt"
  
  table$Dr.reg <- ifelse(abs(table$Dr)>=log2fc.threshold & table$Dr.pvalue<0.05 & !is.na(table$Dr.pvalue), "Dr!=0", "Dr=0")
  table$Dr.reg[table$Dr>=log2fc.threshold & table$Dr.reg=="Dr!=0"] <- "Dr>0, stronger At"
  table$Dr.reg[table$Dr<=-log2fc.threshold & table$Dr.reg=="Dr!=0"] <- "Dr<0, stronger Dt"
  
  return(table)
}

# categorization using A, B, Bp, Hp
effect = classEffects(A.res=A, B.res=B, Bp.res=Bp1, Hp.res=Hp1, A.n=3, B.n=5, Bp.n=6, Hp.n=5, log2fc.threshold=0)
unique(effect[,1:6]==res[,1:6]) 
res=cbind(res, effect[,-c(1:6)])

# make plots
sumT =  rbind(c("Regulation Pattern","Measure","A","D"),
              c("Diploid divergence, A2vsD5", "A", getSig(A,direction="up"), getSig(A,direction="down")),
              c("Homoeolog bias, F1", "B", getSig(B,direction="up"), getSig(B,direction="down")),
              c("Homoeolog bias, AD1", "Bp1", getSig(Bp1,direction="up"), getSig(Bp1,direction="down")),
              c("Hybridization effect direction", "Hr", length(grep("Hr>0",res$Hr.reg)), length(grep("Hr<0",res$Hr.reg))),
              c("Allopolyploidy effect direction", "Pr1", length(grep("Pr>0",res$Pr.reg)), length(grep("Pr<0",res$Pr.reg))),
              c("Genome doubling effect direction", "Wr1", length(grep("Wr>0",res$Wr.reg)), length(grep("Wr<0",res$Wr.reg))),
              c("Hp-Bp", "Dr1", length(grep("Dr>0",res$Dr.reg)),length(grep("Dr<0",res$Dr.reg)))
)
sumT
[,1]                               [,2]      [,3]   [,4]  
[1,] "Regulation Pattern"               "Measure" "A"    "D"   
[2,] "Diploid divergence, A2vsD5"       "A"       "3403" "3914"
[3,] "Homoeolog bias, F1"               "B"       "1793" "1818"
[4,] "Homoeolog bias, AD1"              "Bp1"     "3874" "3982"
[5,] "Hybridization effect direction"   "Hr"      "1619" "1004"
[6,] "Allopolyploidy effect direction"  "Pr1"     "2603" "2164"
[7,] "Genome doubling effect direction" "Wr1"     "1706" "1692"
[8,] "Hp-Bp"                            "Dr1"     "76"   "102" 

pdf("plotEvoImpactAD1.pdf")
textplot(data.frame(sumT),cex=0.6)
# compare impacts
plot( res[,"Hr"], res[,"Pr"],pch=".", main=paste("cor = ", cor(res[,"Hr"],res[,"Pr"],use="complete.obs")))
lines(c(-10,10),c(-10,10),col="blue")
plot( res[,"Hr"], res[,"Wr"],pch=".", main=paste("cor = ", cor(res[,"Hr"],res[,"Wr"],use="complete.obs")))
lines(c(-10,10),c(-10,10),col="blue")
plot( res[,"Wr"], res[,"Pr"],pch=".", main=paste("cor = ", cor(res[,"Wr"],res[,"Pr"],use="complete.obs")))
lines(c(-10,10),c(-10,10),col="blue")
plot( res[,"Dr"], res[,"Hr"],pch=".", main=paste("cor = ", cor(res[,"Dr"],res[,"Hr"],use="complete.obs")))
lines(c(-10,10),c(-10,10),col="blue")
plot( res[,"Dr"], res[,"Pr"],pch=".", main=paste("cor = ", cor(res[,"Dr"],res[,"Pr"],use="complete.obs")))
lines(c(-10,10),c(-10,10),col="blue")
plot( res[,"Dr"], res[,"Wr"],pch=".", main=paste("cor = ", cor(res[,"Dr"],res[,"Wr"],use="complete.obs")))
lines(c(-10,10),c(-10,10),col="blue")
dev.off()

# ouput
write.table(res, file ="AD1regPattern.txt",sep="\t")
save(res, file="AD1regPattern.rdata")

#############################
## Genome Evolution Impact-AD2##
#############################
load("diffExprPattern.rdata")
load("regPattern.rdata")
# Categorization of additional effects according to extended analytic framework
# Hr = B - A
# Pr = Bp - A
# Wr = Bp - B

# Sr = Bp1 - Bp2
# Dr1 = Bp1 - Hp1
# Dr2 = Bp2 - Hp2
classEffects<-function(A.res, B.res, Bp.res, Hp.res, A.n, B.n, Bp.n, Hp.n, log2fc.threshold=0)
{
    # A = log2(A2/D5)
    A <- A.res[,c("log2FoldChange", "lfcSE", "padj")]
    names(A) <- c("A", "A.SE", "A.padj")
    # B = log2(F1_At/F1_Dt)
    B <- B.res[,c("log2FoldChange", "lfcSE", "padj")]
    names(B) <- c("B", "B.SE", "B.padj")
    # Bp = log2(AD1_At/AD1_Dt)
    Bp <- Bp.res[,c("log2FoldChange", "lfcSE", "padj")]
    names(Bp) <- c("Bp", "Bp.SE", "Bp.padj")
    # Hp = log2(hAD1_At/hAD1_Dt)
    Hp <- Hp.res[,c("log2FoldChange", "lfcSE", "padj")]
    names(Hp) <- c("Hp", "Hp.SE", "Hp.padj")
	
    # make sure all rownames identical
    table <- cbind(A,B,Bp,Hp)
  
    # Hr
    table$Hr <- table$B - table$A
    table$Hr.pvalue <- apply(table,1,function(x) t.test2(m1=x[4],m2=x[1],se1=x[5], se2=x[2], n1=B.n, n2=A.n)["p-value"])
  
    # Pr
    table$Pr <- table$Bp - table$A
    table$Pr.pvalue <- apply(table,1,function(x) t.test2(m1=x[7],m2=x[1],se1=x[8], se2=x[2], n1=Bp.n, n2=A.n)["p-value"])
  
    # Wr
    table$Wr <- table$Bp - table$B
    table$Wr.pvalue <- apply(table,1,function(x) t.test2(m1=x[7],m2=x[4],se1=x[8], se2=x[5], n1=Bp.n, n2=B.n)["p-value"])
    
    #Dr
    table$Dr <- table$Hp - table$Bp
    table$Dr.pvalue <- apply(table,1,function(x) t.test2(m1=x[10],m2=x[7],se1=x[11], se2=x[5], n1=Hp.n, n2=Bp.n)["p-value"])
  
    # interpretation
    table$Hr.reg <- ifelse(abs(table$Hr)>=log2fc.threshold & table$Hr.pvalue<0.05 & !is.na(table$Hr.pvalue), "Hr!=0", "Hr=0")
    table$Hr.reg[table$Hr>=log2fc.threshold & table$Hr.reg=="Hr!=0"] <- "Hr>0, stronger At"
    table$Hr.reg[table$Hr<=-log2fc.threshold & table$Hr.reg=="Hr!=0"] <- "Hr<0, stronger Dt"
  
    table$Pr.reg <- ifelse(abs(table$Pr)>=log2fc.threshold & table$Pr.pvalue<0.05 & !is.na(table$Pr.pvalue), "Pr!=0", "Pr=0")
    table$Pr.reg[table$Pr>=log2fc.threshold & table$Pr.reg=="Pr!=0"] <- "Pr>0, stronger At"
    table$Pr.reg[table$Pr<=-log2fc.threshold & table$Pr.reg=="Pr!=0"] <- "Pr<0, stronger Dt"
  
    table$Wr.reg <- ifelse(abs(table$Wr)>=log2fc.threshold & table$Wr.pvalue<0.05 & !is.na(table$Wr.pvalue), "Wr!=0", "Wr=0")
    table$Wr.reg[table$Wr>=log2fc.threshold & table$Wr.reg=="Wr!=0"] <- "Wr>0, stronger At"
    table$Wr.reg[table$Wr<=-log2fc.threshold & table$Wr.reg=="Wr!=0"] <- "Wr<0, stronger Dt"
    
    table$Dr.reg <- ifelse(abs(table$Dr)>=log2fc.threshold & table$Dr.pvalue<0.05 & !is.na(table$Dr.pvalue), "Dr!=0", "Dr=0")
    table$Dr.reg[table$Dr>=log2fc.threshold & table$Dr.reg=="Dr!=0"] <- "Dr>0, stronger At"
    table$Dr.reg[table$Dr<=-log2fc.threshold & table$Dr.reg=="Dr!=0"] <- "Dr<0, stronger Dt"
	
    return(table)
}

# categorization using A, B, Bp, Hp
effect = classEffects(A.res=A, B.res=B, Bp.res=Bp2, Hp.res=Hp2, A.n=3, B.n=5, Bp.n=6, Hp.n=6, log2fc.threshold=0)
unique(effect[,1:6]==res[,1:6])
res=cbind(res, effect[,-c(1:6)])

# make plots
sumT =  rbind(c("Regulation Pattern","Measure","A","D"),
              c("Diploid divergence, A2vsD5", "A", getSig(A,direction="up"), getSig(A,direction="down")),
              c("Homoeolog bias, F1", "B", getSig(B,direction="up"), getSig(B,direction="down")),
              c("Homoeolog bias, AD2", "Bp", getSig(Bp2,direction="up"), getSig(Bp2,direction="down")),
              c("Hybridization effect direction", "Hr", length(grep("Hr>0",res$Hr.reg)), length(grep("Hr<0",res$Hr.reg))),
              c("Allopolyploidy effect direction", "Pr2", length(grep("Pr>0",res$Pr.reg)), length(grep("Pr<0",res$Pr.reg))),
              c("Genome doubling effect direction", "Wr2", length(grep("Wr>0",res$Wr.reg)), length(grep("Wr<0",res$Wr.reg))),
              c("Hp-Bp", "Dr2", length(grep("Dr>0",res$Dr.reg)),length(grep("Dr<0",res$Dr.reg)))
)
sumT
[,1]                               [,2]      [,3]   [,4]  
[1,] "Regulation Pattern"               "Measure" "A"    "D"   
[2,] "Diploid divergence, A2vsD5"       "A"       "3403" "3914"
[3,] "Homoeolog bias, F1"               "B"       "1793" "1818"
[4,] "Homoeolog bias, AD2"              "Bp"      "3779" "3761"
[5,] "Hybridization effect direction"   "Hr"      "1619" "1004"
[6,] "Allopolyploidy effect direction"  "Pr2"     "2665" "2270"
[7,] "Genome doubling effect direction" "Wr2"     "1626" "1653"
[8,] "Hp-Bp"                            "Dr2"     "43"   "39" 

pdf("plotEvoImpactAD2.pdf")
textplot(data.frame(sumT),cex=0.6)
# compare impacts
plot( res[,"Hr"], res[,"Pr"],pch=".", main=paste("cor = ", cor(res[,"Hr"],res[,"Pr"],use="complete.obs")))
lines(c(-10,10),c(-10,10),col="blue")
plot( res[,"Hr"], res[,"Wr"],pch=".", main=paste("cor = ", cor(res[,"Hr"],res[,"Wr"],use="complete.obs")))
lines(c(-10,10),c(-10,10),col="blue")
plot( res[,"Wr"], res[,"Pr"],pch=".", main=paste("cor = ", cor(res[,"Wr"],res[,"Pr"],use="complete.obs")))
lines(c(-10,10),c(-10,10),col="blue")
plot( res[,"Dr"], res[,"Hr"],pch=".", main=paste("cor = ", cor(res[,"Dr"],res[,"Hr"],use="complete.obs")))
lines(c(-10,10),c(-10,10),col="blue")
plot( res[,"Dr"], res[,"Pr"],pch=".", main=paste("cor = ", cor(res[,"Dr"],res[,"Pr"],use="complete.obs")))
lines(c(-10,10),c(-10,10),col="blue")
plot( res[,"Dr"], res[,"Wr"],pch=".", main=paste("cor = ", cor(res[,"Dr"],res[,"Wr"],use="complete.obs")))
lines(c(-10,10),c(-10,10),col="blue")
dev.off()

# ouput
write.table(res, file ="EvoImpactregPattern.txt",sep="\t")
save(res, file="EvoImpactregPattern.rdata")

###################################
## Expression Dominance Analysis ##
###################################
load("diffExprPattern.rdata")
load("regPattern.rdata")
# F1vsMid
# F1vsA2
# F1vsD5
# hAD1vsMid
# hAD1vsA2
# hAD1vsD5
# hAD1vsF1
# hAD2vsMid
# hAD2vsA2
# hAD2vsD5
# hAD2vsF1
classDominance<-function(TvsMid, TvsA2, TvsD5, log2fc.threshold=0)
{
  # Hybrid/polyploid vs Mid parental val
  TvsMid <- data.frame(TvsMid[,c("log2FoldChange", "lfcSE", "padj")])
  names(TvsMid) <- c("TvsMid", "TvsMid.SE", "TvsMid.padj")
  # Hybrid/polyploid vs Parent 1
  TvsA <- data.frame(TvsA2[,c("log2FoldChange", "lfcSE", "padj")])
  names(TvsA) <- c("TvsA", "TvsA.SE", "TvsA.padj")
  # Hybrid/polyploid vs Parent 2
  TvsD <- data.frame(TvsD5[,c("log2FoldChange", "lfcSE", "padj")])
  names(TvsD) <- c("TvsD", "TvsD.SE", "TvsD.padj")
  
  tbl = cbind(TvsMid, TvsA, TvsD)
  
  tbl$TvsMid[is.na(tbl$TvsMid)] =0
  tbl$TvsA[is.na(tbl$TvsA)] =0
  tbl$TvsD[is.na(tbl$TvsD)] =0
  
  # judge
  tbl$additivity <- ifelse(abs(tbl$TvsMid)>=log2fc.threshold & tbl$TvsMid.padj<0.05 & !is.na(tbl$TvsMid.padj), "T!=Mid", "T=Mid")
  tbl$TvsA.reg <- ifelse(abs(tbl$TvsA)>=log2fc.threshold & tbl$TvsA.padj<0.05 & !is.na(tbl$TvsA.padj), "T!=A", "T=A")
  tbl$TvsD.reg <- ifelse(abs(tbl$TvsD)>=log2fc.threshold & tbl$TvsD.padj<0.05 & !is.na(tbl$TvsD.padj), "T!=D", "T=D")
  
  # together
  tbl$class <- paste(tbl$additivity, tbl$TvsA.reg, tbl$TvsD.reg, sep=";")
  
  # assign category
  tbl$category = "Other non-additivity"
  tbl$category[grep("T=Mid",tbl$class)] = "Additivity"
  tbl$category[grep("T!=Mid;T=A;T!=D",tbl$class)] = "A-dominant"
  tbl$category[grep("T!=Mid;T!=A;T=D",tbl$class)] = "D-dominant"
  tbl$category[grepl("T!=Mid;T!=A;T!=D",tbl$class) & tbl$TvsA>0 & tbl$TvsD>0] = "Transgressive Up"
  tbl$category[grepl("T!=Mid;T!=A;T!=D",tbl$class) & tbl$TvsA<0 & tbl$TvsD<0] = "Transgressive Down"
  
  return(tbl)
}

# categorization using A, B and Bp
dominance.F1 = classDominance(TvsMid=F1vsMid, TvsA2=F1vsA2, TvsD5=F1vsD5, log2fc.threshold=0)
dominance.AD1 = classDominance(TvsMid=AD1vsMid, TvsA2=AD1vsA2, TvsD5=AD1vsD5, log2fc.threshold=0)
dominance.AD2 = classDominance(TvsMid=AD2vsMid, TvsA2=AD2vsA2, TvsD5=AD2vsD5, log2fc.threshold=0)
dominance.hAD1 = classDominance(TvsMid=hAD1vsMid, TvsA2=hAD1vsA2, TvsD5=hAD1vsD5, log2fc.threshold=0)
dominance.hAD2 = classDominance(TvsMid=hAD2vsMid, TvsA2=hAD2vsA2, TvsD5=hAD2vsD5, log2fc.threshold=0)

# make plots
sumT =  rbind(table(dominance.AD1$category), table(dominance.F1$category),table(dominance.AD2$category),table(dominance.hAD1$category),table(dominance.hAD2$category))
rownames(sumT)=c("AD1","F1","AD2","hAD1","hAD2")
pdf("plotDominance.pdf")
textplot(data.frame(sumT),cex=0.6)
dev.off()
sumT
#A-dominant Additivity D-dominant Other non-additivity Transgressive Down Transgressive Up
#AD1        3398      11567       2967                 2322               1174             1461
#F1         1476      13942       3172                 2934                569              796
#AD2        3930      10328       2545                 2366               1620             2100
#hAD1       2792      11163       2564                 3148               1774             1448
#hAD2       3332      11104       2545                 2767               1388             1753

# merge
res$dominance.AD1=dominance.AD1$category
res$dominance.F1 =dominance.F1$category
res$dominance.AD2 =dominance.AD2$category
res$dominance.hAD1 =dominance.hAD1$category
res$dominance.hAD2 =dominance.hAD2$category

# ouput
write.table(res, file ="ELDregPattern.txt",sep="\t")
save(res,dominance.AD1,dominance.F1,dominance.AD2,dominance.hAD1,dominance.hAD2, file="ELDregPattern.rdata")

###################################
## Summarize Regulatory Patterns ##
###################################
load("regPattern.rdata")

sumTbl =as.data.frame(rbind(
  c("Diploid divergence, higher expresson", "A", length(which(res$A>0 & res$A.padj<0.05 & !is.na(res$A.padj))), length(which(res$A<0 & res$A.padj<0.05 & !is.na(res$A.padj)))),
  c("Homoeolog bias in F1, higher expression", "B", length(which(res$B>0 & res$B.padj<0.05 & !is.na(res$B.padj))), length(which(res$B<0 & res$B.padj<0.05 & !is.na(res$B.padj)))),
  c("Homoeolog bias in AD1, higher expression", "Bp", length(which(res$Bp>0 & res$Bp.padj<0.05 & !is.na(res$Bp.padj))), length(which(res$Bp<0 & res$Bp.padj<0.05 & !is.na(res$Bp.padj)))),
  c("F1 total expression vs diploids, DE" ,"F1vsA2 & F1vsD5", length(which(dominance.F1$TvsA.reg =="T!=A")), length(which(dominance.F1$TvsD.reg =="T!=D"))),
  c("AD1 total expression vs diploids, DE" ,"AD1vsA2 & AD1vsD5", length(which(dominance.AD1$TvsA.reg =="T!=A")), length(which(dominance.AD1$TvsD.reg =="T!=D"))),
  c("Hybridization effect, up regulate", "Hr", length(grep("Hr>0",res$Hr.reg)), length(grep("Hr<0",res$Hr.reg))),
  c("Allopolyploidy effect, up regulate", "Pr", length(grep("Pr>0",res$Pr.reg)), length(grep("Pr<0",res$Pr.reg))),
  c("Genome doubling effect, up regulate", "Wr", length(grep("Wr>0",res$Wr.reg)), length(grep("Wr<0",res$Wr.reg)))
))
names(sumTbl)=c("Regulation Pattern","Measure","A","D")
sumTbl[,c("A","D")]=apply(sumTbl[,c("A","D")],2,as.numeric)
sumTbl$balance=ifelse(sumTbl$A>sumTbl$D,"A>D","A<D")
sumTbl$chisqTest.pval = round(apply(sumTbl[,c("A","D")],1,function(x)chisq.test(x)$"p.value") ,6)
sumTbl
library(gplots)
pdf("plotRegSummary.pdf")
textplot(sumTbl,cex=0.6)
dev.off()
