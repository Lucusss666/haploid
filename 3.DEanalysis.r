### DE-analysis
library(DESeq2)
combat_Expr.new <- read.delim("combat_Expr.new.txt", row.names=1)
AD1 <-list(count = combat_Expr.new[,15:20])
A2  <-list(count = combat_Expr.new[,4:6])
F1  <-list(count = combat_Expr.new[,10:14])
D5  <-list(count = combat_Expr.new[,1:3])
A2D5 <- list(count = combat_Expr.new[,7:9])
AD2  <-list(count = combat_Expr.new[,21:26])
hAD1  <-list(count = combat_Expr.new[,27:31])
hAD2  <-list(count = combat_Expr.new[,32:37])

# F1 vs AD1
count = cbind(F1$count,AD1$count)
info = data.frame(sample = names(count), genome = c("F1","F1","F1","F1","F1","AD1","AD1","AD1","AD1","AD1","AD1"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","F1","AD1"))
print( summary(res,alpha=.05) ) # higher: AD1 11232, F1 11561;
write.table(res, file="DE.F1vsAD1.txt",  sep="\t")
diff_genes <- subset(res, padj < 0.05)
write.table(diff_genes, file = "filtered.F1vsAD1.txt",  sep="\t")
F1vsAD1 =res

# F1 vs AD2
count = cbind(F1$count,AD2$count)
info = data.frame(sample=names(count), genome = c("F1","F1","F1","F1","F1","AD2","AD2","AD2","AD2","AD2","AD2"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","F1","AD2"))
print( summary(res,alpha=.05) ) # higher: AD2 11855, F1 11367;
write.table(res, file="DE.F1vsAD2.txt",  sep="\t")
diff_genes <- subset(res, padj < 0.05)
write.table(diff_genes, file = "filtered.F1vsAD2.txt",  sep="\t")
F1vsAD2=res

# F1 vs hAD2
count = cbind(F1$count,hAD2$count)
info = data.frame(sample=names(count), genome = c("F1","F1","F1","F1","F1","hAD2","hAD2","hAD2","hAD2","hAD2","hAD2"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","F1","hAD2"))
print( summary(res,alpha=.05) ) # higher: F1 9917, hAD2 10978;
write.table(res, file="DE.F1vshAD2.txt",  sep="\t")
diff_genes <- subset(res, padj < 0.05)
write.table(diff_genes, file = "filtered.F1vshAD2.txt",  sep="\t")
F1vshAD2=res

# F1 vs hAD1
count = cbind(F1$count,hAD1$count)
info = data.frame(sample=names(count), genome = c("F1","F1","F1","F1","F1","hAD1","hAD1","hAD1","hAD1","hAD1"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","F1","hAD1"))
print( summary(res,alpha=.05) ) # higher: F1 8493, hAD1 8297
write.table(res, file="DE.F1vshAD1.txt",  sep="\t")
diff_genes <- subset(res, padj < 0.05)
write.table(diff_genes, file = "filtered.F1vshAD1.txt",  sep="\t")
F1vshAD1=res

# hAD1 vs hAD2
count = cbind(hAD1$count,hAD2$count)
info = data.frame(sample=names(count), genome = c("hAD1","hAD1","hAD1","hAD1","hAD1","hAD2","hAD2","hAD2","hAD2","hAD2","hAD2"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","hAD1","hAD2"))
print( summary(res,alpha=.05) )  # higher: hAD1 4322, hAD2 4425;
write.table(res, file="DE.hAD1vshAD2.txt",  sep="\t")
diff_genes <- subset(res, padj < 0.05)
write.table(diff_genes, file = "filtered.hAD1vshAD2.txt",  sep="\t")
hAD1vshAD2=res

#AD1 vs AD2
count = cbind(AD1$count,AD2$count)
info = data.frame(sample=names(count), genome = c("AD1","AD1","AD1","AD1","AD1","AD1","AD2","AD2","AD2","AD2","AD2","AD2"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","AD1","AD2"))
print( summary(res,alpha=.05) ) # higher: AD1 7179, AD2 7902;
write.table(res, file="DE.AD1vsAD2.txt",  sep="\t")
diff_genes <- subset(res, padj < 0.05)
write.table(diff_genes, file = "filtered.AD1vsAD2.txt",  sep="\t")
AD1vsAD2=res

#AD1 vs hAD1
count = cbind(AD1$count,hAD1$count)
info = data.frame(sample=names(count), genome = c("AD1","AD1","AD1","AD1","AD1","AD1","hAD1","hAD1","hAD1","hAD1","hAD1"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","AD1","hAD1"))
print( summary(res,alpha=.05) )# higher: AD1 2682, hAD1 2964;
write.table(res, file="DE.AD1vshAD1.txt",  sep="\t")
diff_genes <- subset(res, padj < 0.05)
write.table(diff_genes, file = "filtered.AD1vshAD1.txt",  sep="\t")
AD1vshAD1=res

#AD2 vs hAD2
count = cbind(AD2$count,hAD2$count)
info = data.frame(sample=names(count), genome = c("AD2","AD2","AD2","AD2","AD2","AD2","hAD2","hAD2","hAD2","hAD2","hAD2","hAD2"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","AD2","hAD2"))
print( summary(res,alpha=.05) ) # higher: AD2 327, hAD2 59;
write.table(res, file="DE.AD2vshAD2.txt",  sep="\t")
diff_genes <- subset(res, padj < 0.05)
write.table(diff_genes, file = "filtered.AD2vshAD2.txt",  sep="\t")
AD2vshAD2=res

#A2D5 vs F1
count = cbind(A2D5$count,F1$count)
info = data.frame(sample=names(count), genome = c("A2D5","A2D5","A2D5","F1","F1","F1","F1","F1"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","A2D5","F1"))
print( summary(res,alpha=.05) ) # higher: A2D5 3435, F1 3587;
write.table(res, file="DE.A2D5vsF1.txt",  sep="\t")
diff_genes <- subset(res, padj < 0.05)
write.table(diff_genes, file = "filtered.A2D5vsF1.txt",  sep="\t")
A2D5vsF1=res

#A2D5 vs hAD1
count = cbind(A2D5$count,hAD1$count)
info = data.frame(sample=names(count), genome = c("A2D5","A2D5","A2D5","hAD1","hAD1","hAD1","hAD1","hAD1"))
info$genome <- factor(info$genome)
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","A2D5","hAD1"))
print( summary(res,alpha=.05) ) # higher: A2D5 10827, hAD1 10259;
write.table(res, file="DE.A2D5vshAD1.txt",  sep="\t")
diff_genes <- subset(res, padj < 0.05)
write.table(diff_genes, file = "filtered.A2D5vshAD1.txt",  sep="\t")
A2D5vshAD1=res

#A2D5 vs hAD2
count = cbind(A2D5$count,hAD2$count)
info = data.frame(sample=names(count), genome = c("A2D5","A2D5","A2D5","hAD2","hAD2","hAD2","hAD2","hAD2","hAD2"))
info$genome <- factor(info$genome)
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","A2D5","hAD2"))
print( summary(res,alpha=.05) ) # higher: A2D5 10248, hAD2 10965;
write.table(res, file="DE.A2D5vshAD2.txt",  sep="\t")
diff_genes <- subset(res, padj < 0.05)
write.table(diff_genes, file = "filtered.A2D5vshAD2.txt",  sep="\t")
A2D5vshAD2=res

#A2D5 vs AD1
count = cbind(A2D5$count,AD1$count)
info = data.frame(sample=names(count), genome = c("A2D5","A2D5","A2D5","AD1","AD1","AD1","AD1","AD1","AD1"))
info$genome <- factor(info$genome)
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","A2D5","AD1"))
print( summary(res,alpha=.05) ) # higher: A2D5 10986, AD1 11326;
write.table(res, file="DE.A2D5vsAD1.txt",  sep="\t")
#diff_genes <- subset(res, padj < 0.05)
#write.table(diff_genes, file = "filtered.A2D5vshAD1.txt",  sep="\t")
#A2D5vshAD1=res

#A2D5 vs AD2
count = cbind(A2D5$count,AD2$count)
info = data.frame(sample=names(count), genome = c("A2D5","A2D5","A2D5","AD2","AD2","AD2","AD2","AD2","AD2"))
info$genome <- factor(info$genome)
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","A2D5","AD2"))
print( summary(res,alpha=.05) ) # higher: A2D5 11870, AD2 12552;
write.table(res, file="DE.A2D5vsAD2.txt",  sep="\t")
#diff_genes <- subset(res, padj < 0.05)
#write.table(diff_genes, file = "filtered.A2D5vshAD2.txt",  sep="\t")
#A2D5vshAD2=res

# Muti-factor
library(installr)
library(knitr)
library(dplyr)

count = cbind(AD1$count,hAD1$count,AD2$count,hAD2$count)
info = data.frame(sample=names(count), species = c("AD1","AD1","AD1","AD1","AD1","AD1","AD1","AD1","AD1","AD1","AD1","AD2","AD2","AD2","AD2","AD2","AD2","AD2","AD2","AD2","AD2","AD2","AD2"),
                  ploidy_level = c(4,4,4,4,4,4,2,2,2,2,2,4,4,4,4,4,4,2,2,2,2,2,2))
info$species <- factor(info$species)
info$ploidy_level <- factor(info$ploidy_level)
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design =  ~ species + ploidy_level + species:ploidy_level)
dds$ploidy_level = relevel( dds$ploidy_level, "4")
dds$species = relevel( dds$species, "AD1")
dds$ploidy_level
dds$species
dds <- DESeq(dds)
resultsNames(dds)

#Gh vs Gb
res = results(dds, contrast=list(c("species_AD2_vs_AD1")))
print( summary(res,alpha=.05) ) # higher: AD1 5608, AD2 5387;
write.table(res, file="DE.GhvshGb.txt",  sep="\t")
diff_genes <- subset(res, padj < 0.05)
write.table(diff_genes, file = "filtered.GhvshGb.txt",  sep="\t")
GhvshGb=res

# allopolyploids  vs haploids
res = results(dds, contrast=list(c("ploidy_level_2_vs_4")))
print( summary(res,alpha=.05) ) # higher: 4 2611, 2 1939;
write.table(res, file="DE.polyploidvshaploid.txt",  sep="\t")
diff_genes <- subset(res, padj < 0.05)
write.table(diff_genes, file = "filtered.polyploidvshaploid.txt",  sep="\t")
polyploidvshaploid=res

# species-specific
res = results(dds, contrast=list(c("speciesAD2.ploidy_level2")))
print( summary(res,alpha=.05) ) # up: 83, down: 232 ;
write.table(res, file="DE.polyploidvshaploid.txt",  sep="\t")
diff_genes <- subset(res, padj < 0.05)
write.table(diff_genes, file = "filtered.polyploidvshaploid.txt",  sep="\t")
polyploidvshaploid=res