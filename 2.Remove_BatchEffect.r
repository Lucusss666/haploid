BiocManager::install("sva")
count <- A2D5.counts
info = data.frame(sample=names(A2D5.counts), batch = c("batch1","batch1","batch1","batch1","batch1","batch1","batch1","batch1","batch1","batch2","batch2","batch2","batch1","batch1","batch2","batch2","batch2","batch1","batch1","batch1","batch2","batch2","batch2",
                                                "batch1","batch1","batch1","batch2","batch2","batch2","batch1","batch1","batch2","batch2","batch2","batch1","batch1","batch1"))
info$type = as.factor = c("D5","D5","D5","A2","A2","A2","A2D5","A2D5","A2D5","F1","F1","F1","F1","F1","AD1","AD1","AD1","AD1","AD1","AD1","AD2","AD2","AD2","AD2","AD2","AD2","hAD1","hAD1","hAD1","hAD1","hAD1","hAD2","hAD2","hAD2","hAD2","hAD2","hAD2")

#install.packages('FactoMineR')
#install.packages('factoextra')
#install.packages('ggplot2')
library(FactoMineR)
library(ggplot2)
library(factoextra)

# Before
pre.pca <- PCA(t(Expr),graph = FALSE)
ind <- get_pca_ind(pre.pca)
fviz_pca_ind(pre.pca,
             geom= "point",
             col.ind = info$batch,
             addEllipses = TRUE,
             legend.title="Group"  )

# Remove batch effect
model <- model.matrix(~as.factor(info$type))
library(sva)
expr_count_combat <- ComBat_seq(counts = as.matrix(count), 
                                covar_mod = model,
                                batch = info$batch,
                                group = info$type)
expr_count_combat[1:4,1:4]

# After
af.pca <- PCA(t(expr_count_combat),graph = FALSE)
ind <- get_pca_ind(af.pca)
fviz_pca_ind(af.pca,
             geom= "point",
             col.ind = info$batch,
             addEllipses = TRUE,
             legend.title="Group"  )


write.table(expr_count_combat, file="combat_Expr.new.txt",  sep="\t")