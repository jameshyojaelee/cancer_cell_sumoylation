#Pathway Analysis

library(pathview)
library(gage)
library(gageData)
library("AnnotationDbi")
library("org.Hs.eg.db")

columns(org.Hs.eg.db)


col_AUC <- read.csv("C:/Users/james/Desktop/Gene KO/project/inhibition/colorectal_AUC_mut.csv", fileEncoding="UTF-8-BOM")
col_IC50 <- read.csv("C:/Users/james/Desktop/Gene KO/project/inhibition/colorectal_IC50_mut.csv", fileEncoding="UTF-8-BOM")
head(col_AUC)
head(IC50)


col_path_df <- col_AUC
col_path_df[] <- lapply(col_AUC, as.character)

class(col_path_df$Hugo_Symbol)

head(col_path_df)

col_path_df$symbol <- mapIds(org.Hs.eg.db,
                     keys=col_path_df$Hugo_Symbol, # Our genenames
                     keytype="SYMBOL",        # The format of our genenames
                     column="PATH",          # The new format we want to add
                     multiVals="first")


mapIds(org.Hs.eg.db,
       keys="HUWE1", # Our genenames
       keytype="ALIAS",        # The format of our genenames
       column="PATH",          # The new format we want to add
       multiVals="first")

head(path_df)

data(kegg.sets.hs)

head(kegg.sets.hs, 2)

foldchanges = res$log2FoldChange
names(foldchanges) = res$entrez
head(foldchanges)

pathview(col_path_df$AUC, pathway.id="04120")
