#Pathway Analysis

library(BiocManager)
library(DESeq2)
library(pathview)
library(gage)
library(gageData)
library("AnnotationDbi")
library("org.Hs.eg.db")

columns(org.Hs.eg.db)


#import inhibition data with IC50 and AUC values
indf <- read.csv("C:/Users/james/Desktop/Gene_knockout/project/inhibition/inhibition_data.csv", fileEncoding="UTF-8-BOM")

head(indf)



inhibition <- indf
inhibition[] <- lapply(col_AUC, as.character)

class(inhibition$cell)

head(inhibition)

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
