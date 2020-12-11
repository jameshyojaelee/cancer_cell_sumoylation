#import mutation info
ccl_mut <- read.csv("C:/Users/james/Desktop/Gene_knockout/data/CCLE_mutations.csv") 
cclID <- read.csv("C:/Users/james/Desktop/Gene_knockout/data/sample_info.csv") #list of all cell lines, ID and lineage
inhibition <- read.csv("inhibition/inhibition_data.csv", fileEncoding="UTF-8-BOM") #inhibition data 

#dataframe with the cancer cell line names
ccl_df <- data.frame(ID = cclID$DepMap_ID, ccl = cclID$stripped_cell_line_name, lineage = cclID$lineage)


#change the column name of DepMap_ID to ID
colnames(ccl_mut)[which(names(ccl_mut) == "DepMap_ID")] <- "ID"
colnames(ccl_df)[which(names(ccl_df) == "DepMap_ID")] <- "ID"

#subset 3 columns only
mut <- subset(ccl_mut, select = c("ID", "Hugo_Symbol", "Variant_Classification"))

#import the auc data to manipulate
mut_df <- merge(x = ccl_df, 
                y = mut,
                by = "ID")

#eliminate silent mutation
mut_df <- mut_df[mut_df$Variant_Classification != "Silent", ]
mut_df <- mut_df[which(mut_df$lineage == "colorectal" | mut_df$lineage == "pancreas"), ]


#now process inhibition data by eliminating all spaces and special characters
inhibition$cell <- gsub(" ", "", inhibition$cell)
inhibition$cell <- gsub("-", "", inhibition$cell)
inhibition$cell <- gsub("\\.", "", inhibition$cell) #to remove dot, add 2 backslashes

colnames(inhibition)[which(names(inhibition) == "cell")] <- "ccl"

final_df <- merge(x = mut_df,
                  y = inhibition,
                  by = "ccl")

colo_in <- final_df[final_df$lineage == "colorectal", ]
panc_in <- final_df[final_df$lineage == "pancreas", ]


#calculate the mode (most frequent mutation)
calculate_mode <- function(x) {
  uniqx <- unique(na.omit(x))
  uniqx[which.max(tabulate(match(x, uniqx)))]
}

#mode of final_df
calculate_mode(final_df$Hugo_Symbol)

xx <- sort(table(final_df$Hugo_Symbol),decreasing=TRUE)[1:5]
xx


#mode of panc_in
calculate_mode(panc_in$Hugo_Symbol)

xx <- sort(table(panc_in$Hugo_Symbol),decreasing=TRUE)[1:5]
xx

#
below_avg_panc_in <- panc_in[panc_in$Abs.IC50 < mean(panc_in$Abs.IC50), ]


#most frequent mutations
calculate_mode(below_avg_panc_in$Hugo_Symbol)

#or top 3 frequence mutations
xx <- sort(table(below_avg_panc_in$Hugo_Symbol),decreasing=TRUE)[1:5]
xx

nrow(below_avg_panc_in[below_avg_panc_in$Hugo_Symbol == "TP53", ])


below_avg_panc_in <- panc_in[panc_in$AUC < mean(panc_in$AUC), ]


