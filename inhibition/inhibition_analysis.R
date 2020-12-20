#import mutation info
ccl_mut <- read.csv("C:/Users/james/Desktop/Sumoylation_Analysis/data/CCLE_mutations.csv") 
cclID <- read.csv("C:/Users/james/Desktop/Sumoylation_Analysis/data/sample_info.csv") #list of all cell lines, ID and lineage
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

#extract colorectal and pancreatic cancer only
mut_df <- mut_df[which(mut_df$lineage == "colorectal" | mut_df$lineage == "pancreas"), ]


#now process inhibition data by eliminating all spaces and special characters
inhibition$cell <- gsub(" ", "", inhibition$cell)
inhibition$cell <- gsub("-", "", inhibition$cell)
inhibition$cell <- gsub("\\.", "", inhibition$cell) #to remove dot, add 2 backslashes

colnames(inhibition)[which(names(inhibition) == "cell")] <- "ccl"

final_df <- merge(x = mut_df,
                  y = inhibition,
                  by = "ccl")
final_df <- unique(final_df)

colo_in <- final_df[final_df$lineage == "colorectal", ]
panc_in <- final_df[final_df$lineage == "pancreas", ]



#list top 5 common mutations
sort(table(final_df$Hugo_Symbol),decreasing=TRUE)[1:5]

#divide them into two groups: colorectal and pancreatic
sort(table(colo_in$Hugo_Symbol),decreasing=TRUE)[1:5]
sort(table(panc_in$Hugo_Symbol),decreasing=TRUE)[1:5]



#####################################  Pancreatic Cancer ##################################### 

#common mutation with below average IC50
below_avg_panc_IC50 <- panc_in[panc_in$Abs.IC50 < mean(panc_in$Abs.IC50), ]

panc_top5_mut_IC50 <- sort(table(below_avg_panc_IC50$Hugo_Symbol),decreasing=TRUE)[1:5]
panc_top5_mut_IC50
# TP53 KRAS TTN KRT17 MT-ND5 

#common mutation with below average AUC value
below_avg_panc_AUC <- panc_in[panc_in$AUC < mean(panc_in$AUC), ]
panc_top5_mut_AUC <- sort(table(below_avg_panc_AUC$Hugo_Symbol),decreasing=TRUE)[1:5]
panc_top5_mut_AUC
#same as the IC50 result

# statistical significance of difference in IC50 values between cells with that mutation and cells without the mutation
?t.test
TP53 <- below_avg_panc_IC50[below_avg_panc_IC50$Hugo_Symbol=="TP53", ]
NOT_TP53 <- below_avg_panc_IC50[below_avg_panc_IC50$Hugo_Symbol!="TP53", ]

t.test(TP53$Abs.IC50, )





##############################################################################################