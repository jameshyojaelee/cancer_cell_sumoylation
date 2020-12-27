library(dplyr)


#import mutation info
ccl_mut <- read.csv("C:/Users/james/Desktop/Sumoylation_Analysis/data/CCLE_mutations.csv") 
cclID <- read.csv("C:/Users/james/Desktop/Sumoylation_Analysis/data/sample_info.csv") #list of all cell lines, ID and lineage
inhibition <- read.csv("inhibition_data.csv", fileEncoding="UTF-8-BOM") #inhibition data 

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

nrow(distinct(mut_df, ccl, .keep_all = TRUE))
#total of 1694 unique cell lines

#extract colorectal and pancreatic cancer only
#mut_df <- mut_df[which(mut_df$lineage == "colorectal" | mut_df$lineage == "pancreas"), ]
#nrow(distinct(mut_df, ccl, .keep_all = TRUE))
#total of 130 unique colorectal/pancreatic cell lines

#now process inhibition data by eliminating all spaces and special characters
inhibition$cell <- gsub(" ", "", inhibition$cell)
inhibition$cell <- gsub("-", "", inhibition$cell)
inhibition$cell <- gsub("\\.", "", inhibition$cell) #to remove dot, add 2 backslashes

inhibition$cell <- toupper(inhibition$cell)


colnames(inhibition)[which(names(inhibition) == "cell")] <- "ccl"


final_df <- merge(x = mut_df,
                  y = inhibition,
                  by = "ccl")

final_df <- unique(final_df)
distinct(final_df, ccl, .keep_all = TRUE)
nrow(distinct(final_df, ccl, .keep_all = TRUE))
#30 unique cell lines
#DLD1 is the only cell line that doesn't exist in the BROAD institute's dataset
#COLO741 is a colon carcinoma and therefore categorized as "skin" cancer. 

#all the unique cell lines
distinct(final_df, lineage)

colo_in <- final_df[final_df$lineage == "colorectal", ]
nrow(distinct(colo_in, ccl))

#need to include skin cancer because of COLO741 being categorized as colon carcinoma
panc_in <- final_df[final_df$lineage == "pancreas" | final_df$lineage == "skin", ]
nrow(distinct(panc_in, ccl))



#list top 5 common mutations
sort(table(final_df$Hugo_Symbol),decreasing=TRUE)[1:5]

#divide them into two groups: colorectal and pancreatic
sort(table(colo_in$Hugo_Symbol),decreasing=TRUE)[1:10]
sort(table(panc_in$Hugo_Symbol),decreasing=TRUE)[1:10]



##################################################  Pancreatic Cancer ################################################## 


#common mutation with below average AUC value
below_avg_panc_AUC <- panc_in[panc_in$AUC < 400, ]
nrow(distinct(below_avg_panc_AUC, ccl))
#list the 10 most common mutations
panc_below_avg <- sort(table(below_avg_panc_AUC$Hugo_Symbol),decreasing=TRUE)[1:10]
panc_below_avg



#common mutation with above average AUC value
above_avg_panc_AUC <- panc_in[panc_in$AUC >= 400, ]
nrow(distinct(above_avg_panc_AUC, ccl))
#list the 10 most common mutations
panc_above_avg <- sort(table(above_avg_panc_AUC$Hugo_Symbol),decreasing=TRUE)[1:10]
panc_above_avg




# statistical significance of difference in AUC values between cells with that mutation and cells without the mutation

ttest_df <- panc_in
ttest_df$Hugo_Symbol <- as.character(ttest_df$Hugo_Symbol)
ttest_df$Hugo_Symbol[ttest_df$Hugo_Symbol != "TP53"] <- "Lacks Mut"
ttest_df$Hugo_Symbol[ttest_df$Hugo_Symbol == "TP53"] <- "Has Mut"
distinct(ttest_df, Hugo_Symbol)

t <- list()
t[[1]] <- t.test(AUC ~ Hugo_Symbol, data = ttest_df)

tt <- sapply(t, function(x) {
          c(x$estimate[1],
            x$estimate[2],
            ci.lower = x$conf.int[1],
            ci.upper = x$conf.int[2],
            p.value = x$p.value)
        }
      )

#create a function to repeat t-test for other genes

ttest <- function(G){
  ttest_df <- panc_in
  ttest_df$Hugo_Symbol <- as.character(ttest_df$Hugo_Symbol)
  ttest_df$Hugo_Symbol[ttest_df$Hugo_Symbol != G] <- "Lacks Mut"
  ttest_df$Hugo_Symbol[ttest_df$Hugo_Symbol == G] <- "Has Mut"
  t.test(AUC ~ Hugo_Symbol, data = ttest_df)
}

ttest("KRAS")


#create dataframe with top 10 mutations (with below 400 AUC)
panc_below_avg_df <- as.data.frame(panc_below_avg, stringsAsFactors=FALSE)
colnames(panc_below_avg_df) <- c('mutation', 'freq')

panc_below_avg_df$mutated_mean <- NA
panc_below_avg_df$others_mean <- NA
panc_below_avg_df$p.value <- NA


for (i in nrow(panc_below_avg_df)) {
  t <- ttest(panc_below_avg_df$mutation[i])
  panc_below_avg_df$mutated_mean <- 
}



########################################################################################################################





##################################################  Colorectal Cancer ################################################## 

#common mutation with below average IC50
below_avg_colo_IC50 <- colo_in[colo_in$Abs.IC50 < mean(colo_in$Abs.IC50), ]

colo_top5_mut_IC50 <- sort(table(below_avg_colo_IC50$Hugo_Symbol),decreasing=TRUE)[1:10]
colo_top5_mut_IC50
# TP53 KRAS TTN KRT17 MT-ND5  

#number of cell lines with below average IC50
nrow(distinct(below_avg_colo_IC50, ccl, .keep_all = TRUE))

#common mutation with below average AUC value
below_avg_colo_AUC <- colo_in[colo_in$AUC < mean(colo_in$AUC), ]
colo_top5_mut_AUC <- sort(table(below_avg_colo_AUC$Hugo_Symbol),decreasing=TRUE)[1:10]
colo_top5_mut_AUC
#same as the IC50 result

#number of  cell lines with below average AUC
nrow(distinct(below_avg_colo_AUC, ccl, .keep_all = TRUE))

########################################################################################################################