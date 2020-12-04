#PCA Analysis of UBA2

head(UBA2_merged)

UBA2_exp <- data.frame(ID = exp$X, UBA2_exp = exp$UBA2..10054.)


UBA2_exp <- merge(x = UBA2_merged, 
                  y = UBA2_exp,
                  by = "ID")

head(UBA2_exp)
