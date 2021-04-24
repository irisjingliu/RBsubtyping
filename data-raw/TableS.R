x = read.csv("/Users/jing/Google Drive/Data/TableS1_FinalSubtype.csv", sep=";", stringsAsFactors = F)
rownames(x) = x$sample_ID
x$sample_ID
TableS1_FinalSubtype = x
# save(TableS1_FinalSubtype, file = "/Users/jing/Google Drive/Data/TableS1_FinalSubtype.RData")

setwd("/Users/jing/Downloads/myfirstpackage/RBsubtyping")
usethis::use_data(TableS1_FinalSubtype)



x = read.csv("/Users/jing/Desktop/myRpackage/RBsubtyping/data-raw/TableS2_dataFig2D.csv", sep=";", stringsAsFactors = F)
rownames(x) = x$sample_ID
x$sample_ID
TableS2_DataFig2D = x
# save(TableS1_FinalSubtype, file = "/Users/jing/Google Drive/Data/TableS1_FinalSubtype.RData")

setwd("/Users/jing/Desktop/myRpackage/RBsubtyping")
usethis::use_data(TableS2_DataFig2D)
