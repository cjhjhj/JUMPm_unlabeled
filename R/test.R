rm(list = ls())

# args = commandArgs(trailingOnly = TRUE)
# cat(args[1])
# cat(args[2])

file = "C:/Research/Projects/7Metabolomics/JUMPm/IROAsamples/IROA_IS_NEG_1.1.feature"
df = read.table(file, header = T, sep = "\t", row.names = NULL, stringsAsFactors = F, comment.char = "", check.names = F)

for (i in 1:ncol(df)) {
    if (colnames(df)[i] == "minMS1ScanNumber") {
        colnames(df)[i] = "minMS1"
    } else if (colnames(df)[i] == "maxMS1ScanNumber") {
        colnames(df)[i] = "maxMS1"
    } else if (colnames(df)[i] == "Intensity") {
        colnames(df)[i] = "intensity"
    } else if (colnames(df)[i] == "Percentage of TF") {
        colnames(df)[i] = "PercentageTF"
    }
}