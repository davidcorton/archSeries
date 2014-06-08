# Script for extraction and cleaning of London sample data
# Requires following files to be present in working directory:
# 1. Period.csv
# 2. Context.csv
# 3. Samples.csv
# 4. Finds_date.csv
# Original Excel versions need to be saved as CSV files with the option set to
# quote all text variables; otherwise commas in comment fields will be treated as field
# breaks
# nb. reading in large xlsx files is extremely slow, xls is much faster but limited to
# 65535 lines. Hence necessary to save as CSV in advance for larger files.

#Load needed packages
library("data.table")

# Read in data and select relevant columns
period <- read.csv("Period.csv")
period <- data.table(period[,c(1,2,3,5,6,7)])
setkey(period, Site_P)
finds <- read.csv("Finds_date.csv")
finds <- data.table(finds[,c(3,4,6,7,8,9)])
setkey(finds, SITE_C)
context <- read.csv("Context.csv")
context <- data.table(context[,c(1,2,4,15,16)])
setkey(context, SITE_C)
sample <- read.csv("Samples.csv")
sample <- data.table(sample[,c(2,3,4,7,9,11,13,16,17,18,19,20)])
setkey(sample, SITE_S)

# Create columns in period table for range and mid-points
period[,RANGE:=End-Start]
period[,MID:=Start+(RANGE/2)]

# Remove rows with missing data, and save into CSV files with 'REJECTED' in names
# Save remaining rows into (a) CSV files called 'CLEANED' and (b) list called data.cleaned
sampleshort <- sample[,c(1,2,3,5,8), with=FALSE]
table.list <- list(period=period, context=context, finds=finds, sample=sampleshort)
good <- as.matrix(sapply(table.list, function(x) complete.cases(x)))
table.list[4] <- list(sample)
clean <- function(x) {
    data <- table.list[[x]]
    if(sum(!good[[x]]) > 0) {
        reject.name <- paste(names(table.list[x]), "-REJECTED", ".csv", sep="")
        write.csv(data[!good[[x]]], reject.name)
    }
    output.name <- paste(names(table.list[x]), "-CLEANED", ".csv", sep="")
    write.csv(data[good[[x]]], output.name)
    data[good[[x]]]
}
data.cleaned <- lapply(1:4, clean)

# Link contexts to periods; drop 
