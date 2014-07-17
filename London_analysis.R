# Load required packages
library(data.table)

# Read in cleaned data files
period <- data.table(read.csv("period-CLEANED.csv"))
setkey(period, Site_P)
finds <- data.table(read.csv("finds-CLEANED.csv"))
setkey(finds, SITE_C)
context <- data.table(read.csv("context-CLEANED.csv"))
setkey(context, SITE_C)
sample <- data.table(read.csv("sample-CLEANED.csv"))
setkey(sample, SITE_S)

# Link contexts to periods; drop duplicate variables
context.period <- merge(data.cleaned[[1]], data.cleaned[[2]][,list(CONTEXT,SITE_C,Site_P)], "Site_P", all=FALSE)

# Plot histogram of context mid-points
# justifies concentration on last 4000 years
hist(context.period$MID, seq(-10000, 2000, 100))

# Set up break points
start.date = 0
end.date = 2000
bin.width = 100
breaks <- seq(start.date, end.date, bin.width)

#Define function to calculate aoristic sum
#Arguments: 'data' is a 2-column data table (Start, End)
# 'start.date', 'end.date', and 'bin.width' set up the bins to use
aorist <- function(data, start.date, end.date, bin.width) { 
    aoristic.sum <- numeric(((end.date-start.date)/bin.width))
    data <- data[End >= start.date & Start <= end.date]
    data[,a:=ceiling((Start-start.date)/bin.width)+1]
    data[,b:=floor((End-start.date)/bin.width)]
    data[,lead.in:=(ceiling((Start-start.date)/bin.width)-((Start-start.date)/bin.width))*bin.width]
    data[,lead.out:=(((End-start.date)/bin.width)-floor((End-start.date)/bin.width))*bin.width]
    data[,diff:=b-a]
    data[,duration:={ifelse(diff==-2, lead.in+lead.out, End-Start)}]
    data[,full.prob:=bin.width/duration]
    data[,in.prob:=lead.in/duration]
    data[,out.prob:=lead.out/duration]
    for(i in 1:nrow(data)) {
        aoristic.sum[data[i,a]-1] <- aoristic.sum[data[i,a]-1] + data[i,in.prob]
        aoristic.sum[data[i,b]+1] <- aoristic.sum[data[i,b]+1] + data[i,out.prob]
        if(data[i,diff] >= 0) {
            for(j in data[i,a]:data[i,b]) {
                aoristic.sum[j] <- aoristic.sum[j] + data[i,full.prob]
            }
        }
        if(i/1000 == round(i/1000)) {print(paste(i/nrow(data)*100, "percent complete"))}
    }
    aoristic.sum
}

#Testing
data <- context.period[1:2, list(Start, End)]
x <- aorist(data, start.date, end.date, bin.width)
system.time(aorist(data, start.date, end.date, bin.width))