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
context.period <- merge(period, context[,list(CONTEXT,SITE_C,Site_P)], "Site_P", all=FALSE)

# Plot histogram of context mid-points
# justifies concentration on last 2000 years
hist(context.period$MID, seq(-10000, 2000, 100))

# Set up break points
start.date = 0
end.date = 2000
bin.width = 100
breaks <- seq(start.date, end.date, bin.width)

#Define function to calculate aoristic sum
#Arguments: 'data' is a data table with two numeric columns called Start and End
# 'start.date' is a single numeric value for the start of the time period to be analysed
# 'end.date' is a single numric value for the end end of the time period to be analysed
# 'bin.width' is a single numeric value setting the resolution of the analysis, in years
# 'weight' (OPTIONAL) is a numeric vector giving a weight for each context/entity
aorist <- function(data, start.date=0, end.date=2000, bin.width=100, weight=1) { 
    require(data.table)
    aoristic.sum <- numeric(((end.date-start.date)/bin.width))
    data <- cbind(data, weight)
    data <- data[End >= start.date & Start <= end.date]
    data[,a:=ceiling((Start-start.date)/bin.width)+1]
    data[,b:=floor((End-start.date)/bin.width)]
    data[,lead.in:=(ceiling((Start-start.date)/bin.width)-((Start-start.date)/bin.width))*bin.width]
    data[,lead.out:=(((End-start.date)/bin.width)-floor((End-start.date)/bin.width))*bin.width]
    data[,diff:=b-a]
    data[,duration:={ifelse(diff==-2, lead.in+lead.out, End-Start)}]
    data[,full.prob:=(bin.width/duration)*weight]
    data[,in.prob:=(lead.in/duration)*weight]
    data[,out.prob:=(lead.out/duration)*weight]
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

date.simulate <- function(data, start.date=0, end.date=2000, bin.width=100, rep=100, weight=1, save.raw=NULL)) {
    require(data.table)
    data <- cbind(data, weight)
    breaks <- seq(start.date, end.date, bin.width)
    labels <- 1:(length(breaks)-1)
    data.sim.cut <- function(data) {
        x <- runif(nrow(data))
        y <- (x*(data[,End]-data[,Start]))+data[,Start]
        cut(y, breaks, labels=labels)
    }
   # z <- replicate(rep, data.sim(data))
    data[,as.character(1:rep):=replicate(rep, data.sim(data), simplify=FALSE)]
    #if(is.null(save.raw) == FALSE) (write.csv(data, save.raw))
    #data[,as.character(1:rep):=replicate()]
}

#Testing
data <- context.period[1:1000, list(Start, End)]
data[8:10, Start:=0]
data[8:10, End:=500]

x <- aorist(data, start.date, end.date, bin.width)
system.time(aorist(data, start.date, end.date, bin.width))


cut(data$rep2, breaks, labels=labels)

x <- list("rep1", "rep2")
rep=100
data[,as.character(1:rep):=list

    
data.sim <- function(data) {
    x <- runif(nrow(data))
    (x*(data[,End]-data[,Start]))+data[,Start]
}
rep=1000
y <- data.sim(data)

z <- lapply(rep, data.sim, data=data)



