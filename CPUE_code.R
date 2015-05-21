# Load required packages
library(data.table)
# Source custom functions
source("date_functions.R")

#Read in data
period <- data.table(read.csv("period-CLEANED.csv"))
setkey(period, Site_P)
context <- data.table(read.csv("context-CLEANED.csv"))
setkey(context, SITE_C)
sample <- data.table(read.csv("sample-CLEANED.csv"))
setkey(sample, SITE_S)
fish <- data.table(read.csv("Fish/Fish.csv"))

# Merge data and remove un-needed data
context.period <- merge(period[,list(Site_P,Start,End,MID)], context[,list(SITE_C,Site_P)], "Site_P", all=FALSE)
sample.period <- merge(context.period, sample[,list(SITE_S, SITE_C)], "SITE_C", all=FALSE)
fish.period <- merge(fish, context.period, by="SITE_C", all=FALSE)
rm(context, period, fish, sample)

# 1. Illustrate fish frequencies against dummy set
set.seed(1901)
fish <- freq.simulate(sample.period[, list(Start,End)], rep=2000, bin.width=100, RoC=TRUE, summ=TRUE)

#set up axes
with(fish[[1]], plot(bin.no, count, xlab="", xaxt="n", ylab="Estimated frequency density", type="n"))
names <- unique(fish[[1]]$bin)
ticks <- seq(1, length(names),by=1)
axis(1, at=ticks, labels=unique(fish[[1]]$bin)[ticks], las=2)

#plot boxes
with(fish[[1]], boxplot(count~bin.no, outline=FALSE, xaxt="n", add=TRUE))

#plot polygon
x <- c(1:length(names), length(names):1)
y <- c(fish[[2]]$dummy_0.975, rev(fish[[2]]$dummy_0.025))
polygon(x,y,col="firebrick")

#plot median line
with(fish[[2]], lines(count_0.5, col="white"))

#plot all lines
for(i in 1:2000) {
    with(fish[[1]], lines(count[rep.no==i], type="l", col=rgb(46,25,235,20, maxColorValue=255)))
}





# 2. Illustrate rate of change


# Simulate
fresh <- freq.simulate(fish.period, ROC=FALSE, calib, filter.field="Fresh_Marine", filter.values=c("Fresh", "Fresh/Marine"), quant.list=c(0.025,0.05,0.1,0.5,0.9,0.95,0.975), rep=2000)
marine <- freq.simulate(fish.period, ROC=FALSE, calib, filter.field="Fresh_Marine", filter.values="Marine", quant.list=c(0.025,0.05,0.1,0.5,0.9,0.95,0.975), rep=2000)

