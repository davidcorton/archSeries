# Load required packages
library(data.table)
library(reshape2)
# Source custom functions
source("date_functions.R")

# Read in data
period <- data.table(read.csv("period-CLEANED.csv"))
setkey(period, Site_P)
context <- data.table(read.csv("context-CLEANED.csv"))
setkey(context, SITE_C)
fish <- data.table(read.csv("Fish/Fish.csv"))
calib <- data.table(read.csv("aorist_samples_by_50.csv"))

# Merge data and remove un-needed data
context.period <- merge(period[,list(Site_P,Start,End,MID)], context[,list(SITE_C,Site_P)], "Site_P", all=FALSE)
fish.period <- merge(fish, context.period, by="SITE_C", all=FALSE)
rm(context, period, fish, context.period)

# Simulate
fresh <- freq.simulate(fish.period, calib, filter.field="Fresh_Marine", filter.values=c("Fresh", "Fresh/Marine"), quant.list=c(0.025,0.05,0.1,0.5,0.9,0.95,0.975), rep=2000)
marine <- freq.simulate(fish.period, calib, filter.field="Fresh_Marine", filter.values="Marine", quant.list=c(0.025,0.05,0.1,0.5,0.9,0.95,0.975), rep=2000)

# Plot frequencies
names <- fresh[[2]][,bin]
ticks <- seq(1, length(names),by=2)
x <- c(1:length(names), length(names):1)

boxplot(real ~ bin.no, data=fresh[[1]], xaxt="n", outline=FALSE, xlab="", ylab="Frequency")
axis(1, at=ticks, labels=names[ticks], las=2)
y <- c(fresh[[2]]$dummy_0.9, rev(fresh[[2]]$dummy_0.1))
polygon(x,y,col="firebrick")
boxplot(real ~ bin.no, data=fresh[[1]], main="Freshwater & migratory: frequency", xaxt="n", add=TRUE, outline=FALSE, col="darkolivegreen", xlab="", ylab="Frequency")

boxplot(real ~ bin.no, data=marine[[1]], xaxt="n", outline=FALSE, xlab="", ylab="Frequency")
axis(1, at=ticks, labels=names[ticks], las=2)
y <- c(marine[[2]]$dummy_0.9, rev(marine[[2]]$dummy_0.1))
polygon(x,y,col="firebrick")
boxplot(real ~ bin.no, data=marine[[1]], main="Marine: frequency", xaxt="n", add=TRUE, outline=FALSE, col="deepskyblue3", xlab="", ylab="Frequency")

# Plot rates of change
boxplot(ROC.real ~ bin.no, data=fresh[[1]], xaxt="n", outline=FALSE, xlab="", ylab="Rate of change in frequency")
axis(1, at=ticks, labels=names[ticks], las=2)
x <- c(1:(length(names)-1), (length(names)-1):1)
y <- c(fresh[[2]]$ROC.dummy_0.9, rev(fresh[[2]]$ROC.dummy_0.1))
y <- c(y[1:(length(y)/2-1)], y[(length(y)/2+2):length(y)])
polygon(x,y,col="firebrick")
boxplot(ROC.real ~ bin.no, data=fresh[[1]], main="Freshwater & migratory: rate of change", xaxt="n", add=TRUE, col="darkolivegreen", outline=FALSE, xlab="", ylab="Rate of change in frequency")

#boxplot(ROC.dummy ~ bin.no, data=marine[[1]], main="Marine", names=fresh[[2]][!bin=="1950-2000",bin], outline=FALSE, xlab="Time (50yr intervals)", ylab="Rate of change in frequency")
boxplot(ROC.real ~ bin.no, data=marine[[1]], xaxt="n", outline=FALSE, xlab="", ylab="Rate of change in frequency")
axis(1, at=ticks, labels=names[ticks], las=2)
x <- c(1:(length(names)-1), (length(names)-1):1)
y <- c(marine[[2]]$ROC.dummy_0.9, rev(marine[[2]]$ROC.dummy_0.1))
y <- c(y[1:(length(y)/2-1)], y[(length(y)/2+2):length(y)])
polygon(x,y,col="firebrick")
boxplot(ROC.real ~ bin.no, data=marine[[1]], main="Marine: rate of change", xaxt="n", add=TRUE, col="deepskyblue3", outline=FALSE, xlab="", ylab="Rate of change in frequency")

# Calculate bones/litre
x <- 1:length(names)
fresh.upper <- fresh[[2]]$real_0.9/calib$V1
fresh.mid <- fresh[[2]]$real_0.5/calib$V1
fresh.lower <- fresh[[2]]$real_0.1/calib$V1
x.poly <- c(1:length(names), length(names):1)
y.poly <- c(fresh.upper, rev(fresh.lower))
plot(x, fresh.upper, type="n", main="Freshwater & migratory: catch per litre", xaxt="n", ylab="Frags per litre processed", xlab="")
ticks <- seq(1, length(names),by=2)
axis(1, at=ticks, labels=names[ticks], las=2)
polygon(x.poly, y.poly, col="darkolivegreen3")
lines(x, fresh.mid, col="darkgreen")

marine.upper <- marine[[2]]$real_0.9/calib$V1
marine.mid <- marine[[2]]$real_0.5/calib$V1
marine.lower <- marine[[2]]$real_0.1/calib$V1
x.poly <- c(1:length(names), length(names):1)
y.poly <- c(marine.upper, rev(marine.lower))
plot(x, marine.upper, type="n", main="Marine: catch per litre", xaxt="n", ylab="Frags per litre processed", xlab="")
ticks <- seq(1, length(names),by=2)
axis(1, at=ticks, labels=names[ticks], las=2)
polygon(x.poly, y.poly, col="deepskyblue2")
lines(x, marine.mid, col="dodgerblue4")
