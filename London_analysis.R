# Load required packages
library(data.table)
library(reshape2)
# Source custom functions
source("date_functions.R")

# Read in cleaned data files
period <- data.table(read.csv("period-CLEANED.csv"))
setkey(period, Site_P)
context <- data.table(read.csv("context-CLEANED.csv"))
setkey(context, SITE_C)
sample <- data.table(read.csv("sample-CLEANED.csv"))
setkey(sample, SITE_S)

# Link tables together; drop un-needed fields
context.period <- merge(period[,list(Site_P,Start,End,MID)], context[,list(SITE_C,Site_P)], "Site_P", all=FALSE)
sample.period <- merge(context.period, sample[,list(SITE_S, SITE_C)], "SITE_C", all=FALSE)

# Plot histogram of context mid-points
# justifies concentration on last 2000 years
hist(context.period$MID, seq(-10000, 2000, 100))

# Calculate distribution of CONTEXTS across time
# First remove unecessary data to free up RAM
rm(context, period, sample)

# a) using aoristic sum method (should take 15-20 minutes with default arguments)
x <- aorist(context.period[, list(Start,End)])
barplot(x[,aoristic.sum], names.arg=breaks[1:(length(breaks)-1)])
write.csv(x, paste("contexts_aoristic_sum_by_period", params, ".csv",sep=""), row.names=FALSE)

# b) using simulation method (should take 3-5 minutes with default arguments)
set.seed(1901)
y <- date.simulate(context.period[, list(Start,End)], rep=500, bin.width=10)
write.csv(y, paste("contexts_simulated_by_period", params, ".csv",sep=""), row.names=FALSE)
boxplot(V1 ~ bin, data=y)
z <- y[,quantile(V1, probs=c(0.025,0.25,0.5,0.75,0.975)), by=bin]
z[,id:=c(0.025,0.25,0.5,0.75,0.975)]
z <- dcast.data.table(z, bin ~ id, value.var="V1")
write.csv(z, paste("summary_context_sim_by_period", params, ".csv",sep=""), row.names=FALSE)


# Calculate distribution of SAMPLES across time
rm(context.period)

# a) using aoristic sum method (should take 3-5 minutes with default arguments)
x <- aorist(sample.period[, list(Start,End)])
barplot(x[,aoristic.sum], names.arg=breaks[1:(length(breaks)-1)])
write.csv(x, paste("samples_aoristic_sum_by_period", params, ".csv",sep=""), row.names=FALSE)

# b) using simulation method (should take 10-20 seconds with default arguments)
set.seed(1982)
y <- date.simulate(sample.period[, list(Start,End)], rep=500, bin.width=10)
write.csv(y, paste("samples_simulated_by_period", params, ".csv",sep=""), row.names=FALSE)
boxplot(V1 ~ bin, data=y)
z <- y[,quantile(V1, probs=c(0.025,0.25,0.5,0.75,0.975)), by=bin]
z[,id:=c(0.025,0.25,0.5,0.75,0.975)]
z <- dcast.data.table(z, bin ~ id, value.var="V1")
write.csv(z, paste("summary_sample_sim_by_period", params, ".csv",sep=""), row.names=FALSE)


