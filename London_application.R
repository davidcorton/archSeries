# APPLICATION TO COD DATA
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

# Free up RAM
rm(context, context.period, period, sample)

# read in species data
cod <- data.table(read.csv("cod/cod_data_for_simulation.csv")[1:3])
# set resolution
bin.width <- 50


# simulate from actual data
set.seed(1483)
gm <- date.simulate(cod[,list(Start,End)], weight=cod[,weight], bin.width=bin.width, rep=1000)
write.csv(gm, paste("cod/cod_simulated_by_period", params, ".csv",sep=""), row.names=FALSE)
boxplot(V1 ~ bin, data=gm)
gms <- gm[,quantile(V1, probs=c(0.025,0.25,0.5,0.75,0.975)), by=bin]
gms[,id:=c(0.025,0.25,0.5,0.75,0.975)]
gms <- dcast.data.table(gms, bin ~ id, value.var="V1")
write.csv(gms, paste("cod/summary_cod_sim_by_period", params, ".csv",sep=""), row.names=FALSE)

# create calibration distribution
x <- aorist(sample.period[, list(Start,End)], bin.width=bin.width)

# simulate dummy set from calibration distribution
set.seed(593)
gm.dummy <- dummy.simulate(probs=x, weight=cod[,weight], rep=1000)
boxplot(V1 ~ bin, data=gm.dummy)
gms.dummy <- gm.dummy[,quantile(V1, probs=c(0.025,0.25,0.5,0.75,0.975)), by=bin]
gms.dummy[,id:=c(0.025,0.25,0.5,0.75,0.975)]
gms.dummy <- dcast.data.table(gms.dummy, bin ~ id, value.var="V1")
write.csv(gms.dummy, paste("cod/summary_cod_dummy_by_period", params, ".csv",sep=""), row.names=FALSE)