# Load required packages
library(data.table)
# Source custom functions
source("date_functions.R")
source("plot_functions.R")

#Read in data
period <- data.table(read.csv("period-CLEANED.csv"))
setkey(period, Site_P)
context <- data.table(read.csv("context-CLEANED.csv"))
setkey(context, SITE_C)
sample <- data.table(read.csv("sample-CLEANED.csv"))
setkey(sample, SITE_S)
zoo <- data.table(read.csv("maintable.csv"))
setkey(zoo, SITE_C)
species <- data.table(read.csv("species_codes.csv"))

# Merge data and remove un-needed data
context.period <- merge(period[,list(Site_P,Start,End,MID)], context[,list(SITE_C,Site_P)], "Site_P", all=FALSE)
sample.period <- merge(context.period, sample[,list(SITE_S, SITE_C, PROC_VOL)], "SITE_C", all=FALSE)
zoo <- merge(zoo, species, by="SPECIES", all.x=TRUE, all.y=FALSE)
zoo[,SITE_S:=paste(SITECODE, SAMPLE, sep="-")]
zoo.period <- merge(zoo, context.period, by="SITE_C", all=FALSE)
zoo.sample <- merge(zoo, sample.period, by="SITE_S", all=FALSE)
rm(context, period, sample, zoo, species)

#1. Plot distribution of fish bones
all.fish.frags <- freq.simulate(zoo.period[CLASS=="fish"], weight=zoo.period[CLASS=="fish"]$FRAG_COUNT, reps=2000)
box.chron(all.fish.frags, field.list="count", main="Abundance of fish remains by century")

axis.setup(all.fish.frags, main="Abundance of fish remains by century")
poly.chron(all.fish.frags, field.list="dummy", col.list="grey", add=TRUE)
box.chron(all.fish.frags, field.list="count", col.list="darkred", add=TRUE)

lines.chron(all.fish.frags, col.list=c("darkred", "grey"), main="Abundance of fish remains by century")

#1b try ubiquity approach
all.fish <- freq.simulate(zoo.period[CLASS=="fish"], reps=2000)
lines.chron(all.fish, col.list=c("darkred", "grey"), main="Contexts with fish remains by century")

#2. RoC

#3. Compare main mammals
cattle <- freq.simulate(zoo.period[SPECIES=="BOS"], weight=zoo.period[SPECIES=="BOS"]$FRAG_COUNT, reps=2000)
caprines <- freq.simulate(zoo.period[SPECIES%in%c("OVCA", "OVI", "CRA")], weight=zoo.period[SPECIES%in%c("OVCA", "OVI", "CRA")]$FRAG_COUNT, reps=2000)
pigs <- freq.simulate(zoo.period[SPECIES=="SUS"], weight=zoo.period[SPECIES=="SUS"]$FRAG_COUNT, reps=2000)
chicken <- freq.simulate(zoo.period[SPECIES=="CHIK"], weight=zoo.period[SPECIES=="CHIK"]$FRAG_COUNT, reps=2000)
poly.chron(cattle, field.list="count", main="Abundance of main domesticates", legend=FALSE)
poly.chron(caprines, field.list="count", col.list="darkgreen", add=TRUE, legend=FALSE)
poly.chron(pigs, field.list="count", col.list="darkblue", add=TRUE, legend=FALSE)
poly.chron(chicken, field.list="count", col.list="grey30", add=TRUE, legend=FALSE)
legend("topright", legend=c("Cattle", "Caprines", "Pigs", "Chicken"), fill=c("darkred", "darkgreen", "darkblue", "grey30"), bty="n")

#4. Calculate distributions of context and samples
context.sim <- date.simulate(context.period, context.fields="SITE_C", reps=500)
#context <- aorist(context.period)
#write.csv(context, "CONTEXT-AORIST.csv", row.names=FALSE)
context <- data.table(read.csv("CONTEXT-AORIST.csv"))

sample.sim <- date.simulate(sample.period, reps=500)
sample.count <- aorist(sample.period)
write.csv(sample.count, "SAMPLE-COUNT-AORIST.csv", row.names=FALSE)
sample.vol <- aorist(sample.period, weight=sample.period$PROC_VOL)
write.csv(sample.vol, "SAMPLE-VOL-AORIST.csv")
barplot(context$aorist)
barplot(sample.count$aorist)

#5. Use these to create better dummy sets
doms <- freq.simulate(zoo.period[SPECIES%in%c("BOS", "OVCA", "OVI", "CRA", "SUS", "CHIK")], weight=zoo.period[SPECIES%in%c("BOS", "OVCA", "OVI", "CRA", "SUS", "CHIK")]$FRAG_COUNT, probs=context$aorist, reps=1000)
poly.chron(doms, field.list=c("dummy", "count"), col.list=c("darkgreen", "darkred"), main="Abundance of major domesticates")
lines.chron(doms, main="Abundance of major domesticates")

all.fish.frags <- freq.simulate(zoo.period[CLASS=="fish"], weight=zoo.period[CLASS=="fish"]$FRAG_COUNT, probs=sample.count$aorist, reps=2000)
poly.chron(all.fish.frags, field.list=c("dummy", "count"), col.list=c("darkgreen", "darkred"), main="Abundance of fish in samples")
lines.chron(all.fish.frags, main="Abundance of fish in samples")

#5b or use the mammals as a dummy for the fish!
all.mammal.frags <- date.simulate(zoo.period[CLASS=="mammal"], weight=zoo.period[CLASS=="mammal"]$FRAG_COUNT, reps=1000)
poly.chron(all.mammal.frags)

mam.med <- all.mammal.frags[[2]][quantile==0.5, V1]
fish.vs.mammals <- freq.simulate(zoo.period[CLASS=="fish"], weight=zoo.period[CLASS=="fish"]$FRAG_COUNT, probs=mam.med, reps=1000)
poly.chron(fish.vs.mammals, field.list=c("dummy", "count"), col.list=c("darkgreen", "darkred"), main="Fish, with dummy set based on mammals")

#6 catch per unit
fish <- cbind(all.fish.frags[[2]], all.mammal.frags[[2]]$V1)
fish[,per.mam:=V1/V2]
poly.chron(fish, field.list="count", col.list="darkblue", value.field="per.mam", main="Fish per mammal bone", legend=FALSE)

####
Compare context types

#####
context <- aorist(context.period, bin.width=50)
write.csv(context, "CONTEXT-AORIST-50.csv", row.names=FALSE)




