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

volumes <- data.table(read.csv("Sample volume-CLEANED.csv"))
volumes[,SITE_S:=paste(SITECODE, SAMPLE, sep="-")]

# Merge data and remove un-needed data
sample <- merge(sample[,list(SITE_S, SITE_C, PROC_VOL)], volumes[,list(SITE_S, WTS_VOL)], by="SITE_S", all=TRUE)
context.period <- merge(period[,list(Site_P,Start,End,MID)], context[,list(SITE_C,Site_P)], "Site_P", all=FALSE)
sample.period <- merge(context.period, sample, by="SITE_C", all=FALSE)
zoo <- merge(zoo, species, by="SPECIES", all.x=TRUE, all.y=FALSE)
zoo[,SITE_S:=paste(SITECODE, SAMPLE, sep="-")]
zoo.period <- merge(zoo, context.period, by="SITE_C", all=FALSE)
zoo.sample <- merge(zoo, sample.period, by="SITE_S", all=FALSE)
rm(context, period, sample, zoo, species)

# Sample sizes
sum(zoo.period[CLASS=="fish", FRAG_COUNT])
length(unique(zoo.period[CLASS=="fish", SITECODE]))
length(unique(zoo.period[CLASS=="fish", SITE_C]))
length(unique(zoo.period[CLASS=="fish", SITE_S]))

length(unique(zoo.period[, SITECODE]))
length(unique(zoo.period[, SITE_C]))
length(unique(zoo.period[, SITE_S]))


#1. Plot distribution of fish bones
fish.aorist <- aorist(zoo.period[CLASS=="fish"], weight=zoo.period[CLASS=="fish", FRAG_COUNT], bin.width=50)
with(fish.aorist, barplot(aorist, ylab="Estimated frequency density", main="Aoristic sum of London fish (50-yr bins)", col="darkred", names.arg=bin, las=2, cex.lab=1.2, cex.names=1.05, cex.axis=1.05, cex.main=1.2))

all.fish.frags <- freq.simulate(zoo.period[CLASS=="fish"], weight=zoo.period[CLASS=="fish"]$FRAG_COUNT, bin.width=50, reps=1000)
box.chron(all.fish.frags, field.list="count", main="Simulated distribution of London fish (50-yr bins)", lab.sp=2)

axis.setup(all.fish.frags, main="Simulated distribution of London fish (50-yr bins)", lab.sp=2)
poly.chron(all.fish.frags, field.list="dummy", col.list="grey", add=TRUE)
box.chron(all.fish.frags, field.list="count", col.list="darkred", add=TRUE)

lines.chron(all.fish.frags, col.list=c("darkred", "grey"), main="Simulated distribution of London fish (50-yr bins)", lab.sp=2)

#1b try ubiquity approach
all.fish <- freq.simulate(zoo.period[CLASS=="fish"], reps=2000)
lines.chron(all.fish, col.list=c("darkred", "grey"), main="Contexts with fish remains by century")

#2. RoC
all.fish.frags <- freq.simulate(zoo.period[CLASS=="fish"], weight=zoo.period[CLASS=="fish"]$FRAG_COUNT, bin.width=50, reps=1000, RoC=TRUE)
axis.setup(all.fish.frags, field.list=c("RoC.count", "RoC.dummy"), main="Simulated rates of change in fish frequency per 50yrs", ylab="Estimated rate of change (%)", lab.sp=1)
poly.chron(all.fish.frags, field.list="RoC.dummy", col.list="grey", add=TRUE)
box.chron(all.fish.frags, field.list="RoC.count", col.list="darkred", add=TRUE)


#3. Compare main mammals
cattle <- freq.simulate(zoo.period[SPECIES=="BOS"], weight=zoo.period[SPECIES=="BOS"]$FRAG_COUNT, reps=2000, bin.width=50)
caprines <- freq.simulate(zoo.period[SPECIES%in%c("OVCA", "OVI", "CRA")], weight=zoo.period[SPECIES%in%c("OVCA", "OVI", "CRA")]$FRAG_COUNT, reps=2000, bin.width=50)
pigs <- freq.simulate(zoo.period[SPECIES=="SUS"], weight=zoo.period[SPECIES=="SUS"]$FRAG_COUNT, reps=2000, bin.width=50)
chicken <- freq.simulate(zoo.period[SPECIES=="CHIK"], weight=zoo.period[SPECIES=="CHIK"]$FRAG_COUNT, reps=2000, bin.width=50)
poly.chron(cattle, field.list="count", main="Simulated distributions for main domesticates (50-yr bins)", legend=FALSE, lab.sp=2)
poly.chron(caprines, field.list="count", col.list="darkgreen", add=TRUE, legend=FALSE)
poly.chron(pigs, field.list="count", col.list="darkblue", add=TRUE, legend=FALSE)
poly.chron(chicken, field.list="count", col.list="grey30", add=TRUE, legend=FALSE)
legend("topright", legend=c("Cattle", "Caprines", "Pigs", "Chicken"), fill=c("darkred", "darkgreen", "darkblue", "grey30"), bty="n")

#4. Calculate distributions of context and samples
context.sim <- date.simulate(context.period, context.fields="SITE_C", reps=500)
#context <- aorist(context.period)
#write.csv(context, "CONTEXT-AORIST.csv", row.names=FALSE)
context <- data.table(read.csv("CONTEXT-AORIST-50.csv"))
con <- data.table(read.csv("CONTEXT-AORIST.csv"))

sample.sim <- date.simulate(sample.period, reps=500)
sample.count <- aorist(sample.period, bin.width=50)
write.csv(sample.count, "SAMPLE-COUNT-AORIST-50.csv", row.names=FALSE)
good.vols <- sample.period[is.na(WTS_VOL)==FALSE]
sample.vol <- aorist(good.vols, weight=good.vols$WTS_VOL, bin.width=50)
write.csv(sample.vol, "SAMPLE-VOL-AORIST-50-GOOD.csv", row.names=FALSE)
barplot(context$aorist)
with(context, barplot(aorist, ylab="Estimated frequency density", main="Aoristic sum of London contexts (50-yr bins)", col="grey70", names.arg=bin, las=2, cex.lab=1.2, cex.names=1.05, cex.axis=0.95, cex.main=1.2))
with(sample.count, barplot(aorist, ylab="Estimated frequency density", main="Aoristic sum of London environmental samples (50-yr bins)", col="grey70", names.arg=bin, las=2, cex.lab=1.2, cex.names=1.05, cex.axis=0.95, cex.main=1.2))
with(sample.vol, barplot(aorist, ylab="Estimated volume density", main="Aoristic sum of London environmental sample volumes (50-yr bins)", col="grey70", names.arg=bin, las=2, cex.lab=1.2, cex.names=1.05, cex.axis=0.95, cex.main=1.2))

#5. Use these to create better dummy sets
doms <- freq.simulate(zoo.period[SPECIES%in%c("BOS", "OVCA", "OVI", "CRA", "SUS", "CHIK")], weight=zoo.period[SPECIES%in%c("BOS", "OVCA", "OVI", "CRA", "SUS", "CHIK")]$FRAG_COUNT, probs=context$aorist, reps=1000)
poly.chron(doms, field.list=c("dummy", "count"), col.list=c("grey", "darkred"), main="Abundance of major domesticates (50-yr bins)")
lines.chron(doms, main="Abundance of major domesticates")

#all.fish.frags.samples <- freq.simulate(zoo.period[CLASS=="fish"], weight=zoo.period[CLASS=="fish"]$FRAG_COUNT, probs=sample.count$aorist, reps=1000, RoC=TRUE)
#poly.chron(all.fish.frags.samples, field.list=c("dummy", "count"), col.list=c("grey", "darkred"), main="Abundance of fish in samples (50-yr bins; dummy based on counts)")
#poly.chron(all.fish.frags.samples, field.list=c("RoC.dummy", "RoC.count"), col.list=c("grey", "darkred"), main="Simulated rates of change in fish from samples (dummy<-counts)", ylab="Estimated rate of change (%)")

all.fish.frags.vol <- freq.simulate(zoo.period[CLASS=="fish"], weight=zoo.period[CLASS=="fish"]$FRAG_COUNT, probs=sample.vol$aorist, reps=1000, RoC=TRUE)
poly.chron(all.fish.frags.vol, field.list=c("dummy", "count"), col.list=c("grey", "darkred"), main="Abundance of fish in samples (50-yr bins; dummy based on volumes)", lab.sp=2)
poly.chron(all.fish.frags.vol, field.list=c("RoC.dummy", "RoC.count"), col.list=c("grey", "darkred"), main="Simulated rates of change in fish from samples (dummy<-volume)", ylab="Estimated rate of change (%)")

#fresh.samples <- freq.simulate(fish.period[Fresh_Marine%in%c("Fresh", "Fresh/Marine")], weight=fish.period[Fresh_Marine%in%c("Fresh", "Fresh/Marine")]$Frag, probs=sample.count$aorist, reps=1000, RoC=TRUE)
#poly.chron(fresh.samples, field.list=c("dummy", "count"), col.list=c("grey", "darkgreen"), main="Freshwater/migratory fish in samples (dummy<-samples)", lab.sp=2)
#poly.chron(fresh.samples, field.list=c("RoC.dummy", "RoC.count"), col.list=c("grey", "darkgreen"), main="Freshwater/migratory fish in samples (dummy<-samples)", ylab="Estimated rate of change (%)", lab.sp=2)

fresh.vols <- freq.simulate(fish.period[Fresh_Marine%in%c("Fresh", "Fresh/Marine")], weight=fish.period[Fresh_Marine%in%c("Fresh", "Fresh/Marine")]$Frag, probs=sample.vol$aorist, reps=1000, RoC=TRUE)
poly.chron(fresh.vols, field.list=c("dummy", "count"), col.list=c("grey", "darkgreen"), main="Freshwater/migratory fish in samples (dummy<-volumes)", lab.sp=2)
poly.chron(fresh.vols, field.list=c("RoC.dummy", "RoC.count"), col.list=c("grey", "darkgreen"), main="Freshwater/migratory fish in samples (dummy<-volumes)", ylab="Estimated rate of change (%)", lab.sp=2)

#marine.samples <- freq.simulate(fish.period[Fresh_Marine=="Marine"], weight=fish.period[Fresh_Marine=="Marine"]$Frag, probs=sample.count$aorist, reps=1000, RoC=TRUE)
#poly.chron(marine.samples, field.list=c("dummy", "count"), col.list=c("grey", "darkblue"), main="Marine fish in samples (dummy<-samples)", lab.sp=2)
#poly.chron(marine.samples, field.list=c("RoC.dummy", "RoC.count"), col.list=c("grey", "darkblue"), main="Marine fish in samples (dummy<-samples)", ylab="Estimated rate of change (%)", lab.sp=2)

marine.vols <- freq.simulate(fish.period[Fresh_Marine=="Marine"], weight=fish.period[Fresh_Marine=="Marine"]$Frag, probs=sample.vol$aorist, reps=1000, RoC=TRUE)
poly.chron(marine.vols, field.list=c("dummy", "count"), col.list=c("grey", "darkblue"), main="Marine fish in samples (dummy<-volumes)", lab.sp=2)
poly.chron(marine.vols, field.list=c("RoC.dummy", "RoC.count"), col.list=c("grey", "darkblue"), main="Marine fish in samples (dummy<-volumes)", ylab="Estimated rate of change (%)", lab.sp=2)


#6 catch per unit
all.mammal.frags <- freq.simulate(zoo.period[CLASS=="mammal"], weight=zoo.period[CLASS=="mammal"]$FRAG_COUNT, bin.width=50, reps=1000, RoC=TRUE)

fresh <- cbind(fresh.samples[[2]], mam=all.mammal.frags[[2]]$V1)
fresh <- fresh[order(quantile)]
fresh <- cbind(fresh, litres=sample.vol$aorist, n=sample.count$aorist)
fresh[,cpl:=V1/litres]
fresh[,cps:=V1/n]
fresh[,cpm:=V1/mam]
poly.chron(fresh, field.list="count", col.list="darkgreen", value.field="cpl", main="Freshwater/migratory fish bones per litre", ylab="Estimated frequency per litre", legend=FALSE, lab.sp=2)
poly.chron(fresh, field.list="count", col.list="darkgreen", value.field="cps", main="Freshwater/migratory fish bones per sample", ylab="Estimated frequency per sample", legend=FALSE, lab.sp=2)
poly.chron(fresh, field.list="count", col.list="darkgreen", value.field="cpm", main="Freshwater/migratory fish per mammal bone", ylab="Estimated fish / mammal", legend=FALSE, lab.sp=2)


marine <- cbind(marine.samples[[2]], mam=all.mammal.frags[[2]]$V1)
marine <- marine[order(quantile)]
marine <- cbind(marine, litres=sample.vol$aorist, n=sample.count$aorist)
marine[,cpl:=V1/litres]
marine[,cps:=V1/n]
marine[,cpm:=V1/mam]
poly.chron(marine, field.list="count", col.list="darkblue", value.field="cpl", main="Marine fish bones per litre", ylab="Estimated frequency per litre", legend=FALSE, lab.sp=2)
poly.chron(marine, field.list="count", col.list="darkblue", value.field="cps", main="Marine fish bones per sample", ylab="Estimated frequency per sample", legend=FALSE, lab.sp=2)
poly.chron(marine, field.list="count", col.list="darkblue", value.field="cpm", main="Marine fish per mammal bone", ylab="Estimated fish / mammal", legend=FALSE, lab.sp=2)

col2rgb("darkgreen")
legend.blue <- rgb(0,0,139,126, maxColorValue=255)
legend.green <- rgb(0,100,0,126, maxColorValue=255)

axis.setup(marine, field.list=c("count"), main="Fish bones per litre", value.field="cpl", ylab="Estimated frequency per litre", lab.sp=1, type="summary")
poly.chron(fresh, field.list="count", col.list="darkgreen", value.field="cpl", add=TRUE, legend=FALSE)
poly.chron(marine, field.list="count", col.list="darkblue", value.field="cpl", add=TRUE, legend=FALSE)
legend("topleft", legend=c("Freshwater/migratory", "Marine"), fill=c(legend.green, legend.blue), bty="n")

axis.setup(marine, field.list=c("count"), main="Fish bones per sample", value.field="cps", ylab="Estimated frequency per sample", lab.sp=1, type="summary")
poly.chron(fresh, field.list="count", col.list="darkgreen", value.field="cps", add=TRUE, legend=FALSE)
poly.chron(marine, field.list="count", col.list="darkblue", value.field="cps", add=TRUE, legend=FALSE)
legend("topleft", legend=c("Freshwater/migratory", "Marine"), fill=c(legend.green, legend.blue), bty="n")

axis.setup(marine, field.list=c("count"), main="Fish bones per mammal", value.field="cpm", ylab="Estimated fish / mammal", lab.sp=1, type="summary")
poly.chron(fresh, field.list="count", col.list="darkgreen", value.field="cpm", add=TRUE, legend=FALSE)
poly.chron(marine, field.list="count", col.list="darkblue", value.field="cpm", add=TRUE, legend=FALSE)
legend("topleft", legend=c("Freshwater/migratory", "Marine"), fill=c(legend.green, legend.blue), bty="n")


####
Compare context types

#####
context <- aorist(context.period, bin.width=50)
write.csv(context, "CONTEXT-AORIST-50.csv", row.names=FALSE)




