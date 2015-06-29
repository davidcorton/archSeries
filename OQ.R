# Load required packages
library(data.table)
# Source custom functions
source("date_functions.R")
source("plot_functions.R")

# Read in data
period <- data.table(read.csv("period-CLEANED.csv"))
setkey(period, SITE_P)
context <- data.table(read.csv("context-RAW.csv"))
setkey(context, SITE_C)
sample <- data.table(read.csv("sample-CLEANED.csv"))
setkey(sample, SITE_S)
zoo <- data.table(read.csv("zoo-CLEANED.csv"))
setkey(zoo, SITE_C)
species <- data.table(read.csv("species_codes.csv"))
landuse <- data.table(read.csv("landuse_codes.csv"))
interpretation <- data.table(read.csv("int_codes.csv"))

# Merge species information and select fish
zoo <- merge(zoo, species, by="SPECIES", all.x=TRUE, all.y=FALSE)
fish <- zoo[CLASS=="fish"]
other.classes <- zoo[CLASS%in%c("mammal", "bird", "amphib", "reptile")]

# Merge context information
context.period <- merge(period[,list(SITE_P,Start,End)], context[,list(SITE_C,SITE_P,BASIC_INT,LU_INT)], by="SITE_P", all=FALSE)
fish.period <- merge(fish[!SAMPLE==0], context.period, by="SITE_C", all=FALSE)
other.period <- merge(other.classes[!SAMPLE==0], context.period, by="SITE_C", all=FALSE)
sample.period <- merge(sample, context.period, by="SITE_C", all=FALSE)
sample.period <- merge(sample.period, landuse, by="LU_INT", all.x=TRUE, all.y=FALSE)
sample.period <- merge(sample.period, interpretation, by="BASIC_INT", all.x=TRUE, all.y=FALSE)
rm(context, period, context.period, sample, fish, species)

# Summarise distribution of fish remains (FIGURE 3)
fish.sim <- freq.simulate(fish.period, weight=fish.period$FRAG_COUNT, bin.width=50, reps=2000)
axis.setup(fish.sim, lab.sp=2)
poly.chron(fish.sim, field.list="dummy", col.list="grey", add=TRUE, legend=FALSE)
box.chron(fish.sim, field.list="count", col.list="darkred", add=TRUE)

# Calculate distribution of sample volumes (Figure 4a)
sample.volumes <- aorist(sample.period[!is.na(WTS_VOL)], weight=sample.period[!is.na(WTS_VOL), WTS_VOL], bin.width=50)
with(sample.volumes, barplot(aorist, ylab="Estimated frequency density", names.arg=bin, las=2, cex.lab=1.2, cex.names=1.05, cex.axis=1.05, cex.main=1.2))

# Plot all fish with sample-based dummy (Figure 4b)
fish.by.volumes <- freq.simulate(fish.period, weight=fish.period$FRAG_COUNT, probs=sample.volumes$aorist, bin.width=50, reps=2000)
poly.chron(fish.by.volumes, field.list=c("dummy","count"), col.list=c("grey","darkred"), legend=FALSE)

# Demonstrate reason for wide confidence band (Figure 4c)
lines.chron(fish.by.volumes, col.list=c("darkred", "grey"), lab.sp=2, legend=FALSE)

# Plot by sample instead (Figure 4d)
sample.counts <- aorist(sample.period, bin.width=50)
fish.by.counts <- freq.simulate(fish.period, probs=sample.counts$aorist, bin.width=50, reps=2000)
lines.chron(fish.by.counts, col.list=c("darkred", "grey"), lab.sp=2, legend=FALSE)

# Compare freshwater and marine (FIGURE 5)
fresh <- fish.period[FRESH_MARINE%in%c("fresh", "fresh/marine", "migratory")]
marine <- fish.period[FRESH_MARINE=="marine"]
# Figure 5a - FW with dummy
fresh.vols <- freq.simulate(fresh, weight=fresh$FRAG_COUNT, probs=sample.volumes$aorist, reps=2000)
poly.chron(fresh.vols, field.list=c("dummy", "count"), col.list=c("grey", "darkgreen"), lab.sp=2)
# Figure 5b - Marine with dummy
marine.vols <- freq.simulate(marine, weight=marine$FRAG_COUNT, probs=sample.volumes$aorist, reps=2000)
poly.chron(marine.vols, field.list=c("dummy", "count"), col.list=c("grey", "darkblue"), lab.sp=2)
# Figure 5c - both by litre
fresh.cpue <- cpue(x=fresh, y=sample.period[!is.na(WTS_VOL)], fresh$FRAG_COUNT, small.n=c(2000,1000), sample.period[!is.na(WTS_VOL), WTS_VOL], quant.list=c(0.025,0.5,0.975),bin.width=50, reps=5000)
marine.cpue <- cpue(x=marine, y=sample.period[!is.na(WTS_VOL)], marine$FRAG_COUNT, sample.period[!is.na(WTS_VOL), WTS_VOL], quant.list=c(0.025,0.5,0.975),bin.width=50, reps=5000)
poly.chron(fresh.cpue, col="darkgreen", quant=c(0.025, 0.975), small.n=c("grey", "grey20"), small.n.op=200, ylim=1.1, legend=FALSE)
poly.chron(marine.cpue, col="darkblue", quant=c(0.025, 0.975), legend=FALSE, add=TRUE)
legend("topright", legend=c("Est. volume < 2000 litres", "Est. volume < 1000 litres"), fill=c("#BEBEBEBE", "#888888EE"), cex=0.9)
#Figure 5d - both by other taxa
fresh.rel <- cpue(x=fresh, y=other.period, fresh$FRAG_COUNT, other.period$FRAG_COUNT, small.n=c(500,250), quant.list=c(0.025,0.5,0.975),bin.width=50, reps=5000)
marine.rel <- cpue(x=marine, y=other.period, marine$FRAG_COUNT, other.period$FRAG_COUNT, quant.list=c(0.025,0.5,0.975),bin.width=50, reps=5000)
poly.chron(fresh.rel, col="darkgreen", quant=c(0.025, 0.975), small.n=c("grey","grey20"), small.n.op=200,ylim=1.8, legend=FALSE)
poly.chron(marine.rel, col="darkblue", quant=c(0.025, 0.975), legend=FALSE, add=TRUE)
legend("topright", legend=c("Est. other remains < 500", "Est. other remains < 250"), fill=c("#BEBEBEBE", "#888888EE"), cex=0.85)
#setnames(fresh.cpue[[1]], old="cpue", new="Freshwater/migratory")
#combine.cpue <- cbind(fresh.cpue[[1]], Marine=marine.cpue[[1]]$cpue)
#lines.chron(combine.cpue, field.list=field.list, col.list=c("darkgreen", "darkblue"), legend=FALSE, add=TRUE)

# Figure 6a
fresh.ubiq <- cpue(x=fresh, y=sample.period, quant.list=c(0.025,0.5,0.975), bin.width=50, reps=5000)
marine.ubiq <- cpue(x=marine, y=sample.period, quant.list=c(0.025,0.5,0.975), bin.width=50, reps=5000)
poly.chron(fresh.ubiq, col="darkgreen", quant=c(0.025, 0.975), legend=FALSE, ylab="Estimated ubiquity")
poly.chron(marine.ubiq, col="darkblue", quant=c(0.025, 0.975), legend=FALSE, add=TRUE)
# Figure 6b
fresh.hits.sample <- cpue(x=fresh, y=sample.period, UoA="SPECIES", quant.list=c(0.05,0.5,0.95), bin.width=50, reps=5000)
marine.hits.sample <- cpue(x=marine, y=sample.period, UoA="SPECIES", quant.list=c(0.05,0.5,0.95), bin.width=50, reps=5000)
poly.chron(fresh.hits.sample, col="darkgreen", quant=c(0.05, 0.95), legend=FALSE, ylab="Estimated mean nTaxa")
poly.chron(marine.hits.sample, col="darkblue", quant=c(0.05, 0.95), legend=FALSE, add=TRUE)


### Plot different FW taxa...
FW.taxa <- comp.simulate(fish.period, weight=fish.period$FRAG_COUNT, comp.field="FRESH_MARINE", comp.values=c("migratory", "fresh/marine", "fresh"), reps=2000, bin.width=50)
poly.chron(FW.taxa, col.list=c("goldenrod", "turquoise", "darkgreen"),legend=TRUE)
lines.chron(FW.taxa, col.list=c("goldenrod", "turquoise", "darkgreen"),legend=TRUE)

# Explore different context types (--> SI?)
context.types <- comp.simulate(sample.period[!is.na(INT_GROUP)], weight=sample.period[!is.na(INT_GROUP), WTS_VOL], comp.field="INT_GROUP", comp.values=c("pit", "ditch etc.", "dump", "occupation", "other cut"), reps=2000, bin.width=50)
ticks <- seq(1, length(names))
ylim=1000
plot(ticks, rep(1,length(ticks)), ylim=c(0,ylim), xlab="", xaxt="n", ylab="Estimated fish/other vertebrates", type="n")
axis(1, at=ticks, labels=names[ticks], las=2)
poly.chron(context.types, legend=FALSE, add=TRUE)
legend("topright", legend=c("Pit", "Ditch etc.", "Dump", "Occupation", "Other cut"), fill=c("darkred", "darkgreen", "blue", "grey", "goldenrod"), bty="n")


