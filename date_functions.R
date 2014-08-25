# Define function to calculate aoristic sum
# Arguments: 'data' is a data table with two numeric columns called Start and End
# 'start.date' is a single numeric value for the start of the time period to be analysed
# 'end.date' is a single numric value for the end end of the time period to be analysed
# 'bin.width' is a single numeric value setting the resolution of the analysis, in years
# 'weight' is a numeric vector giving a weight for each context/entity
# Returns: a two-column data table with the aoristic sum itself (numeric) and bin labels (character)
# Also outputs: 'breaks', a numeric vector of breaks points,
# 'params', a character value summarising the arguments, for use in naming output files

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
    breaks <<- seq(start.date, end.date, bin.width)
    labels <- numeric(length(breaks)-1)
    for(i in 1:length(labels)) {
        labels[i] <- paste(breaks[i], breaks[i+1], sep="-")
    }
    params <<- paste("_", start.date, "-", end.date, "_by_", bin.width, sep="")
    data.table(aoristic.sum[1:length(labels)], labels)
}

# Define function to simulate distribution of dates
# Arguments: 'data' is a data table with two numeric columns called Start and End
# 'start.date' is a single numeric value for the start of the time period to be analysed
# 'end.date' is a single numric value for the end end of the time period to be analysed
# 'bin.width' is a single numeric value setting the resolution of the analysis, in years
# 'rep' is the number of times the simulation will be run
# 'weight' is a numeric vector giving a weight for each context/entity, defaulting to 1
# Returns: a long-format data.table giving the sum of weight for each bin in each rep
# Also outputs: 'breaks', a numeric vector of breaks points,
# 'params', a character value summarising the arguments, for use in naming output files

date.simulate <- function(data, start.date=0, end.date=2000, bin.width=100, rep=100, weight=1) {
    require(data.table)
    data <- cbind(data, weight)
    data <- data[End >= start.date & Start <= end.date]
    breaks <<- seq(start.date, end.date, bin.width)
    params <<- paste("_", start.date, "-", end.date, "_by_", bin.width, "_x", rep, sep="")
    labels <- numeric(length(breaks)-1)
    for(i in 1:length(labels)) {
        labels[i] <- paste(breaks[i], breaks[i+1], sep="-")
    }
    data <- cbind(rep(1:rep, each=nrow(data)), data)
    data[,bin:={x<-runif(nrow(data)); (x*(data[,End]-data[,Start]))+data[,Start]}]
    data[,bin:=cut(bin,breaks,labels=labels)]
    data <- data[is.na(bin)==FALSE, sum(weight), by=list(V1,bin)]
    setnames(data, old=1, new="rep.no")
    data[order(rep.no, bin)]
}

# Define function to simulate a dummy set by sampling from within an aoristic sum output
# Arguments: 'probs' is a normally a data.table (the output of an aorist call) consisting of
#   a numeric column ('aoristic.sum') to be used as relative probabilities, and a character
#   column ('labels') containing bin labels.
#   Alternatively, for a uniform dummy set, pass a uniform numeric vector whose length
#   matches the desired number of bins - e.g. rep(1, 100), where 100 bins are required.
# 'weight' is a numeric vector represented (weighted) instances to be simulated
# 'start.date' and 'end.date' are single numeric values. Only required where a single vector
#   is passed to 'probs' and the range under study is not 0-2000AD.

dummy.simulate <- function(probs, weight, start.date=0, end.date=2000, rep=500) {
    require(data.table)
    probs <- data.table(probs)
    if(ncol(probs)==1) {
        bin.width <- (end.date-start.date)/nrow(probs)
        breaks <- seq(start.date, end.date, bin.width)
        labels <- numeric(nrow(probs))
        for(i in 1:length(labels)) {
            labels[i] <- paste(breaks[i], breaks[i+1], sep="-")
        } 
        probs[,labels:=labels]
    }
    setnames(probs, c(1,2), c("aoristic.sum", "labels"))
    dummy <- data.table(weight)
    a.sum <- sum(probs$aoristic.sum)
    a.breaks <- c(0, cumsum(probs$aoristic.sum))
    dummy <- cbind(rep(1:rep, each=nrow(dummy)), dummy)
    setnames(dummy, old=1, new="rep.no")
    dummy[,bin:=runif(nrow(dummy), 0, a.sum)]
    dummy[,bin:=cut(bin, a.breaks, labels=probs$labels)]
    dummy <- dummy[is.na(bin)==FALSE, sum(weight), by=list(rep.no,bin)]
    dummy[order(rep.no, bin)]
}



