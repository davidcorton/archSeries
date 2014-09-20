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
# Arguments: 'data' is a data table with at least two numeric columns called Start and End.
# If 'data' also has a column called taxon, this can be used with the 'species' argument to
# select the rows to be included in the simulation.
# 'species' is a single character value indicating which rows should be included in analysis
# It will be ignored if no taxon column is provided in 'data', and it defaults to NULL.
# 'start.date' is a single numeric value for the start of the time period to be analysed
# 'end.date' is a single numric value for the end end of the time period to be analysed
# 'bin.width' is a single numeric value setting the resolution of the analysis, in years
# 'rep' is the number of times the simulation will be run
# 'weight' is a numeric vector giving a weight for each context/entity, defaulting to 1
# Returns: a long-format data.table giving the sum of weight for each bin in each rep
# Also outputs: 'breaks', a numeric vector of breaks points,
# 'params', a character value summarising the arguments, for use in naming output files

date.simulate <- function(data, species=NULL, start.date=0, end.date=2000, bin.width=100, rep=100, weight=1) {
    require(data.table)
    data <- cbind(data, weight)
    if(length(species)>0 & "taxon" %in% colnames(data)) {data <- data[taxon==species,]}
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
# Arguments: 'probs' is normally a data.table (the output of an aorist call) consisting of
#   a numeric column ('aoristic.sum') to be used as relative probabilities, and a character
#   column ('labels') containing bin labels.
#   Alternatively, for a uniform dummy set, pass a uniform numeric vector whose length
#   matches the desired number of bins - e.g. rep(1, 100), where 100 bins are required.
# 'weight' is a numeric vector (or data frame/data.table) representing (weighted) instances
# to be simulated. If given an additional character column called taxon, this can be used
# to select rows for analysis using the 'species' argument.
# 'species' is a single character value indicating rows to be included in analysis. Ignored
# unless 'weight' has a taxon colum. Defaults to NULL.
# 'start.date' and 'end.date' are single numeric values. Only required where a single vector
#   is passed to 'probs' and the range under study is not 0-2000AD.

dummy.simulate <- function(probs, weight, species=NULL, start.date=0, end.date=2000, rep=500) {
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
    if(length(species)>0 & "taxon" %in% colnames(weight)) {weight <- weight[taxon==species,]}
    a.sum <- sum(probs$aoristic.sum)
    a.breaks <- c(0, cumsum(probs$aoristic.sum))
    dummy <- cbind(rep(1:rep, each=nrow(dummy)), dummy)
    setnames(dummy, old=1, new="rep.no")
    dummy[,bin:=runif(nrow(dummy), 0, a.sum)]
    dummy[,bin:=cut(bin, a.breaks, labels=probs$labels)]
    dummy <- dummy[is.na(bin)==FALSE, sum(weight), by=list(rep.no,bin)]
    dummy[order(rep.no, bin)]
}

# Define function that performs both 'real' and dummy simulation on target bone data
# Arguments: 'data' is a data.frame or data.table with columns including Start, End, and
# Frag. If additional factor columns are provided, these can be used to filter the data
# using the 'filter.field' and 'filter.values' arguments, below. The function saves both
# full and summary simulation results to .csv files, and also returns the latter.
# 'probs' is normally a data.table (the output of an aorist call) consisting of
#   a numeric column ('aoristic.sum') to be used as relative probabilities, and a character
#   column ('labels') containing bin labels.
#   Alternatively, for a uniform dummy set, pass a uniform numeric vector whose length
#   matches the desired number of bins - e.g. rep(1, 100), where 100 bins are required.
# 'filter.field' is a single character value denoting the name of a column that will be
#   used to filter the data. This defaults to "Species" but will be ignored unless 
#   'filter.values' is set.
# 'filter.values'is a character vector containing all values of the filter column that will
#   be included in the analysis. Defaults to NULL.
# 'quant.list' is a numeric vector of quantiles to be included in the summary output.
# 'start.date' and 'end.date' are the chronological limits of the overall analysis.
# 'rep' is the number of times that both 'real' and dummy simulations will be repeated.

freq.simulate <- function(data, probs, filter.field="Species", filter.values=NULL, quant.list=c(0.025,0.25,0.5,0.75,0.975), start.date=0, end.date=2000, rep=1000) {
    require(data.table)
    require(reshape2)
    data <- data.table(data)  #just in case it isn't already in this format
    if(is.null(filter.values)==FALSE) {  #if criteria have been provided...
        setnames(data, old=filter.field, new="FILTER")   #...sets the filter column...
        data <- data[FILTER %in% filter.values,]  #...and applies the criteria
    } else {filter.values <- "ALL"}
    data <- data[End >= start.date & Start <= end.date]  #drops records outside the date range FROM BOTH SIMULATION SETS
    bin.width <- (end.date-start.date)/nrow(probs)  #sets bin widths (and hence no. of bins) to match the calibration dataset
    
    # set up list of all bins and rep no.s
    breaks <- seq(start.date, end.date, bin.width)
    labels <- numeric(nrow(probs))
    for(i in 1:length(labels)) {
        labels[i] <- paste(breaks[i], breaks[i+1], sep="-")
    } 
    frame <- data.table(rep(1:rep, each=length(labels)), rep(1:length(labels), rep), rep(labels, rep))
    setnames(frame, old=c("V1", "V2", "V3"), new=c("rep.no", "bin.no", "bin"))
    
    # simulate from real data
    real <- date.simulate(data=data[,list(Start, End)], weight=data[,Frag], bin.width=bin.width, start.date=start.date, end.date=end.date, rep=rep)
    setnames(real, old="V1", new="real")
    
    # simulate dummy set
    dummy <- dummy.simulate(probs=probs, weight=data[,Frag], start.date=start.date, end.date=end.date, rep=rep)    
    setnames(dummy, old="V1", new="dummy")
    
    # merge the above three data.tables together
    results <- merge(frame, dummy, by=c("rep.no", "bin"), all=TRUE)
    results <- merge(results, real, by=c("rep.no", "bin"), all=TRUE)
    results[is.na(real)==TRUE, real:=0] 
    results[is.na(dummy)==TRUE, dummy:=0]
    
    # calculate rate of change variables
    for(i in 1:(nrow(results)-1)) {
        results[i,ROC.real:=(results[i+1,real]-results[i,real])/bin.width]
        results[i,ROC.dummy:=(results[i+1,dummy]-results[i,dummy])/bin.width]
    }
    results[bin==labels[length(labels)], ROC.real:=NA]
    results[bin==labels[length(labels)], ROC.dummy:=NA]
    
    # save full dataset
    write.csv(results, paste("TEST_", filter.values[1], "_simulated_by_period", params, ".csv",sep=""), row.names=FALSE)

    # create summary dataset
    real.summary <- results[,quantile(real, probs=quant.list, na.rm=TRUE), by=bin]
    real.summary[,id:=paste(rep("real", length(quant.list)), quant.list, sep="_")]
    real.summary <- dcast.data.table(real.summary, bin ~ id, value.var="V1")
   
    dummy.summary <- results[,quantile(dummy, probs=quant.list, na.rm=TRUE), by=bin]
    dummy.summary[,id:=paste(rep("dummy", length(quant.list)), quant.list, sep="_")]
    dummy.summary <- dcast.data.table(dummy.summary, bin ~ id, value.var="V1")
    
    ROC.real.summary <- results[, quantile(ROC.real, probs=quant.list, na.rm=TRUE), by=bin]
    ROC.real.summary[,id:=paste(rep("ROC.real", length(quant.list)), quant.list, sep="_")]
    ROC.real.summary <- dcast.data.table(ROC.real.summary, bin ~ id, value.var="V1")
    
    ROC.dummy.summary <- results[, quantile(ROC.dummy, probs=quant.list, na.rm=TRUE), by=bin]
    ROC.dummy.summary[,id:=paste(rep("ROC.dummy", length(quant.list)), quant.list, sep="_")]
    ROC.dummy.summary <- dcast.data.table(ROC.dummy.summary, bin ~ id, value.var="V1")
    
    summary <- cbind(real.summary, dummy.summary, ROC.real.summary, ROC.dummy.summary)
    
    #save summary dataset
    write.csv(summary, paste("TEST_summary_", filter.values[1], "_simulated_by_period", params, ".csv", sep=""), row.names=FALSE)

    #output list with full and summary datasets
    list(results, summary)
}

anat.simulate <- function(data, probs, filter.field="group", filter.values=c("cranial", "postcranial"), quant.list=c(0.025,0.25,0.5,0.75,0.975), start.date=0, end.date=2000, rep=1000) {
    require(data.table)
    require(reshape2)
    data <- data.table(data) 
    setnames(data, old=filter.field, new="FILTER")
    data <- data[End >= start.date & Start <= end.date]
    bin.width <- (end.date-start.date)/length(probs)
    comparison <<- filter.values
    
    # set up list of all bins and rep no.s
    breaks <- seq(start.date, end.date, bin.width)
    labels <- numeric(length(probs))
    for(i in 1:length(labels)) {
        labels[i] <- paste(breaks[i], breaks[i+1], sep="-")
    } 
    frame <- data.table(rep(1:rep, each=length(labels)), rep(1:length(labels), rep), as.character(rep(labels, rep)))
    setnames(frame, old=c("V1", "V2", "V3"), new=c("rep.no", "bin.no", "bin"))
    
    # run simulation on each group
    sim.1 <- date.simulate(data[FILTER==filter.values[1],list(Start, End)], start.date=start.date, end.date=end.date, bin.width=bin.width, rep=rep, weight=data[FILTER==filter.values[1],Frag])
    sim.2 <- date.simulate(data[FILTER==filter.values[2],list(Start, End)], start.date=start.date, end.date=end.date, bin.width=bin.width, rep=rep, weight=data[FILTER==filter.values[2],Frag]) 
    real <- merge(sim.1, sim.2, by=c("rep.no", "bin"), all=TRUE)
    rm(sim.1, sim.2)
    setnames(real, old=c("V1.x", "V1.y"), new=c("real.1", "real.2"))
    
    # run dummy simulation on each group
    dummy.1 <- dummy.simulate(probs, data[FILTER==filter.values[1], Frag], start.date=start.date, end.date=end.date, rep=rep)
    dummy.2 <- dummy.simulate(probs, data[FILTER==filter.values[2], Frag], start.date=start.date, end.date=end.date, rep=rep) 
    dummy <- merge(dummy.1, dummy.2, by=c("rep.no", "bin"), all=TRUE)
    rm(dummy.1, dummy.2)
    setnames(dummy, old=c("V1.x", "V1.y"), new=c("dummy.1", "dummy.2"))
        
    # combine results
    combined <- merge(real, dummy, by=c("rep.no", "bin"), all=TRUE)
    combined <- merge(combined, frame, by=c("rep.no", "bin"), all=TRUE)
    
    # replace NAs with 0s
    combined[is.na(real.1)==TRUE, real.1:=0]
    combined[is.na(real.2)==TRUE, real.2:=0]
    combined[is.na(dummy.1)==TRUE, dummy.1:=0]
    combined[is.na(dummy.2)==TRUE, dummy.2:=0]
    
    # calculate differences
    combined[,diff:=(real.1-real.2)/(real.1+real.2)]
    combined[,dummy.diff:=dummy.1-dummy.2]
    
    # calculate rates of change
    for(i in 1:(nrow(combined)-1)) {
        combined[i,ROC.real.1:=(combined[i+1,real.1]-combined[i,real.1])/bin.width]
        combined[i,ROC.dummy.1:=(combined[i+1,dummy.1]-combined[i,dummy.1])/bin.width]
        combined[i,ROC.real.2:=(combined[i+1,real.2]-combined[i,real.2])/bin.width]
        combined[i,ROC.dummy.2:=(combined[i+1,dummy.2]-combined[i,dummy.2])/bin.width]
        combined[i,ROC.diff:=(combined[i+1,diff]-combined[i,diff])/bin.width]
        combined[i,ROC.dummy.diff:=(combined[i+1,dummy.diff]-combined[i,dummy.diff])/bin.width]
        if(i/500 == round(i/500)) {print(paste(i/nrow(combined)*100, "percent complete"))}
    }
    combined[bin==labels[length(labels)], ROC.real.1:=NA]
    combined[bin==labels[length(labels)], ROC.dummy.1:=NA]
    combined[bin==labels[length(labels)], ROC.real.2:=NA]    
    combined[bin==labels[length(labels)], ROC.dummy.2:=NA]
    combined[bin==labels[length(labels)], ROC.diff:=NA]
    combined[bin==labels[length(labels)], ROC.dummy.diff:=NA]
    
    
    # create summary dataset
    real.1.summary <- combined[,quantile(real.1, probs=quant.list, na.rm=TRUE), by=bin]
    real.1.summary[,id:=paste(rep("real.1", length(quant.list)), quant.list, sep="_")]
    real.1.summary <- dcast.data.table(real.1.summary, bin ~ id, value.var="V1")
    
    real.2.summary <- combined[,quantile(real.2, probs=quant.list, na.rm=TRUE), by=bin]
    real.2.summary[,id:=paste(rep("real.2", length(quant.list)), quant.list, sep="_")]
    real.2.summary <- dcast.data.table(real.2.summary, bin ~ id, value.var="V1")
    
    dummy.1.summary <- combined[,quantile(dummy.1, probs=quant.list, na.rm=TRUE), by=bin]
    dummy.1.summary[,id:=paste(rep("dummy.1", length(quant.list)), quant.list, sep="_")]
    dummy.1.summary <- dcast.data.table(dummy.1.summary, bin ~ id, value.var="V1")
  
    dummy.2.summary <- combined[,quantile(dummy.2, probs=quant.list, na.rm=TRUE), by=bin]
    dummy.2.summary[,id:=paste(rep("dummy.2", length(quant.list)), quant.list, sep="_")]
    dummy.2.summary <- dcast.data.table(dummy.2.summary, bin ~ id, value.var="V1")
    
    ROC.dummy.1.summary <- combined[, quantile(ROC.dummy.1, probs=quant.list, na.rm=TRUE), by=bin]
    ROC.dummy.1.summary[,id:=paste(rep("ROC.dummy.1", length(quant.list)), quant.list, sep="_")]
    ROC.dummy.1.summary <- dcast.data.table(ROC.dummy.1.summary, bin ~ id, value.var="V1")
    
    ROC.dummy.2.summary <- combined[, quantile(ROC.dummy.2, probs=quant.list, na.rm=TRUE), by=bin]
    ROC.dummy.2.summary[,id:=paste(rep("ROC.dummy.2", length(quant.list)), quant.list, sep="_")]
    ROC.dummy.2.summary <- dcast.data.table(ROC.dummy.2.summary, bin ~ id, value.var="V1")
    
    ROC.dummy.diff.summary <- combined[, quantile(ROC.dummy.diff, probs=quant.list, na.rm=TRUE), by=bin]
    ROC.dummy.diff.summary[,id:=paste(rep("ROC.dummy.diff", length(quant.list)), quant.list, sep="_")]
    ROC.dummy.diff.summary <- dcast.data.table(ROC.dummy.diff.summary, bin ~ id, value.var="V1")
    
    summary <- cbind(real.1.summary, real.2.summary, dummy.1.summary, dummy.2.summary, ROC.dummy.1.summary, ROC.dummy.2.summary, ROC.dummy.diff.summary)
    for(i in 1:6){summary[,bin:=NULL]}
    summary <- merge(combined[rep.no==1, list(bin, bin.no)], summary, by="bin")
    
    list(combined[order(rep.no, bin.no)], summary[order(bin.no)])
}

