# Define function to calculate aoristic sum

aorist <- function(data, start.date=0, end.date=2000, bin.width=100, weight=1, progress=TRUE) { 
    require(data.table)
    aoristic.sum <- numeric(((end.date-start.date)/bin.width)) #sets up empty output variable of correct length
    data <- cbind(data, weight) #appends weights to list of date ranges, recycling if necessary (e.g. for uniform weight)
    data <- data[End >= start.date & Start <= end.date] #excludes ranges that fall entirely outside the study period
    data[,first.full:=ceiling((Start-start.date)/bin.width)+1] #finds number of first bin fully after Start
    data[,last.full:=floor((End-start.date)/bin.width)] #finds the number of the last bin fully before End
    data[,lead.in:=(ceiling((Start-start.date)/bin.width)-((Start-start.date)/bin.width))*bin.width]
    data[,lead.out:=(((End-start.date)/bin.width)-floor((End-start.date)/bin.width))*bin.width]
    data[,diff:=last.full-first.full]
    data[,duration:={ifelse(diff==-2, lead.in+lead.out, End-Start)}] #finds duration of range (special case where w/in single bin)
    data[,full.prob:=(bin.width/duration)*weight] #sets prob. mass to be assigned to bins fully within range (then applies weight)
    data[,in.prob:=(lead.in/duration)*weight] #sets prob. mass to be assigned to the lead-in bin (then applies weight)
    data[,out.prob:=(lead.out/duration)*weight] #sets prob. mass to be assigned to the lead-out bin (then applies weight)
    for(i in 1:nrow(data)) { #cycles through data, first assigning lead-in & lead-out probs to relevant bins in output variable
        aoristic.sum[data[i,first.full]-1] <- aoristic.sum[data[i,first.full]-1] + data[i,in.prob]
        aoristic.sum[data[i,last.full]+1] <- aoristic.sum[data[i,last.full]+1] + data[i,out.prob]
        if(data[i,diff] >= 0) { #for ranges spanning complete bins, assigns relevant probabilities to those bins 
            for(j in data[i,first.full]:data[i,last.full]) {
                aoristic.sum[j] <- aoristic.sum[j] + data[i,full.prob]
            }
        }
        if(progress==TRUE & i/1000 == round(i/1000)) {print(paste(i/nrow(data)*100, "percent complete"))} #progress monitor
    }
    breaks <<- seq(start.date, end.date, bin.width) #saves vector of breaks
    labels <- numeric(length(breaks)-1)
    for(i in 1:length(labels)) {
        labels[i] <- paste(breaks[i], breaks[i+1], sep="-") #sets bin labels
    }
    params <<- paste("_", start.date, "-", end.date, "_by_", bin.width, sep="") #saves character value with key parameters 
    data.table(aoristic.sum[1:length(labels)], labels) #returns aoristic sum with labels appended
}

# Define function to simulate distribution of dates

date.simulate <- function(data, filter=NULL, start.date=0, end.date=2000, bin.width=100, reps=100, weight=1) {
    require(data.table)
    data <- cbind(data, weight) #appends weights to list of date ranges, recycling if necessary (e.g. for uniform weight)
    if(length(filter)>0 & "group" %in% colnames(data)) {data <- data[group%in%filter,]} #filters data, if appropriate
    data <- data[End >= start.date & Start <= end.date] #excludes ranges that fall entirely outside the study period
    breaks <<- seq(start.date, end.date, bin.width) #sets breaks and saves them externally
    labels <- numeric(length(breaks)-1)
    for(i in 1:length(labels)) {
        labels[i] <- paste(breaks[i], breaks[i+1], sep="-") #sets bin labels
    }
    params <<- paste("_", start.date, "-", end.date, "_by_", bin.width, "_x", reps, sep="") #saves char value with key parameters 
    rep.no <- rep(1:reps, each=nrow(data))
    data <- cbind(rep.no, data) #recycles input data 'reps' times to provide frame for simulation
    data[,sim:={x<-runif(nrow(data)); (x*(data[,End]-data[,Start]))+data[,Start]}] #simulates a date for each row
    data[,bin.no:=cut(sim,breaks, labels=FALSE)] #records the relevant bin for each simulated date
    data[,bin:=cut(sim,breaks,labels=labels)] #records the relevant bin labels
    data <- data[is.na(bin)==FALSE, j=list(count=sum(weight)), by=list(rep.no,bin,bin.no)] #sums weights by bin and rep number
    data[order(rep.no, bin.no)]
}

# Define function to simulate a dummy set by sampling from within a specified distribution

dummy.simulate <- function(weight, probs=1, breaks=NULL, filter=NULL, start.date=0, end.date=2000, bin.width=100, reps=100) {
    require(data.table)
    
    if(breaks==NULL) {breaks <<- seq(start.date, end.date, bin.width)} #if breaks not specified, sets them based on other arguments
    labels <- numeric(length(breaks)-1)
    for(i in 1:length(labels)) {
        labels[i] <- paste(breaks[i], breaks[i+1], sep="-") #sets bin labels based on breaks
    } 
    probs <- cbind(probs, labels) #append labels to relative probs, recycling the latter if necessary
    
    dummy <- data.table(weight) #set 
    if(nrow(dummy)==1) {dummy <- rep(1, dummy)} #if weight is a single value, use as number of entities

    if(length(filter)>0 & "group"%in%colnames(dummy)) {dummy <- dummy[group==filter,]} #filter if appropriate
    rep.no <- rep(1:reps, each=nrow(dummy))
    dummy <- cbind(rep.no, dummy) #recycles input data 'reps' times to provide frame for simulation 
    p.sum <- sum(probs$probs) #'stacks up' all the relative probabilities
    p.breaks <- c(0, cumsum(probs$probs)) #uses cumulative sum of relative probabilities to set breaks
    
    dummy[,sim:=runif(nrow(dummy), 0, p.sum)] #samples from within p.sum
    
    
    dummy[,bin.no:=cut(sim,breaks, labels=FALSE)] #records the relevant bin for each simulated date
    dummy[,bin:=cut(sim, p.breaks, labels=probs$labels)] #finds the relevant bin labels
    dummy <- dummy[is.na(bin)==FALSE, j=list(count=sum(weight)), by=list(rep.no,bin,bin.no)] #sums weights by bin and rep number
    dummy[order(rep.no, bin.no)]
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
# 'ROC' is a logical value indicating whether rates-of-change should be calculated and
#   appended to the output
# 'start.date' and 'end.date' are the chronological limits of the overall analysis.
# 'rep' is the number of times that both 'real' and dummy simulations will be repeated.

freq.simulate <- function(data, probs, filter.field="Species", filter.values=NULL, quant.list=c(0.025,0.25,0.5,0.75,0.975), ROC=FALSE, start.date=0, end.date=2000, rep=100) {
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
    if(ROC==TRUE) {
        for(i in 1:(nrow(results)-1)) {
            results[i,ROC.real:=(results[i+1,real]-results[i,real])/bin.width]
            results[i,ROC.dummy:=(results[i+1,dummy]-results[i,dummy])/bin.width]
        }
        results[bin==labels[length(labels)], ROC.real:=NA]
        results[bin==labels[length(labels)], ROC.dummy:=NA]
    }
        
    # save full dataset
    write.csv(results, paste("TEST_", filter.values[1], "_simulated_by_period", params, ".csv",sep=""), row.names=FALSE)

    # create summary dataset
    real.summary <- results[,quantile(real, probs=quant.list, na.rm=TRUE), by=bin]
    real.summary[,id:=paste(rep("real", length(quant.list)), quant.list, sep="_")]
    real.summary <- dcast.data.table(real.summary, bin ~ id, value.var="V1")
   
    dummy.summary <- results[,quantile(dummy, probs=quant.list, na.rm=TRUE), by=bin]
    dummy.summary[,id:=paste(rep("dummy", length(quant.list)), quant.list, sep="_")]
    dummy.summary <- dcast.data.table(dummy.summary, bin ~ id, value.var="V1")
    
    summary <- cbind(real.summary, dummy.summary)
   
    if(ROC==TRUE) {
        ROC.real.summary <- results[, quantile(ROC.real, probs=quant.list, na.rm=TRUE), by=bin]
        ROC.real.summary[,id:=paste(rep("ROC.real", length(quant.list)), quant.list, sep="_")]
        ROC.real.summary <- dcast.data.table(ROC.real.summary, bin ~ id, value.var="V1")
        ROC.dummy.summary <- results[, quantile(ROC.dummy, probs=quant.list, na.rm=TRUE), by=bin]
        ROC.dummy.summary[,id:=paste(rep("ROC.dummy", length(quant.list)), quant.list, sep="_")]
        ROC.dummy.summary <- dcast.data.table(ROC.dummy.summary, bin ~ id, value.var="V1")
        summary <- cbind(summary, ROC.real.summary, ROC.dummy.summary)
    }
    
    #save summary dataset
    write.csv(summary, paste("TEST_summary_", filter.values[1], "_simulated_by_period", params, ".csv", sep=""), row.names=FALSE)

    #output list with full and summary datasets
    list(results, summary)
}

anat.simulate <- function(data, probs, filter.field="group", filter.values=c("cranial", "postcranial"), quant.list=c(0.025,0.25,0.5,0.75,0.975), ROC=FALSE, start.date=0, end.date=2000, rep=1000) {
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
    combined[,dummy.diff:=(dummy.1-dummy.2)/(real.1+real.2)]
    
    # calculate rates of change
    if(ROC==TRUE) {
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
    }
    
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
    
    summary <- cbind(real.1.summary, real.2.summary, dummy.1.summary, dummy.2.summary)
    binvarcount <- 3
    
    if(ROC==TRUE) {
        ROC.dummy.1.summary <- combined[, quantile(ROC.dummy.1, probs=quant.list, na.rm=TRUE), by=bin]
        ROC.dummy.1.summary[,id:=paste(rep("ROC.dummy.1", length(quant.list)), quant.list, sep="_")]
        ROC.dummy.1.summary <- dcast.data.table(ROC.dummy.1.summary, bin ~ id, value.var="V1")
        ROC.dummy.2.summary <- combined[, quantile(ROC.dummy.2, probs=quant.list, na.rm=TRUE), by=bin]
        ROC.dummy.2.summary[,id:=paste(rep("ROC.dummy.2", length(quant.list)), quant.list, sep="_")]
        ROC.dummy.2.summary <- dcast.data.table(ROC.dummy.2.summary, bin ~ id, value.var="V1")
        ROC.dummy.diff.summary <- combined[, quantile(ROC.dummy.diff, probs=quant.list, na.rm=TRUE), by=bin]
        ROC.dummy.diff.summary[,id:=paste(rep("ROC.dummy.diff", length(quant.list)), quant.list, sep="_")]
        ROC.dummy.diff.summary <- dcast.data.table(ROC.dummy.diff.summary, bin ~ id, value.var="V1")
        summary <- cbind(summary, ROC.dummy.1.summary, ROC.dummy.2.summary, ROC.dummy.diff.summary)
        binvarcount <- 6
    }
    
    for(i in 1:binvarcount){summary[,bin:=NULL]}
    summary <- merge(combined[rep.no==1, list(bin, bin.no)], summary, by="bin")
    
    list(combined[order(rep.no, bin.no)], summary[order(bin.no)])
}

