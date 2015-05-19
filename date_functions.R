# Define function to calculate aoristic sum

aorist <- function(data, start.date=0, end.date=2000, bin.width=100, weight=1, progress=TRUE) { 
    #Load required package
    require(data.table)
    
    #Tidies up input data
    data <- cbind(data, weight) #appends weights to list of date ranges, recycling if necessary (e.g. for uniform weight)
    data <- data[End >= start.date & Start <= end.date] #excludes ranges that fall entirely outside the study period
    
    #Detect first full bin after start & last full bin before end (yes, this seems wrong for short ranges - bear with me.)
    data[,first.full:=ceiling((Start-start.date)/bin.width)+1] #finds number of first bin fully after Start
    data[,last.full:=floor((End-start.date)/bin.width)] #finds the number of the last bin fully before End
    
    #Find extent of overlap with lead-in and lead-out bins for each date range
    data[,lead.in:=(ceiling((Start-start.date)/bin.width)-((Start-start.date)/bin.width))*bin.width]
    data[,lead.out:=(((End-start.date)/bin.width)-floor((End-start.date)/bin.width))*bin.width]
    
    #Work out probability mass to be assigned to whole bins and to lead-in and lead-out bins
    data[,diff:=last.full-first.full]
    data[,duration:={ifelse(diff==-2, lead.in+lead.out, End-Start)}] #finds duration of range (special case where w/in single bin)
    data[,full.prob:=(bin.width/duration)*weight] #sets prob. mass to be assigned to bins fully within range (then applies weight)
    data[,in.prob:=(lead.in/duration)*weight] #sets prob. mass to be assigned to the lead-in bin (then applies weight)
    data[,out.prob:=(lead.out/duration)*weight] #sets prob. mass to be assigned to the lead-out bin (then applies weight)
    
    #Cycle through data assigning probs to bins (this is the slow bit)
    aoristic.sum <- numeric(((end.date-start.date)/bin.width)) #sets up empty output variable of correct length
    for(i in 1:nrow(data)) { #First assigns lead-in & lead-out probs to relevant bins
        aoristic.sum[data[i,first.full]-1] <- aoristic.sum[data[i,first.full]-1] + data[i,in.prob]
        aoristic.sum[data[i,last.full]+1] <- aoristic.sum[data[i,last.full]+1] + data[i,out.prob]
        if(data[i,diff] >= 0) { #Then for ranges spanning complete bins, assigns relevant probs to those bins 
            for(j in data[i,first.full]:data[i,last.full]) {
                aoristic.sum[j] <- aoristic.sum[j] + data[i,full.prob]
            }
        }
        if(progress==TRUE & i/1000 == round(i/1000)) {print(paste(i/nrow(data)*100, "percent complete"))} #progress monitor
    }
    
    #Set up breaks and labels
    breaks <<- seq(start.date, end.date, bin.width) #creates and saves vector of breaks
    labels <- numeric(0)
    for(i in 1:(length(breaks)-1)) {
        labels[i] <- paste(breaks[i], breaks[i+1], sep="-") #sets bin labels
    }
    params <<- paste("_", start.date, "-", end.date, "_by_", bin.width, sep="") #saves character value with key parameters 
    
    #Prepare and return results table
    data.table(aorist=aoristic.sum[1:length(labels)], bin=labels, bin.no=1:length(labels)) #returns aoristic sum with labels appended
}

# Define function to simulate distribution of dates

date.simulate <- function(data, weight=1, filter=NULL, start.date=0, end.date=2000, bin.width=100, reps=100, RoC=FALSE) {
    #Load required package
    require(data.table)
    
    #Tidy up input data and apply filters
    data <- data.table(cbind(data, weight)) #appends weights to list of date ranges, recycling if necessary (e.g. for uniform weight)
    if(length(filter)>0 & "group" %in% colnames(data)) {data <- data[group%in%filter,]} #filters data, if appropriate
    data <- data[End >= start.date & Start <= end.date] #excludes ranges that fall entirely outside the study period
    
    #Set up breaks and labels
    breaks <<- seq(start.date, end.date, bin.width) #sets breaks and saves them externally
    labels <- numeric(0)
    for(i in 1:(length(breaks)-1)) {
        labels[i] <- paste(breaks[i], breaks[i+1], sep="-") #sets bin labels
    }
    params <<- paste("_", start.date, "-", end.date, "_by_", bin.width, "_x", reps, sep="") #saves char value with key parameters 
       
    #Perform simulation
    rep.no <- rep(1:reps, each=nrow(data))
    data <- cbind(rep.no, data) #recycles input data 'reps' times to provide frame for simulation 
    data[,sim:={x<-runif(nrow(data)); (x*(data[,End]-data[,Start]))+data[,Start]}] #simulates a date for each row
    data[,bin:=cut(sim,breaks,labels=labels)] #finds the relevant bin for each simulated date
    data <- data[is.na(bin)==FALSE, j=list(count=sum(as.numeric(weight))), by=list(rep.no,bin)] #sums weights by bin and rep number
    
    #Prepare results table
    frame <- data.table(rep.no=rep(1:reps, each=length(labels)), bin.no=rep(1:length(labels), reps), bin=rep(labels, reps))
    results <- merge(frame, data, by=c("rep.no", "bin"), all=TRUE)
    results[is.na(results)] <- 0
    
    #Calculate rates of change, if necessary (this is the slow bit, so default is to skip)
    if(RoC==TRUE) {
        for(i in 1:(nrow(results)-1)) {
            results[i,RoC:=(results[i+1,count]-results[i,count])/bin.width]
        }
        results[bin==labels[length(labels)], RoC:=NA]
    }
    
    #Return results
    results
}

# Define function to simulate a dummy set by sampling from within a specified distribution

dummy.simulate <- function(weight, probs=1, breaks=NULL, filter=NULL, start.date=0, end.date=2000, bin.width=100, reps=100, RoC=FALSE) {
    #Load required package
    require(data.table)
    
    #Tidy up input data and apply filters
    if(is.vector(weight)==1 & length(weight)==1) {weight <- rep(1, weight)} #if weight is a single value, use as number of entities
    dummy <- data.table(weight) #convert weights to data table format, if necessary.
    if(length(filter)>0 & "group"%in%colnames(dummy)) {dummy <- dummy[group==filter,]} #filter if appropriate
    
    #Set up breaks and labels
    if(is.null(breaks)==TRUE) {breaks <- seq(start.date, end.date, bin.width)} #if breaks not specified, sets them based on other arguments
    labels <- numeric(0)
    for(i in 1:(length(breaks)-1)) {
        labels[i] <- paste(breaks[i], breaks[i+1], sep="-") #sets bin labels based on breaks
    }
    probs <- data.table(cbind(probs, labels)) #append labels to relative probs, recycling the latter if necessary
    
    #Perform simulation
    rep.no <- rep(1:reps, each=nrow(dummy))
    dummy <- cbind(rep.no, dummy) #recycles input data 'reps' times to provide frame for simulation 
    p.sum <- sum(as.numeric(probs$probs)) #'stacks up' all the relative probabilities
    p.breaks <- c(0, cumsum(as.numeric(probs$probs))) #uses cumulative sum of relative probabilities to set breaks
    dummy[,sim:=runif(nrow(dummy), 0, p.sum)] #samples from within p.sum
    dummy[,bin:=cut(sim, p.breaks, labels=probs$labels)] #finds the relevant bin for each simulated date
    dummy <- dummy[is.na(bin)==FALSE, j=list(dummy=sum(as.numeric(weight))), by=list(rep.no,bin)] #sums weights by bin and rep number
    
    #Prepares and returns results table
    frame <- data.table("rep.no"=rep(1:reps, each=length(labels)), "bin.no"=rep(1:length(labels), reps), "bin"=rep(labels, reps))
    results <- merge(frame, dummy, by=c("rep.no", "bin"), all=TRUE)
    results[is.na(results)] <- 0
    
    #Calculate rates of change, if necessary (this is the slow bit, so default is to skip)
    if(RoC==TRUE) {
        for(i in 1:(nrow(results)-1)) {
            results[i,RoC:=(results[i+1,dummy]-results[i,dummy])/bin.width]
        }
        results[bin==labels[length(labels)], RoC:=NA]
    }
    
    #Return results
    results
}

# Define function that performs both 'real' and dummy simulation on target bone data

freq.simulate <- function(data, probs=1, weight=1, filter.field="group", filter.values=NULL, quant.list=c(0.025,0.25,0.5,0.75,0.975), start.date=0, end.date=2000, bin.width=100, reps=100, RoC=FALSE, save.full=FALSE, save.summ=FALSE) {
    #Load required packages
    require(data.table)
    require(reshape2)
    
    #Tidy up input data; apply filters
    data <- data.table(data)  #just in case it isn't already in this format
    if(is.null(filter.values)==FALSE) {  #if criteria have been provided...
        setnames(data, old=filter.field, new="FILTER")   #...sets the filter column...
        data <- data[FILTER %in% filter.values,]  #...and applies the criteria
    } else {filter.values <- "ALL"}
    data <- data[End >= start.date & Start <= end.date]  #drops records outside the date range FROM BOTH SIMULATION SETS
    if(length(weight)==1) {weight=rep(weight, nrow(data))} #if weight set as constant, repeats to length of data
    
    #Reset bin.width based on probs, if necessary
    if(is.vector(probs)==TRUE & length(probs)>1) {bin.width <- (end.date-start.date)/length(probs)}  #if probs supplied, use to set bin.widths
    if(sum(class(probs)=="data.frame")==1) {bin.width <- (end.date-start.date)/nrow(probs)}  #likewise if supplied as data.frame/data.table
   
    #Simulate from real data, then generate dummy set. Merge the two together.
    real <- date.simulate(data=data[,list(Start, End)], weight=weight, bin.width=bin.width, start.date=start.date, end.date=end.date, reps=reps)
    dummy <- dummy.simulate(weight=weight, probs=probs, breaks=breaks, start.date=start.date, end.date=end.date, bin.width=bin.width, reps=reps)    
    results <- merge(real, dummy, by=c("rep.no", "bin", "bin.no"), all=TRUE)
    
    #Calculate rate of change variables (could be done within core functions, but faster to loop through together here)
    if(RoC==TRUE) {
        for(i in 1:(nrow(results)-1)) {
            results[i,RoC.count:=(results[i+1,count]-results[i,count])/bin.width]
            results[i,RoC.dummy:=(results[i+1,dummy]-results[i,dummy])/bin.width]
        }
        results[bin==unique(bin)[length(unique(bin))], RoC.count:=NA]
        results[bin==unique(bin)[length(unique(bin))], RoC.dummy:=NA]
    }
        
    #Save full dataset, if requested
    if(save.full==TRUE) {write.csv(results, paste("FULL_", filter.values[1], "_simulated_by_period", params, ".csv",sep=""), row.names=FALSE)}

    #Create summary dataset
    real.summary <- sim.summ(results[,list(rep.no, bin, bin.no, count)], quant.list)
    dummy.summary <- sim.summ(results[,list(rep.no, bin, bin.no, dummy)], quant.list)  
    summary <- merge(real.summary, dummy.summary, by="bin")
   
    #Summarise and append RoC results, if requested
    if(RoC==TRUE) {
        RoC.count.summary <- sim.summ(results[,list(rep.no, bin, bin.no, RoC.count)], quant.list)
        RoC.dummy.summary <- sim.summ(results[,list(rep.no, bin, bin.no, RoC.dummy)], quant.list)
        RoC.summary <- merge(RoC.count.summary, RoC.dummy.summary, by="bin")
        summary <- merge(summary, RoC.summary, by="bin")
    }
    
    #Save summary dataset, if requested
    if(save.summ==TRUE) {write.csv(summary, paste("SUMMARY_", filter.values[1], "_simulated_by_period", params, ".csv", sep=""), row.names=FALSE)}

    #Return list with full and summary datasets
    list(results, summary)
}

#Define function to create summary table from results of data.simulate or dummy.simulate

sim.summ <- function(results, quant.list=c(0.025,0.25,0.5,0.75,0.975)) {
    #Load required packages
    require(data.table)
    require(reshape2)
    
    #Save name of data field then replace with "data" (easier than trying to refer to it by number below)
    data.field <- colnames(results)[4]
    setnames(results, 4, "count")    

    #Create and return summary table
    summary <- results[,quantile(count, probs=quant.list, na.rm=TRUE), by=bin]
    summary[,id:=paste(rep(data.field, length(quant.list)), quant.list, sep="_")]
    dcast.data.table(summary, bin ~ id, value.var="V1") 
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

