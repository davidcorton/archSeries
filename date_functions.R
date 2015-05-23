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

date.simulate <- function(data, weight=1, filter=NULL, start.date=0, end.date=2000, bin.width=100, reps=100, RoC=FALSE, summ=TRUE, ...) {
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
    
    #Create summary dataset if required
    if(summ==TRUE) {
        summary <- sim.summ(results)
        results <- list(results, summary)
    }
    
    #Return results
    results
}

# Define function to simulate a dummy set by sampling from within a specified distribution

dummy.simulate <- function(weight, probs=1, breaks=NULL, filter.values=NULL, start.date=0, end.date=2000, bin.width=100, reps=100, RoC=FALSE, summ=TRUE, ...) {
    #Load required package
    require(data.table)
    
    #Tidy up input data and apply filters
    if(is.vector(weight)==1 & length(weight)==1) {weight <- rep(1, weight)} #if weight is a single value, use as number of entities
    dummy <- data.table(weight) #convert weights to data table format, if necessary.
    if(length(filter.values)>0 & "group"%in%colnames(dummy)) {dummy <- dummy[group%in%filter.values,]} #filter if appropriate
    
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
    
    #Create summary dataset if required
    if(summ==TRUE) {
        summary <- sim.summ(results)
        results <- list(results, summary)
    }
    
    #Return results
    results
}

# Define function that performs both 'real' and dummy simulation on target bone data

freq.simulate <- function(data, probs=1, weight=1, filter.field="group", filter.values=NULL, quant.list=c(0.025,0.25,0.5,0.75,0.975), start.date=0, end.date=2000, bin.width=100, reps=100, RoC=FALSE, summ=TRUE, ...) {
    #Load required packages
    require(data.table)
    require(reshape2)
    
    #Tidy up input data; apply filters
    data <- data.table(data)  #just in case it isn't already in this format
    if(is.null(filter.values)==FALSE) {  #if criteria have been provided...
        data <- data[get(filter.field) %in% filter.values,]  #...applies the criteria
    } else {filter.values <- "ALL"} #this is just for the sake of setting 'params' in date.simulate
    data <- data[End >= start.date & Start <= end.date]  #drops records outside the date range FROM BOTH SIMULATION SETS
    if(length(weight)==1 & is.vector(weight)==TRUE) {weight=rep(weight, nrow(data))} #if weight set as constant, repeats to length of data
    
    #Reset bin.width based on probs, if necessary
    if(is.vector(probs)==TRUE & length(probs)>1) {bin.width <- (end.date-start.date)/length(probs)}  #if probs supplied, use to set bin.widths
    if(sum(class(probs)=="data.frame")==1) {bin.width <- (end.date-start.date)/nrow(probs)}  #likewise if supplied as data.frame/data.table
   
    #Simulate from real data, then generate dummy set. Merge the two together.
    real <- date.simulate(data=data[,list(Start, End)], weight=weight, bin.width=bin.width, start.date=start.date, end.date=end.date, reps=reps, summ=FALSE)
    dummy <- dummy.simulate(weight=weight, probs=probs, breaks=breaks, start.date=start.date, end.date=end.date, bin.width=bin.width, reps=reps, summ=FALSE)    
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
           
    #Create summary dataset if required
    if(summ==TRUE) {
        summary <- sim.summ(results)
        results <- list(results, summary)
    }
    
    #Return results
    results
}

#Define function to create summary table from results of data.simulate or dummy.simulate

sim.summ <- function(results, summ.col=NULL, quant.list=c(0.025,0.25,0.5,0.75,0.975)) {
    #Load required packages
    require(data.table)
    require(reshape2)
    
    if(is.null(summ.col)==TRUE) {summ.col <- colnames(results)[4:ncol(results)]}
    
    #Create summary tables
    for(i in 1:length(summ.col)) {    
        x <- results[,quantile(get(summ.col[i]), probs=quant.list, na.rm=TRUE), by=bin] #calculate quantiles
        x[,quantile:=quant.list] #create column to specify quantiles        
        x[,id:=summ.col[i]] #create column to specify variable       
        if(i==1) {summary <- data.table()}
        summary <- rbind(summary, x)
    }   
    
    #Return summary table
    summary
}

#Function for comparisons - WORK IN PROGRESS
comp.simulate <- function(data, probs=1, weight=1, comp.values=NULL, comp.field="group", comp.fun=date.simulate, filter.field="group", filter.values=NULL, quant.list=c(0.025,0.25,0.5,0.75,0.975), start.date=0, end.date=2000, bin.width=100, reps=100, RoC=FALSE, summ=TRUE) {
    #Load required packages
    require(data.table)
    
    #Deal with comparison and filter fields
    if(is.null(comp.values)==TRUE) {comp.values <- unique(data[,get(comp.field)])} #compare all values of comp.field if not specified otherwise
    if(filter.field==comp.field) {comp.values <- comp.values[comp.values%in%filter.values]} #removes filtered values from comp.values
    
    #Run the requested function for each group
    if(length(weight)==1 & is.vector(weight)==TRUE) {weight=rep(weight, nrow(data))} #if weight set as constant, repeats to length of data
    for(i in 1:length(comp.values)) {
        y <- comp.fun(data=data[get(comp.field)==comp.values[i],], probs=probs, weight=weight[data[,get(comp.field)]==comp.values[i]], filter.field=filter.field, filter.values=filter.values, quant.list=quant.list,
                 start.date=start.date, end.date=end.date, bin.width=bin.width, reps=reps, RoC=RoC, summ=FALSE)
        if(class(y)[1]=="list") {y <- y[[1]]} #takes only full results for freq.simulate
        if(i==1) {results <- y[,list(rep.no, bin, bin.no)]}
        setnames(y, old=4:ncol(y), new=paste(comp.values[i], colnames(y)[4:ncol(y)], sep=".")) 
        results <- merge(results, y, by=c("rep.no", "bin", "bin.no"))
    }
    
    #Create summary dataset if required
    if(summ==TRUE) {
        summary <- sim.summ(results)
        results <- list(results, summary)
    }
    
    #Return results
    results
}


