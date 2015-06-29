# Define function to calculate aoristic sum

aorist <- function(data, start.date=0, end.date=2000, bin.width=100, weight=1) { 
    #Load required package
    require(data.table)
    
    #Tidies up input data
    data <- cbind(data, weight) #appends weights to list of date ranges, recycling if necessary (e.g. for uniform weight)
    data <- data[End >= start.date & Start <= end.date] #excludes ranges that fall entirely outside the study period

    #Set up columns for duration and for weight per year
    data[,duration:=End-Start]
    data[,weight.per.unit:=weight/duration]
    
    #Set up breaks and labels
    breaks <<- seq(start.date, end.date, bin.width) #creates and saves vector of breaks
    labels <- numeric(0)
    for(i in 1:(length(breaks)-1)) {
        labels[i] <- paste(breaks[i], breaks[i+1], sep="-") #sets bin labels
    }
    params <<- paste("_", start.date, "-", end.date, "_by_", bin.width, sep="") #saves character value with key parameters 
    
    #Set frame for results
    aorist <- data.table(bin=labels, bin.no=1:length(labels), aorist=0)
    
    #Cycle through bins, assigning probability mass to cases where appropriate
    for(i in 1:((end.date-start.date)/bin.width)) {
        bin.2 <- start.date + i*bin.width #Find end date of bin
        bin.1 <- bin.2 - bin.width #Find start date of bin
        data[Start>=bin.1 & Start<bin.2, assign("a", labels[i]):=(bin.2-Start)*weight.per.unit] #
        data[End>bin.1 & End<=bin.2, assign("a", labels[i]):=(End-bin.1)*weight.per.unit]
        data[Start<bin.1 & End>bin.2, assign("a", labels[i]):=bin.width*weight.per.unit]
        data[Start>=bin.1 & End<=bin.2, assign("a", labels[i]):=as.double(weight)]
        aorist$aorist[i] <- sum(data[,get(labels[i])], na.rm=TRUE)
    }
    aorist
}

# Define function to simulate distribution of dates

date.simulate <- function(data, weight=1, UoA=NULL, context.fields=c("SITE_C"), start.date=0, end.date=2000, bin.width=100,
                          reps=100, RoC=FALSE, summ=TRUE, ...) {
    #Load required package
    require(data.table) 
    
    #Tidy up input data and apply filters
    data <- data.table(cbind(data, weight)) #appends weights to list of date ranges, recycling if necessary (e.g. for uniform weight) 
    data <- data[End >= start.date & Start <= end.date] #excludes ranges that fall entirely outside the study period
        
    #Aggregate data
    agg.list <- c(context.fields, "Start", "End", UoA)
    data <- data[, j=list(weight=sum(as.numeric(weight))), by=agg.list]
        
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
    
    #Calculate rates of change, if necessary
    if(RoC==TRUE) {
        results[, RoC:=c(diff(count), NA)]
        results[bin==labels[length(labels)], RoC:=NA]
    }
    
    #Create summary dataset if required, and return results
    if(summ==TRUE) {
        summary <- sim.summ(results)
        results <- list(results, summary)
    }
    results
}

# Define function to simulate a dummy set by sampling from within a specified distribution

dummy.simulate <- function(weight, probs=1, breaks=NULL, start.date=0, end.date=2000, bin.width=100, reps=100, RoC=FALSE, summ=TRUE, ...) {
    #Load required package
    require(data.table)
    
    #Tidy up input data and apply filters
    if(is.vector(weight)==1 & length(weight)==1) {weight <- rep(1, weight)} #if weight is a single value, use as number of entities
    dummy <- data.table(weight) #convert weights to data table format, if necessary.
    
    #Set up breaks and labels
    if(is.null(breaks)==TRUE) {breaks <- seq(start.date, end.date, bin.width)} #if breaks not specified, sets them based on other arguments
    labels <- numeric(0)
    for(i in 1:(length(breaks)-1)) {
        labels[i] <- paste(breaks[i], breaks[i+1], sep="-") #sets bin labels based on breaks
    }
    probs <- cbind(1:length(labels), probs)
    
    #Perform simulation
    rep.no <- rep(1:reps, each=nrow(dummy))
    dummy <- cbind(rep.no, dummy) #recycles input data 'reps' times to provide frame for simulation 
    dummy[,bin:=sample(labels, size=nrow(dummy), replace=TRUE, prob=probs[,2])]
    dummy <- dummy[, j=list(dummy=sum(as.numeric(weight))), by=list(rep.no,bin)] #sums weights by bin and rep number
    
    #Prepares and returns results table
    frame <- data.table("rep.no"=rep(1:reps, each=length(labels)), "bin.no"=rep(1:length(labels), reps), "bin"=rep(labels, reps))
    results <- merge(frame, dummy, by=c("rep.no", "bin"), all=TRUE)
    results[is.na(dummy), dummy:=0]
    results <- results[order(rep.no, bin.no)]
    
    #Calculate rates of change, if necessary
    if(RoC==TRUE) {
        results[, RoC:=c(diff(dummy), NA)]
        results[bin==labels[length(labels)], RoC:=NA]
    }
    
    #Create summary dataset if required, and return results
    if(summ==TRUE) {
        summary <- sim.summ(results)
        results <- list(results, summary)
    }
    results
}

# Define function that performs both 'real' and dummy simulation on target bone data

freq.simulate <- function(data, probs=1, weight=1, UoA=NULL, context.fields=c("SITE_C"), quant.list=c(0.025,0.25,0.5,0.75,0.975),
                          start.date=0, end.date=2000, bin.width=100, reps=100, RoC=FALSE, ...) {
    #Load required packages
    require(data.table)

    #Tidy up input data; apply filters
    data <- data.table(cbind(data, weight)) #appends weights to list of date ranges, recycling if necessary (e.g. for uniform weight) 
    data <- data[End >= start.date & Start <= end.date]  #drops records outside the date range FROM BOTH SIMULATION SETS
    
    #Aggregate data
    agg.list <- c(context.fields, "Start", "End", UoA)
    data <- data[, j=list(weight=sum(as.numeric(weight))), by=agg.list]
    if(length(weight)==1) {data[,weight:=1]}
    
    #Reset bin.width based on probs, if necessary
    if(is.vector(probs)==TRUE & length(probs)>1) {bin.width <- (end.date-start.date)/length(probs)}  #if probs supplied, use to set bin.widths
    if(sum(class(probs)=="data.frame")==1) {bin.width <- (end.date-start.date)/nrow(probs)}  #likewise if supplied as data.frame/data.table
   
    #Simulate from real data, then generate dummy set. Merge the two together.
    real <- date.simulate(data=data, weight=data$weight, context.fields=context.fields, UoA=UoA, bin.width=bin.width, start.date=start.date, end.date=end.date, reps=reps, summ=FALSE)
    dummy <- dummy.simulate(weight=data$weight, probs=probs, breaks=breaks, start.date=start.date, end.date=end.date, bin.width=bin.width, reps=reps, summ=FALSE)    
    results <- merge(real, dummy, by=c("rep.no", "bin.no", "bin"), all=TRUE)
    
    #Calculate rate of change variables (could be done within core functions, but faster to loop through together here)
    if(RoC==TRUE) {
        results[, RoC.count:=c(diff(count), NA)]
        results[, RoC.dummy:=c(diff(dummy), NA)]
        results[bin==unique(bin)[length(unique(bin))], RoC.count:=NA]
        results[bin==unique(bin)[length(unique(bin))], RoC.dummy:=NA]
    }
           
    #Create summary dataset and return results
    summary <- sim.summ(results)
    list(results, summary)
}

#Define function to create summary table from results of data.simulate or dummy.simulate

sim.summ <- function(results, summ.col=NULL, quant.list=c(0.025,0.25,0.5,0.75,0.975)) {
    #Load required packages
    require(data.table)

    if(is.null(summ.col)==TRUE) {summ.col <- colnames(results)[!colnames(results)%in%c("rep.no", "bin", "bin.no")]}
    
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

#Function for comparisons
comp.simulate <- function(data, probs=1, weight=1, comp.values=NULL, comp.field="group", context.fields=c("SITE_C", "SITE_S"), UoA=NULL, 
                    comp.fun=date.simulate, quant.list=c(0.025,0.25,0.5,0.75,0.975), start.date=0, end.date=2000, bin.width=100, reps=100, RoC=FALSE) {
    #Load required packages
    require(data.table)
        
    #Deal with comparison field
    data <- data.table(cbind(data, weight)) #appends weights to list of date ranges, recycling if necessary (e.g. for uniform weight) 
    if(is.null(comp.values)==TRUE) {comp.values <- unique(data[,get(comp.field)])} #compare all values of comp.field if not specified otherwise
    
    #Aggregate data
    if(is.null(UoA)==FALSE) {
        if(UoA==comp.field) {UoA <- NULL}
    }
    agg.list <- c(context.fields, "Start", "End", comp.field, UoA)
    data <- data[, j=list(weight=sum(as.numeric(weight))), by=agg.list]
#    form <- as.formula(paste(paste(context.fields, collapse='+'), "+Start+End ~ ", comp.field))
#    data <- dcast.data.table(data=data, formula=form, value.var="weight",fun.aggregate=sum)
    
    #Run the requested function for each group
    for(i in 1:length(comp.values)) {
        y <- comp.fun(data=data[get(comp.field)==comp.values[i],], probs=probs, weight=data[get(comp.field)==comp.values[i],weight], filter.field=filter.field, filter.values=filter.values, quant.list=quant.list,
                 UoA=UoA, context.field=context.fields, start.date=start.date, end.date=end.date, bin.width=bin.width, reps=reps, RoC=RoC, summ=FALSE)
        if(i==1) {results <- y[,list(rep.no, bin, bin.no)]}
        setnames(y, old=4:ncol(y), new=paste(comp.values[i], colnames(y)[4:ncol(y)], sep=".")) 
        results <- merge(results, y, by=c("rep.no", "bin.no", "bin"))
    }
    
    #Create summary dataset and return results
    summary <- sim.summ(results)
    list(results, summary)
}

##CPUE function

cpue <- function(x, y, weight.x=1, weight.y=1, context.fields=c("SITE_C"), UoA=NULL, quant.list=c(0.025,0.25,0.5,0.75,0.975), start.date=0,
                 end.date=2000, bin.width=100, reps=100, RoC=FALSE, small.n=NULL, ...) {
    #Load required package
    require(data.table)
    
    #Tidy up input data and apply filters
    x <- data.table(cbind(x, weight.x)) #appends weights to list of date ranges, recycling if necessary (e.g. for uniform weight) 
    x <- x[End >= start.date & Start <= end.date] #excludes ranges that fall entirely outside the study period
    y <- data.table(cbind(y, weight.y)) #appends weights to list of date ranges, recycling if necessary (e.g. for uniform weight) 
    y <- y[End >= start.date & Start <= end.date] #excludes ranges that fall entirely outside the study period
    
    #Aggregate data
    x <- x[, j=list(weight.x=sum(as.numeric(weight.x))), by=c(context.fields, "Start", "End", UoA)]
    y <- y[, j=list(weight.y=sum(as.numeric(weight.y))), by=c(context.fields, "Start", "End")]
    if(length(weight.x)==1) {x[,weight.x:=1]}
    if(length(weight.y)==1) {y[,weight.y:=1]}
    data <- merge(x,y,by=c(context.fields, "Start", "End"), all=TRUE)
    
    #Set up breaks and labels
    breaks <<- seq(start.date, end.date, bin.width) #sets breaks and saves them externally
    labels <- numeric(0)
    for(i in 1:(length(breaks)-1)) {
        labels[i] <- paste(breaks[i], breaks[i+1], sep="-") #sets bin labels
    }
    
    #Perform simulation
    rep.no <- rep(1:reps, each=nrow(data))
    data <- cbind(rep.no, data) #recycles input data 'reps' times to provide frame for simulation 
    data[,sim:={x<-runif(nrow(data)); (x*(data[,End]-data[,Start]))+data[,Start]}] #simulates a date for each row
    data[,bin:=cut(sim,breaks,labels=labels)] #finds the relevant bin for each simulated date
    
    #Prepare results table
    x <- data[is.na(weight.x)==FALSE&is.na(bin)==FALSE, j=list(x=sum(as.numeric(weight.x))), by=list(rep.no,bin)] #sums weights by bin and rep number
    y <- data[is.na(weight.y)==FALSE&is.na(bin)==FALSE, j=list(y=sum(as.numeric(weight.y))), by=list(rep.no,bin)] #sums weights by bin and rep number
    frame <- data.table(rep.no=rep(1:reps, each=length(labels)), bin.no=rep(1:length(labels), reps), bin=rep(labels, reps))
    results <- merge(frame, x=x, by=c("rep.no", "bin"), all=TRUE)
    results <- merge(results, y=y, by=c("rep.no", "bin"), all=TRUE)
    results[,cpue:=x/y]
    results[is.na(results)] <- 0
    
    #Calculate rates of change, if necessary
    if(RoC==TRUE) {
        results[, RoC:=c(diff(cpue), NA)]
        results[bin==labels[length(labels)], RoC:=NA]
    }
    
    #Create summary dataset
    summary <- sim.summ(results, summ.col="cpue", quant.list=quant.list)
    
    #Define small-n polygons if required
    if(length(small.n > 0)) {
        n <- results[,median(y), by="bin"]
        ylim <- max(results$cpue*1.1)
        boxes <- list()
        for(i in 1:length(small.n)) {
            weak <- c(FALSE, n$V1 < small.n[i], FALSE)
            weak.x <- rep((0:length(labels))[diff(weak)==1|diff(weak)==-1], each=2)+0.5
            weak.y <- rep(c(-1,ylim, ylim, -1), length(weak.x)/4)
            coords <- data.frame(cbind(weak.x, weak.y))
            coords <- split(coords, rep(1:(nrow(coords)/4), each=4))
            boxes[[i]] <- coords
        }
    }

    #Return results
    results <- list(full=results, summary=summary)
    if(length(small.n)>0) {results$small.n <- boxes}
    results
    
}

# Function to convert results to survivorship format
surv.convert <- function(results, field.list=NULL, summ=TRUE) {
    #Load required package
    require(data.table)
    
    #Select only full results (if applicable); establish field names to use
    if(class(results)[1]=="list") {results <- results[[1]]}
    if(is.null(field.list)==TRUE) {field.list <- colnames(results)[!colnames(results)%in%c("bin","bin.no", "rep.no")]}
    
    #Build new data table for results
    frame <- data.table(rep.no=rep(1:max(results$rep.no), each=(max(results$bin.no)+1)), bin.no=rep(1:(max(results$bin.no)+1)), bin=rep(c("Start", unique(as.character(results$bin))), max(results$rep.no)))
    column.names <- paste("survive.", field.list, sep="")
    
    #Calculate survivorship for each column, run, and bin (nested in that order)
    for(k in 1:length(field.list)) {
        n <- results[rep.no==1&!is.na(get(field.list[k])), sum(get(field.list[k]))] #Find total sample size for current column
        frame[bin.no==1,assign("a", column.names[k]):=n] #Create column and fill in start value for each run
        frame[!bin.no==1, assign("a", column.names[k]):=results[,get(field.list[k])]]    
        for(i in 1:max(frame$rep.no)) {
            frame[rep.no==i,assign("a", column.names[k]):=results[rep.no==i, n-diffinv(get(field.list[k]))]]
        }
        frame[,assign("a", column.names[k]):=get(column.names[k])/n]
    }
    
    #Create summary dataset if required
    if(summ==TRUE) {
        summary <- sim.summ(frame)
        frame <- list(frame, summary)
    }

    #Return results
    frame
}


