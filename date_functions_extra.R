# Variation of freq.simulate. Renamed freq.simulate.species so not to be confused with the elements freq.simulate. Boxplot graphs included for read and dummy data

freq.simulate.species <- function(data, probs, filter.field="Species", filter.values=NULL, quant.list=c(0.025,0.25,0.5,0.75,0.975), start.date=0, end.date=2000, rep=1000) {
    require(data.table)
    require(reshape2)
    data <- data.table(data)  #just in case it isn't already in this format
    if(is.null(filter.values)==FALSE) {  #if criteria have been provided...
        setnames(data, old=filter.field, new="FILTER")   #...sets the filter column...
        data <- data[FILTER %in% filter.values,]  #...and applies the criteria
    } else {filter.values <- "ALL"}
    data <- data[End >= start.date & Start <= end.date]  #drops records outside the date range FROM BOTH SIMULATION SETS
    bin.width <- (end.date-start.date)/nrow(probs)  #sets bin widths (and hence no. of bins) to match the calibration dataset
    
    # simulate from real data
    real <- date.simulate(data=data[,list(Start, End)], weight=data[,Frag], bin.width=bin.width, start.date=start.date, end.date=end.date, rep=rep)
    setnames(real, old="V1", new="real")
    boxplot(real ~ bin, data=real, main=paste(filter.values, "Frequency simulate"), xlab="Time (50yr intervals)", ylab="Frequency")
    
    
    # simulate dummy set
    dummy <- dummy.simulate(probs=probs, weight=data[,Frag], start.date=start.date, end.date=end.date, rep=rep)    
    setnames(dummy, old="V1", new="dummy")
    boxplot(dummy ~ bin, data=dummy, main=paste(filter.values,"Frequency dummy"), xlab="Time (50yr intervals)", ylab="Frequency")
    
    # set up list of all bins and rep no.s
    labels <- numeric(nrow(probs))
    for(i in 1:length(labels)) {
        labels[i] <- paste(breaks[i], breaks[i+1], sep="-")
    } 
    frame <- data.table(rep(1:rep, each=length(labels)), rep(labels, rep))
    setnames(frame, old=c("V1", "V2"), new=c("rep.no", "bin"))
    
    # merge the above three data.tables together
    results <- merge(frame, dummy, by=c("rep.no", "bin"), all=TRUE)
    results <- merge(results, real, by=c("rep.no", "bin"), all=TRUE)
    results[is.na(real)==TRUE, real:=0] 
    results[is.na(dummy)==TRUE, dummy:=0]
    
    # save full dataset
    write.csv(results, paste("TEST_", filter.values[1], "_simulated_by_period", params, ".csv",sep=""), row.names=FALSE)
    
    # create summary dataset
    real.summary <- results[,quantile(real, probs=quant.list), by=bin]
    real.summary[,id:=paste(rep("real", length(quant.list)), quant.list, sep="_")]
    real.summary <- dcast.data.table(real.summary, bin ~ id, value.var="V1")
    dummy.summary <- results[,quantile(dummy, probs=quant.list), by=bin]
    dummy.summary[,id:=paste(rep("dummy", length(quant.list)), quant.list, sep="_")]
    dummy.summary <- dcast.data.table(dummy.summary, bin ~ id, value.var="V1")
    summary <- cbind(real.summary, dummy.summary)
    
    #save summary dataset
    write.csv(summary, paste("TEST_summary_", filter.values[1], "_simulated_by_period", params, ".csv", sep=""), row.names=FALSE)
    summary
}

# A variation of freq.simulate used for filtering by Elements.. Renamed freq.simulate.element so not to be confused with the species freq.simulate. Boxplot graphs included for read and dummy data.
# Further work needed on boxplot and CSV files to include speices name as well as the filter.values. 
freq.simulate.element <- function(data, probs, filter.field="Elements", filter.values=NULL, quant.list=c(0.025,0.25,0.5,0.75,0.975), start.date=0, end.date=2000, rep=1000) {
    require(data.table)
    require(reshape2)
    data <- data.table(data)  #just in case it isn't already in this format
    if(is.null(filter.values)==FALSE) {  #if criteria have been provided...
        setnames(data, old=filter.field, new="FILTER")   #...sets the filter column...
        data <- data[FILTER %in% filter.values,]  #...and applies the criteria
    } else {filter.values <- "ALL"}
    data <- data[End >= start.date & Start <= end.date]  #drops records outside the date range FROM BOTH SIMULATION SETS
    bin.width <- (end.date-start.date)/nrow(probs)  #sets bin widths (and hence no. of bins) to match the calibration dataset
    
    # simulate from real data
    real <- date.simulate(data=data[,list(Start, End)], weight=data[,Frag], bin.width=bin.width, start.date=start.date, end.date=end.date, rep=rep)
    setnames(real, old="V1", new="real")
    boxplot(real ~ bin, data=real, main=paste(filter.values[1], "Frequency simulate"), xlab="Time (50yr intervals)", ylab="Frequency")
    
    # simulate dummy set
    dummy <- dummy.simulate(probs=probs, weight=data[,Frag], start.date=start.date, end.date=end.date, rep=rep)    
    setnames(dummy, old="V1", new="dummy")
    boxplot(dummy ~ bin, data=dummy, main=paste(filter.values[1],"Frequency dummy"), xlab="Time (50yr intervals)", ylab="Frequency")
    
    # set up list of all bins and rep no.s
    labels <- numeric(nrow(probs))
    for(i in 1:length(labels)) {
        labels[i] <- paste(breaks[i], breaks[i+1], sep="-")
    } 
    frame <- data.table(rep(1:rep, each=length(labels)), rep(labels, rep))
    setnames(frame, old=c("V1", "V2"), new=c("rep.no", "bin"))
    
    # merge the above three data.tables together
    results <- merge(frame, dummy, by=c("rep.no", "bin"), all=TRUE)
    results <- merge(results, real, by=c("rep.no", "bin"), all=TRUE)
    results[is.na(real)==TRUE, real:=0] 
    results[is.na(dummy)==TRUE, dummy:=0]
    
    # save full dataset
    write.csv(results, paste(filter.values[1], "_simulated_by_period", params, ".csv",sep=""), row.names=FALSE)
    
    # create summary dataset
    real.summary <- results[,quantile(real, probs=quant.list), by=bin]
    real.summary[,id:=paste(rep("real", length(quant.list)), quant.list, sep="_")]
    real.summary <- dcast.data.table(real.summary, bin ~ id, value.var="V1")
    dummy.summary <- results[,quantile(dummy, probs=quant.list), by=bin]
    dummy.summary[,id:=paste(rep("dummy", length(quant.list)), quant.list, sep="_")]
    dummy.summary <- dcast.data.table(dummy.summary, bin ~ id, value.var="V1")
    summary <- cbind(real.summary, dummy.summary)
    
    #save summary dataset
    write.csv(summary, paste("Species_element_summary_", filter.values[1], "_simulated_by_period", params, ".csv", sep=""), row.names=FALSE)
    summary
}
