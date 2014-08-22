# Define function to calculate aoristic sum
# Arguments: 'data' is a data table with two numeric columns called Start and End
# 'start.date' is a single numeric value for the start of the time period to be analysed
# 'end.date' is a single numric value for the end end of the time period to be analysed
# 'bin.width' is a single numeric value setting the resolution of the analysis, in years
# 'weight' is a numeric vector giving a weight for each context/entity

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
    aoristic.sum
}

# Define function to simulate distribution of dates
# Arguments: 'data' is a data table with two numeric columns called Start and End
# 'start.date' is a single numeric value for the start of the time period to be analysed
# 'end.date' is a single numric value for the end end of the time period to be analysed
# 'bin.width' is a single numeric value setting the resolution of the analysis, in years
# 'rep' is the number of times the simulation will be run
# 'weight' is a numeric vector giving a weight for each context/entity

date.simulate <- function(data, start.date=0, end.date=2000, bin.width=100, rep=100, weight=1) {
    require(data.table)
    data <- cbind(data, weight)
    data <- data[End >= start.date & Start <= end.date]
    breaks <<- seq(start.date, end.date, bin.width)
    labels <- 1:(length(breaks)-1)
    data <- cbind(rep(1:rep, each=nrow(data)), data)
    data[,bin:={x<-runif(nrow(data)); (x*(data[,End]-data[,Start]))+data[,Start]}]
    data[,bin:=cut(bin,breaks,labels=labels)]
    data <- data[is.na(bin)==FALSE, sum(weight), by=list(V1,bin)]
    setnames(data, old=1, new="rep.no")
    data[order(rep.no, bin)]
}