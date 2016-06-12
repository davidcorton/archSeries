#' Catch Per Unit Effort simulator
#' 
#' Estimates the relationship between a denominator and a numerator variable over time 
#' 
#' @param catch Data table (or object that can be coerced to one) with, minimally, two numeric columns called Start and End.
#' @param effort Data table (or object that can be coerced to one) with, minimally, two numeric columns called Start and End.
#' @param wt.catch Numeric vector: the weight to be applied to each row in `catch', or a constant weight to be applied to all. 
#'      Defaults to 1.
#' @param wt.effort Numeric vector: the weight to be applied to each row in `effort`, or a constant weight to be applied to all.
#'      Defaults to 1.
#' @param context.fields Character vector specifying the column(s) in data which define the minimal stratigraphic entities to analyse. 
#'      Defaults to "SITE_C".
#' @param UoA Unit of Analysis: character vector of names of additional columns by which to group data when aggregating weights,
#'      on top of those specified in context.field. For example, should different taxa be lumped together when analysing bone remains 
#'      from a table of contexts? Defaults to NULL.
#' @param quant.list Numeric vector of quantiles to be calculated in a summary table. Defaults to c(0.025,0.25,0.5,0.75,0.975).
#' @param start.date Numeric: the start of time period to be considered. Defaults to 0.
#' @param end.date Numeric: the end of time period to be considered. Defaults to 2000.
#' @param bin.width Numeric: the resolution of the analysis, in units of time. Defaults to 100.
#' @param reps Integer: the number of times the simulation will be run. Defaults to 100.
#' @param RoC Rate of Change. Logical: should rates of change between adjacent bins be calculated alongside the raw counts?
#' @param small.n Numeric: vector specifying one or more cut-off values to be used to flag periods of where the underlying sample size (i.e. 
#'      magnitude of effort) is small. If multiple values are passed they should be in descending order.
#' @return Minimally, a list with two named elements:
#'      "full" is a data table with six or seven columns: 'rep.no', integer specifying simulation run; 'bin', character specifying 
#'      chronological bin in terms of date range; 'bin.no' integer specifying number of bin, counting from earliest; 'catch', giving the 
#'      simulated frequency of 'catch'; and 'effort', giving the simulated magnitude of 'effort'. If RoC=TRUE there will be an additional 
#'      column, 'RoC', giving the rate of change in cpue between the given bin and the next.
#'      "summary" is a second long format data table with four named columns: 'bin', as above; 'V1', the relevant value for the given bin 
#'      at a given quantile; 'quantile', the quantile at which V1 is calculated; 'id', character specifying which column from "full" V1 is 
#'      based on: "cpue" or "RoC" (catch and effort are ignored when summarising).
#'      If one or more values are passed to 'small.n', then this list will have a third element, "boxes". This is a list of length equal 
#'      to the length of small.n. Each component is a list of four-row data frames giving coordinates that define boxes around periods 
#'      where the simulated value of effort is below the corresponding value in small.n, set up for plotting using grey.zones (either 
#'      alone or within any of the main archSeries plotting functions).
#' @export
#' @examples
#' dates <- data.table(unit=c(1, 2, 3, 4), Start=c(450, 450, 600, 1000), End=c(700, 800, 650, 1200), frags=c(3, 6, 2, 1), vol=c(40, 40, 40, 40))
#' x <- cpue(dates, dates, dates$frags, dates$vol, context.fields=NULL, small.n=1, reps=1000)

cpue <- function(catch, effort, wt.catch=1, wt.effort=1, context.fields=c("SITE_C"), UoA=NULL, quant.list=c(0.025, 0.25, 0.5, 0.75, 0.975), 
                 start.date=0, end.date=2000, bin.width=100, reps=100, RoC=FALSE, small.n=NULL, ...) {
    #Load required package
    require(data.table)
    
    #Tidy up input data and apply filters
    catch <- data.table(cbind(catch, wt.catch)) #appends weights to list of date ranges, recycling if necessary (e.g. for uniform weight) 
    catch <- catch[End >= start.date & Start <= end.date] #excludes ranges that fall entirely outside the study period
    effort <- data.table(cbind(effort, wt.effort)) #appends weights to list of date ranges, recycling if necessary (e.g. for uniform weight) 
    effort <- effort[End >= start.date & Start <= end.date] #excludes ranges that fall entirely outside the study period
    
    #Aggregate data
    catch <- catch[, j=list(wt.catch=sum(as.numeric(wt.catch))), by=c(context.fields, "Start", "End", UoA)]
    effort <- effort[, j=list(wt.effort=sum(as.numeric(wt.effort))), by=c(context.fields, "Start", "End")]
    if(length(wt.catch)==1) {catch[, wt.catch := 1]}
    if(length(wt.effort)==1) {effort[, wt.effort := 1]}
    data <- merge(catch, effort, by=c(context.fields, "Start", "End"), all=TRUE)
    
    #Set up breaks and labels
    breaks <<- seq(start.date, end.date, bin.width) #sets breaks and saves them externally
    labels <- numeric(0)
    for(i in 1:(length(breaks)-1 )) {
        labels[i] <- paste(breaks[i], breaks[i + 1], sep="-") #sets bin labels
    }
    
    #Perform simulation
    rep.no <- rep(1:reps, each=nrow(data))
    data <- cbind(rep.no, data) #recycles input data 'reps' times to provide frame for simulation 
    data[, sim := {x <- runif(nrow(data)); (x * (data[, End] - data[, Start])) + data[, Start]}] #simulates a date for each row
    data[, bin := cut(sim, breaks, labels=labels)] #finds the relevant bin for each simulated date
    
    #Prepare results table
    catch <- data[is.na(wt.catch)==FALSE & is.na(bin)==FALSE, j=list(catch=sum(as.numeric(wt.catch))), by=list(rep.no, bin)] #sums weights by bin and rep number
    effort <- data[is.na(wt.effort)==FALSE & is.na(bin)==FALSE, j=list(effort=sum(as.numeric(wt.effort))), by=list(rep.no, bin)] #sums weights by bin and rep number
    frame <- data.table(rep.no=rep(1:reps, each=length(labels)), bin.no=rep(1:length(labels), reps), bin=rep(labels, reps))
    results <- merge(frame, catch, by=c("rep.no", "bin"), all=TRUE)
    results <- merge(results, effort, by=c("rep.no", "bin"), all=TRUE)
    results[, cpue := catch / effort]
    results[is.na(results)] <- 0
    results <- results[order(rep.no, bin.no)]
    
    #Calculate rates of change, if necessary
    if(RoC==TRUE) {
        results[, RoC := c(diff(cpue), NA)]
        results[bin==labels[length(labels)], RoC := NA]
        RoC <- "RoC"
    } else {RoC <- NULL}
    
    #Create summary dataset
    summary <- sim.summ(results, summ.col=c("cpue", RoC), quant.list=quant.list)
    
    #Define small-n polygons if required
    if(length(small.n > 0)) {
        n <- results[,median(effort), by="bin"]
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