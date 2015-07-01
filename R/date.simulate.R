#' Simulate chronological frequency distribution.
#' 
#' Simulates chronological distributions from a table of entities with defined date ranges, based on assumption of uniform probability
#'       between limits.
#' @param data Data table (or object that can be coerced to one) with, minimally, two numeric columns called Start and End.
#' @param weight Numeric vector: the weight to be applied to each row in `data`, or a constant weight to be applied to all.
#'      Defaults to 1.
#' @param context.fields Character vector specifying the column(s) in data which define the minimal stratigraphic entities to analyse. 
#'      Defaults to "SITE_C".
#' @param UoA Unit of Analysis: character vector of names of additional columns by which to group data when aggregating weights,
#'      on top of those specified in context.field. For example, should different taxa be lumped together when analysing bone remains 
#'      from a table of contexts? Defaults to NULL.
#' @param start.date Numeric: the start of time period to be considered. Defaults to 0.
#' @param end.date Numeric: the end of time period to be considered. Defaults to 2000.
#' @param bin.width Numeric: the resolution of the analysis, in units of time. Defaults to 100.
#' @param reps Integer: the number of times the simulation will be run. Defaults to 100.
#' @param summ Logical: should a summary results table be appended to the output? Defaults to TRUE. No real reason to change, except
#'      when called as part of a more complex function. Nb. the quantiles to use when summarising can be specified using a numeric
#'      vector called `quant.list`, which will be passed straight to the summarising function (`sim.summ`; see below) and defaults to
#'      c(0.025,0.25,0.5,0.75,0.975).
#' @param RoC Rate of Change. Logical: should rates of change between adjacent bins be calculated alongside the raw counts?
#' @return If summ=FALSE, a long-format data table with four or five named columns: 'rep.no', integer specifying simulation run; 'bin',
#'      character specifying chronological bin in terms of date range; 'bin.no' integer specifying number of bin, counting from earliest; 
#'      'count', numeric giving the number of entities (or total weight) assigned to the given bin in the given simulation run; 'RoC' 
#'      numeric giving the rate of change in count between this bin and the next (only present if RoC=TRUE).
#'      If summ=TRUE, a list with two named elements: "full" is as above; "summary" is a second long format data table with four named
#'      columns: 'bin', as above; 'V1', the relevant value for the given bin at a given quantile; 'quantile', 
#'      the quantile at which V1 is calculated; 'id', character specifying which column from "full" V1 is based on - either
#'      "count" or "RoC".
#' @export
#' @examples
#' date.ranges <- data.table(unit=c(1, 2, 3, 4), Start=c(450, 450, 600, 1000), End=c(700, 800, 650, 1200), frags=c(3, 6, 25, 1))
#' x <- date.simulate(date.ranges, weight=date.ranges$frags, context.fields="unit", start.date=500, end.date=1500, bin.width=50, reps=200)

date.simulate <- function(data, weight=1, context.fields=c("SITE_C"), UoA=NULL, start.date=0, end.date=2000, bin.width=100,
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
    for(i in 1:(length(breaks) - 1)) {
        labels[i] <- paste(breaks[i], breaks[i + 1], sep="-") #sets bin labels
    }
    params <<- paste("_", start.date, "-", end.date, "_by_", bin.width, "_x", reps, sep="") #saves char value with key parameters 
    
    #Perform simulation
    rep.no <- rep(1:reps, each = nrow(data))
    data <- cbind(rep.no, data) #recycles input data 'reps' times to provide frame for simulation 
    data[, sim := {x <- runif(nrow(data)); (x * (data[, End] - data[, Start])) + data[, Start]}] #simulates a date for each row
    data[, bin := cut(sim, breaks, labels=labels)] #finds the relevant bin for each simulated date
    data <- data[is.na(bin)==FALSE, j=list(count=sum(as.numeric(weight))), by=list(rep.no, bin)] #sums weights by bin and rep number
    
    #Prepare results table
    frame <- data.table(rep.no=rep(1:reps, each=length(labels)), bin.no=rep(1:length(labels), reps), bin=rep(labels, reps))
    results <- merge(frame, data, by=c("rep.no", "bin"), all=TRUE)
    results[is.na(results)] <- 0
    
    #Calculate rates of change, if necessary
    if(RoC==TRUE) {
        results[, RoC := c(diff(count), NA)]
        results[bin==labels[length(labels)], RoC := NA]
    }
    
    #Create summary dataset if required, and return results
    if(summ==TRUE) {
        summary <- sim.summ(results)
        results <- list(results, summary)
    }
    results
}