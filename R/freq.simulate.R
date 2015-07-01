#' Simulate chronological frequency distribution plus dummy set.
#' 
#' Simulates chronological distributions from a table of entities with defined date ranges, based on assumption of uniform probability
#'       between limits, then simulates a dummy set of the same size drawing from a specified distribution (i.e. a convenient wrapper
#'       for date.simulate and dummy.simulate).
#' @param data Data table (or object that can be coerced to one) with, minimally, two numeric columns called Start and End.
#' @param probs Numeric vector defining a null model from which to sample the dummy set. Will be recycled up to nrow(data), so passing a 
#'      single value results in a uniform null model. Defaults to 1.
#' @param weight Numeric vector: the weight to be applied to each row in `data` (and to its counterpart in the dummy set), or a constant 
#'      weight to be applied to all. Defaults to 1.
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
#' @return A list with two named elements:
#'      "full" is a long-format data table with at least five named columns: 'rep.no', integer specifying simulation run; 'bin',
#'      character specifying chronological bin in terms of date range; 'bin.no' integer specifying number of bin, counting from earliest; 
#'      'count', numeric giving the number of entities (or total weight) assigned to the given bin in the given simulation run; 'dummy', 
#'      giving the number of entities (or total weight) assigned to a bin in the dummy version of a given simulation run. If RoC=TRUE there 
#'      will be two more columns: 'RoC.count' and 'RoC.dummy' give the rate of change between this bin and the next for 'count' and 'dummy' 
#'      respectively.
#'      "summary" is a second long format data table with four named columns: 'bin', as above; 'V1', the relevant value for the given bin 
#'      at a given quantile; 'quantile', the quantile at which V1 is calculated; 'id', character specifying which column from "full" V1 is 
#'      based on: "count", "dummy", "RoC.count", or "RoC.dummy".
#' @export
#' @examples
#' date.ranges <- data.table(unit=c(1, 2, 3, 4), Start=c(450, 450, 600, 1000), End=c(700, 800, 650, 1200), frag.count=c(3, 6, 25, 1))
#' x <- freq.simulate(date.ranges, weight=data$frag.count, context.fields="unit", start.date=500, end.date=1500, bin.width=50, reps=200)

freq.simulate <- function(data, probs=1, weight=1, context.fields=c("SITE_C"), UoA=NULL, quant.list=c(0.025,0.25,0.5,0.75,0.975),
                          start.date=0, end.date=2000, bin.width=100, reps=100, RoC=FALSE, ...) {
    #Load required packages
    require(data.table)
    
    #Tidy up input data; apply filters
    data <- data.table(cbind(data, weight)) #appends weights to list of date ranges, recycling if necessary (e.g. for uniform weight) 
    data <- data[End >= start.date & Start <= end.date]  #drops records outside the date range FROM BOTH SIMULATION SETS
    
    #Aggregate data
    agg.list <- c(context.fields, "Start", "End", UoA)
    data <- data[, j=list(weight=sum(as.numeric(weight))), by=agg.list]
    if(length(weight)==1) {data[, weight := 1]}
    
    #Reset bin.width based on probs, if necessary
    if(is.vector(probs)==TRUE & length(probs)>1) {bin.width <- (end.date - start.date) / length(probs)}  #if probs supplied, use to set bin.widths
    if(sum(class(probs)=="data.frame")==1) {bin.width <- (end.date - start.date) / nrow(probs)}  #likewise if supplied as data.frame/data.table
    
    #Simulate from real data, then generate dummy set. Merge the two together.
    real <- date.simulate(data=data, weight=data$weight, context.fields=context.fields, UoA=UoA, bin.width=bin.width, start.date=start.date, end.date=end.date, reps=reps, summ=FALSE)
    dummy <- dummy.simulate(weight=data$weight, probs=probs, breaks=breaks, start.date=start.date, end.date=end.date, bin.width=bin.width, reps=reps, summ=FALSE)    
    results <- merge(real, dummy, by=c("rep.no", "bin.no", "bin"), all=TRUE)
    
    #Calculate rate of change variables (could be done within core functions, but faster to loop through together here)
    if(RoC==TRUE) {
        results[, RoC.count := c(diff(count), NA)]
        results[, RoC.dummy := c(diff(dummy), NA)]
        results[bin==unique(bin)[length(unique(bin))], RoC.count := NA]
        results[bin==unique(bin)[length(unique(bin))], RoC.dummy := NA]
    }
    
    #Create summary dataset and return results
    summary <- sim.summ(results)
    list(results, summary)
}