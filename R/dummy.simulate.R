#' Simulate frequency distribution based on null model.
#' 
#' Simulates a 'dummy' chronological distribution within specified date limits by sampling from within a distribution defined by an input
#' vector (typically a null model). Designed for use with date.simulate, particularly within wrapper functions like freq.simulate. The
#' idea here is to simulate chronological distributions based on the same number of entities (with the same weights) as a date.simulate
#' call, but unconstrained by the known date ranges.
#' @param weight Numeric vector whose length indicates the number of entities to be drawn from the null model, and whose values indicate
#'      their relative weights. Typically matches the argument of the same name in a date.simulate call (e.g. within a freq.simulate call). 
#'      Alternatively, if a single integer is passed this is taken as the number entities to draw, with equal weights.
#' @param probs Numeric vector defining a null model from which to sample. Will be recycled up to the length of 'weight' (or the value
#'      of weight if the latter has length 1) so passing a single value results in a uniform null model. Defaults to 1.
#' @param breaks Numeric vector of breaks between bins. Designed to allow more complex functions to pass breaks between date.simulate and
#'      dummy.simulate. Where provided, will override 'bin.width'.
#' @param start.date Numeric: the start of time period to be considered. Defaults to 0.
#' @param end.date Numeric: the end of time period to be considered. Defaults to 2000.
#' @param bin.width Numeric: the resolution of the analysis, in units of time. Ignored if 'breaks' is passed. Defaults to 100.
#' @param reps Integer: the number of times the simulation will be run. Defaults to 100.
#' @param summ Logical: should a summary results table be appended to the output? Defaults to TRUE. No real reason to change, except
#'      when called as part of a more complex function. Nb. the quantiles to use when summarising can be specified using a numeric
#'      vector called `quant.list`, which will be passed straight to the summarising function (`sim.summ`; see below) and defaults to
#'      c(0.025,0.25,0.5,0.75,0.975).
#' @param RoC Rate of Change. Logical: should rates of change between adjacent bins be calculated alongside the raw counts?
#' @return If summ=FALSE, a long-format data table with four or five named columns: 'rep.no', integer specifying simulation run; 'bin',
#'      character specifying chronological bin in terms of date range; 'bin.no' integer specifying number of bin, counting from earliest; 
#'      'dummy', numeric giving the number of entities (or total weight) assigned to the given bin in the given simulation run; 'RoC' 
#'      numeric giving the rate of change in count between this bin and the next (only present if RoC=TRUE).
#'      If summ=TRUE, a list with two named elements: "full" is as above; "summary" is a second long format data table with four named
#'      columns: 'bin', as above; 'V1', the relevant value for the given bin at a given quantile; 'quantile', 
#'      the quantile at which V1 is calculated; 'id', character specifying which column from "full" V1 is based on - either
#'      "dummy" or "RoC".
#' @export
#' @examples
#' date.ranges <- data.table(unit=c(1, 2, 3, 4), Start=c(450, 450, 600, 1000), End=c(700, 800, 650, 1200), frag.count=c(3, 6, 25, 1))
#' x <- dummy.simulate(date.ranges$frag.count, start.date=500, end.date=1500, bin.width=50, reps=200)

dummy.simulate <- function(weight, probs=1, breaks=NULL, start.date=0, end.date=2000, bin.width=100, reps=100, RoC=FALSE, summ=TRUE, ...) {
    #Load required package
    require(data.table)
    
    #Tidy up input data and apply filters
    if(is.vector(weight)==1 & length(weight)==1) {weight <- rep(1, weight)} #if weight is a single value, use as number of entities
    dummy <- data.table(weight) #convert weights to data table format, if necessary.
    
    #Set up breaks and labels
    if(is.null(breaks)==TRUE) {breaks <- seq(start.date, end.date, bin.width)} #if breaks not specified, sets them based on other arguments
    labels <- numeric(0)
    for(i in 1:(length(breaks) - 1)) {
        labels[i] <- paste(breaks[i], breaks[i + 1], sep="-") #sets bin labels based on breaks
    }
    probs <- cbind(1:length(labels), probs)
    
    #Perform simulation
    rep.no <- rep(1:reps, each=nrow(dummy))
    dummy <- cbind(rep.no, dummy) #recycles input data 'reps' times to provide frame for simulation 
    dummy[, bin := sample(labels, size=nrow(dummy), replace=TRUE, prob=probs[, 2])]
    dummy <- dummy[, j=list(dummy=sum(as.numeric(weight))), by=list(rep.no, bin)] #sums weights by bin and rep number
    
    #Prepares and returns results table
    frame <- data.table("rep.no"=rep(1:reps, each=length(labels)), "bin.no"=rep(1:length(labels), reps), "bin"=rep(labels, reps))
    results <- merge(frame, dummy, by=c("rep.no", "bin"), all=TRUE)
    results[is.na(dummy), dummy := 0]
    results <- results[order(rep.no, bin.no)]
    
    #Calculate rates of change, if necessary
    if(RoC==TRUE) {
        results[, RoC := c(diff(dummy), NA)]
        results[bin==labels[length(labels)], RoC := NA]
    }
    
    #Create summary dataset if required, and return results
    if(summ==TRUE) {
        summary <- sim.summ(results)
        results <- list(results, summary)
    }
    results
}