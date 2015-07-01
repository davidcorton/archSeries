#' Convert simulated mortality data to survivorship format.
#' 
#' Takes simulated frequencies for a series of ordinal classes and calculates survivorship across that series. Intended for mortality data. 
#' Survivorship currently defined as proportion surviving to *end* of given class; in future plan to add an argument allowing this to be set 
#' to *start* of class, as is standard practice in human demographics.
#' 
#' @param mortality Output from a date.simulate, dummy.simulate, or freq.simulate call: either a data table with columns 'bin', 'bin.no' 
#'      and 'rep.no' (plus at least one column of frequencies) or a list whose first item is such a data table.
#' @param field.list Character vector of columns in 'mortality' which contain the mortality data to convert. Defaults to NULL, in which 
#'      case all frequency columns are used.
#' @param quant.list Numeric vector of quantiles to be calculated in a summary table. Defaults to c(0.025,0.25,0.5,0.75,0.975).
#' @return A list with two named elements:
#'      "full" is a long-format data table with at least four named columns: 'rep.no', integer specifying simulation run; 'bin',
#'      character specifying chronological bin in terms of date range; 'bin.no' integer specifying number of bin, counting from earliest; 
#'      then for each input frequency column a column called 'survive.[input column name], giving survivorship (out of 1) to the end of the
#'      given bin in the given simulation run.
#'      "summary" is a second long format data table with four named columns: 'bin', as above; 'V1', the relevant value for the given bin 
#'      at a given quantile; 'quantile', the quantile at which V1 is calculated; 'id', character specifying which column from "full" V1 is 
#'      based upon.
#' @export
#' @examples
#' # Simulating a sample of 50 mandibles from an ideal dairy herd model, then calculating survivorship
#' dairy.model <- c(0.53, 0.05, 0.03, 0.04, 0.07, 0.05, 0.04, 0.09, 0.10)
#' payne.breaks <- c(0, 2, 6, 12, 24, 36, 48, 2, 96, 20)
#' sim.ages <- dummy.simulate(50, probs=dairy.model, breaks=payne.breaks, reps=1000)
#' sim.survive <- surv.convert(sim.ages)

surv.convert <- function(mortality, field.list=NULL, quant.list=c(0.025, 0.25, 0.5, 0.75, 0.975)) {
    #Load required package
    require(data.table)
    
    #Select only full simulation output (if applicable); establish field names to use
    if(class(mortality)[1]=="list") {mortality <- mortality[[1]]}
    if(is.null(field.list)==TRUE) {field.list <- colnames(mortality)[!colnames(mortality) %in% c("bin", "bin.no", "rep.no")]}
    
    #Build new data table for results
    frame <- data.table(rep.no=rep(1:max(mortality$rep.no), each=(max(mortality$bin.no) + 1)), bin.no=rep(1:(max(mortality$bin.no) + 1)), 
                        bin=rep(c("Start", unique(as.character(mortality$bin))), max(mortality$rep.no)))
    column.names <- paste("survive.", field.list, sep="")
    
    #Calculate survivorship for each column, run, and bin (nested in that order)
    for(k in 1:length(field.list)) {
        n <- mortality[rep.no==1 & !is.na(get(field.list[k])), sum(get(field.list[k]))] #Find total sample size for current column
        frame[bin.no==1, assign("a", column.names[k]) := n] #Create column and fill in start value for each run
        frame[!bin.no==1, assign("a", column.names[k]) := mortality[, get(field.list[k])]]    
        for(i in 1:max(frame$rep.no)) {
            frame[rep.no==i, assign("a", column.names[k]) := mortality[rep.no==i, n - diffinv(get(field.list[k]))]]
        }
        frame[, assign("a", column.names[k]) := get(column.names[k]) / n]
    }
    
    #Create summary dataset and return results
    summary <- sim.summ(frame, quant.list=quant.list)
    list(frame, summary)
}


