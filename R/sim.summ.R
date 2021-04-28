#' Summarise basic output from simulation functions.
#'
#' A utility function designed to be used within the simulation functions in this package, but which can also be used alone if useful.
#' @param results A data table resembling the full-form output from date.simulate or a related function.
#' @param summ.col Character: which columns(s) in 'results' should be summarised? Defaults to NULL, in which case all columns not called
#'      "bin", "bin.no" or "rep.no" are summarised.
#' @param quant.list Numeric vector of quantiles to be calculated in a summary table. Defaults to c(0.025,0.25,0.5,0.75,0.975).
#' @return A long-format data table with four named columns: 'bin', character specifying chronological bin in terms of date range; 'V1',
#'      the relevant value for the given bin at a given quantile; 'quantile', the quantile at which V1 is calculated; 'id', character
#'      specifying which column from 'results' V1 is based on: potentially "count", "dummy", "RoC.count", or "RoC.dummy" depending on the
#'      function call which created 'results'.
#' @export
#' @examples
#' date.ranges <- data.table(Start=c(450, 450, 600), End=c(700, 800, 650), frag.count=c(3, 6, 25))
#' x <- date.simulate(date.ranges, weight=date.ranges$frag.count, bin.width=50, reps=200, summ=FALSE)
#' x.summary <- sim.summ(x)

sim.summ <- function(results, summ.col=NULL, quant.list=c(0.025, 0.25, 0.5, 0.75, 0.975)) {
    bin <- id <- NULL

    if(is.null(summ.col)==TRUE) {summ.col <- colnames(results)[!colnames(results) %in% c("rep.no", "bin", "bin.no")]}

    #Create summary tables
    for(i in 1:length(summ.col)) {
        x <- results[, quantile(get(summ.col[i]), probs=quant.list, na.rm=TRUE), by=bin] #calculate quantiles
        x[, quantile := rep(quant.list, length.out = nrow(x))] #create column to specify quantiles
        x[, id := summ.col[i]] #create column to specify variable
        if(i==1) {summary <- data.table()}
        summary <- rbind(summary, x)
    }

    #Return summary table
    summary
}
