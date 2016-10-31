#' Calculates Rate of Change in simulated variables.
#'
#' Takes the standard "full" output from one of the archSeries simulation functions and adds columns giving rates of change between
#'      bins, for all value columns or a specified subset. Called within simulation functions by setting their 'RoC' arguments to
#'      TRUE, but can also be called directly after the fact.
#' @param results Data table resembling the first element ("full") of the output from one of the archSeries simulation functions.
#' @param field.list Character vector of columns for which RoC should be calculated. Defaults to NULL, in which case all
#'      suitable columns are used.
#' @param type Character: "a" for absolute change between bins; "r" for rate as a proportion of current value. Defaults to "a".
#' @return The input data table with the addition of one new RoC column for each original value column.
#' @export
#' @examples
#' date.ranges <- data.table(ID=c(1, 2, 3, 4), Start=c(450, 450, 600, 900), End=c(700, 800, 650, 1200))
#' results <- date.simulate(date.ranges, start.date=400, end.date=1200, bin.width=200)
#' x <- RoC.fun(results[[1]])

RoC.fun <- function(results, field.list=NULL, type="a") {
    bin.no <- NULL

    # Find column names, if not supplied
    if(is.null(field.list)) {field.list <- colnames(results)[!colnames(results) %in% c("rep.no", "bin", "bin.no", "catch", "effort")]}

    # Set up names for RoC columns
    new.fields <- paste("RoC", field.list, sep=".")

    # Set what to divide results by (either 1, or the column in question)
    abs.denom <- c(1, 1)
    if(type=="r") {denom <- field.list} else {denom <- rep("abs.denom", length(field.list))}

    # Calulate rates of change for each column in turn
    for(i in 1:length(field.list)) {
        results[, assign("a", new.fields[i]) := c(diff(get(field.list[i])) / get(denom[i])[1:(length(get(denom[i])) - 1)], 0)]
        results[bin.no==max(results$bin.no), assign("a", new.fields[i]) := 0] #Set values for final bin to 0
        results[(get(new.fields[i]))==Inf, assign("a", new.fields[i]) := max(results[is.finite(get(new.fields[i])), get(new.fields[i])])]
        results[(get(new.fields[i]))==-Inf, assign("a", new.fields[i]) := min(results[is.finite(get(new.fields[i])), get(new.fields[i])])]
        results[is.na(get(new.fields[i])), assign("a", new.fields[i]) := 0]
    }
    results
}
