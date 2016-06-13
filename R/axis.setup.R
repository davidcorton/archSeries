#' Set up axes for archSeries plotting functions.
#' 
#' A utility function designed to be used within the various plot functions in this package, but which can also be used alone to set up 
#'      axes based on simulation data prior to plotting anything.
#' @param results A list resembling the output from date.simulate or a related function, or a data table resembling one of the components
#'      thereof.
#' @param field.list A character vector of names of variables that are intended to be plotted. Depending on 'results' and 'type' this may 
#'      mean column names in a full simulation output or values of 'id' in a summary table. Defaults to NULL, in which case all suitable 
#'      variables in 'results' are used. Nb. variables entitled "catch" or "effort" will be ignored, due to their roles in output from 
#'      cpue.
#' @param lab.sp Integer: intervals at which to place bin labels. Defaults to 1, i.e. labelling every bin.
#' @param ylab Character: label for the y-axis. Defaults to "Estimated frequency density" (which is the only reason for making it a 
#'      formal argument here rather than just allowing it to be passed to plot via ...).
#' @param ylim Numeric: an easy way to override the built-in scaling in plot - if a vector of length 1 is passed it will 
#'      be converted into c(0, ylim) to be passed to the ylim argument in plot. Alternatively a vector of length 2 will be passed 
#'      straight to plot as is. Defaults to NULL, in which case the built-in scaling in plot takes over.
#' @param type Integer: the type of simulation output to be plotted. 1 = full results (for lines.chron or box.chron), 2 = summary results 
#'      (for poly.chron).
#' @param axis.lab Logical. Should labels be plotted for the x axis?
#' @param ... other graphical arguments to be passed to plot.
#' @return None.
#' @export
#' @examples
#' date.ranges <- data.table(unit=c(1, 2, 3, 4), Start=c(450, 450, 600, 1000), End=c(700, 800, 650, 1200), frags=c(3, 6, 25, 1))
#' x <- date.simulate(date.ranges, weight=date.ranges$frags, context.fields="unit", bin.width=50, reps=200, summ=FALSE)
#' axis.setup(x, lab.sp=2, type=1)

axis.setup <- function(results, field.list=NULL, lab.sp=1, ylab="Estimated frequency density", ylim=NULL,
                       type=1, axis.lab=TRUE,...) {
    if(class(results)[1]=="list") {results <- results[[type]]}
    minmaxer <- numeric(0)
    if(type==1) {
        if(is.null(field.list)) {field.list <- colnames(results)[!colnames(results) %in% c("bin", "bin.no", "rep.no", "catch", "effort")]}
        for(i in 1:length(field.list)) {minmaxer <- c(minmaxer, results[, get(field.list[i])])}
        minmaxer <- cbind(minmaxer, unique(results$bin.no))
    } else {
        if(is.null(field.list)) {field.list <- unique(results$id)}
        minmaxer <- as.matrix(cbind(results[id %in% field.list, V1], 1:length(unique(results$bin))))
    }
    if(!is.null(ylim)) {
        if(length(ylim)==1) {ylim <- c(0, ylim)}
    }
    plot(minmaxer[, 2], minmaxer[, 1], xlab="", xaxt="n", ylab=ylab, type="n", ylim=ylim, ...)
    names <- unique(results$bin)
    ticks <- seq(1, length(names), by=lab.sp)
    if(axis.lab==TRUE) {labels <- names[ticks]} else {labels <- FALSE}
    axis(1, at=ticks, labels=labels, las=2, ...)
}