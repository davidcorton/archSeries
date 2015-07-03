#' Plot medians and confidence polygons for simulation results.
#' 
#' Plots defined confidence intervals (as polygons) and medians (as lines) for output from date.simulate or an associated function.
#' @param results A list resembling the output from date.simulate or a related function, or a data table resembling the second component 
#'      thereof - i.e. the summary simulation results.
#' @param field.list A character vector of values in results$id to be plotted. Defaults to NULL, in which case all suitable 
#'      values are used.
#' @param quant Numeric vector of length 2: lower and upper quantiles to define confidence zone. Defaults to c(0.025, 0.975), i.e. 95%.
#' @param col.list Character vector: colours to be used for each column plotted. Defaults to c("darkred", "darkgreen", "blue", "grey", 
#'      "goldenrod"). Nb. if more than five groups are plotted then the default colours will start to recycle.
#' @param opacity Numeric: opacity of each line. Defaults to 80.
#' @param ylim Numeric: an easy way to override the built-in scaling in plot - if a vector of length 1 is passed it will 
#'      be converted into c(0, ylim) to be passed to the ylim argument in plot. Alternatively a vector of length 2 will be passed 
#'      straight to plot as is. Defaults to NULL, in which case the built-in scaling in plot takes over.
#' @param med.line Logical: should a median line be plotted on top of eachconfidence polygon? Defaults to TRUE.
#' @param border Character: what type of border should the polygons have? Defaults to NA, i.e. no border.
#' @param small.n Character vector of colours to be used for boxes marking periods of low "effort", when plotting cpue results. Has no 
#'      effect if 'results' doesn't have a third element called "small.n". Defaults to NULL, in which case no boxes are plotted.
#' @param small.n.op Numeric vector of length equal to 'small.n' (or otherwise recycled to that length) specifying opacity for small.n 
#'      boxes. Defaults to 126.
#' @param add Logical: should data be added to current plot, or should axis.setup be called to start a new plot? Defaults to FALSE.
#' @param legend Logical: should an automatic legend of column names and corresponding colours be plotted? Defaults to TRUE.
#' @param ... Other graphical arguments to be passed to plot. Nb. (a) includes special arguments for axis.setup (currently just 'lab.sp'), 
#'      (b) 'ylab' will default to "Estimated frequency density", as per axis.setup, unless specified here.
#' @return None.
#' @export
#' @examples
#' date.ranges <- data.table(unit=c(1, 2, 3, 4), Start=c(450, 450, 600, 1000), End=c(700, 800, 650, 1200), frags=c(3, 6, 25, 1))
#' x <- freq.simulate(date.ranges, weight=date.ranges$frags, context.fields="unit", bin.width=50, reps=200)
#' poly.chron(x, field.list=c("dummy", "count"))
#' 
poly.chron <- function(results, field.list=NULL, quant=c(0.025, 0.975), col.list=c("darkred", "darkgreen", "blue", "grey", "goldenrod"),
                       opacity=80, ylim=NULL, med.line=TRUE, border=NA, small.n=NULL, small.n.op=126, add=FALSE, legend=TRUE, ...) {
    boxes <- NULL
    if(class(results)[1]=="list") {
        if(length(results)==3) {boxes <- results$small.n}
        results <- results[[2]]
    }
    if(is.null(field.list)==TRUE) {field.list <- unique(results$id)}
    if(add==FALSE) {axis.setup(results, field.list=field.list, ylim=ylim, type=2, ...)}
    if(!is.null(small.n) & !is.null(boxes)) {grey.zones(boxes, small.n, small.n.op, ylim[length(ylim)])} #Sets up boxes to highlight small n
    plist <- data.frame(field.list, col.list[1:length(field.list)], opacity)
    a <- col2rgb(plist[,2])
    b <- character()
    for(i in 1:nrow(plist)) {
        b[i] <- rgb(a[1, i], a[2, i], a[3, i], plist[i, 3], maxColorValue=255)
    }
    x <- c(1:length(unique(results$bin)), length(unique(results$bin)):1)
    for(i in 1:length(field.list)) {
        y <- c(results[id==field.list[i] & quantile==quant[1], V1], rev(results[id==field.list[i] & quantile==quant[2], V1]))
        skip <- length(x) + 1       
        if(substr(field.list[i], 1, 3)=="RoC") {skip <- length(unique(results$bin))}
        polygon(x[!x==skip], y[!x==skip], col=b[i], border=border)
        if(med.line==TRUE) {
            with(results[id==field.list[i] & quantile==0.500], lines(1:length(unique(bin)), V1, col=col.list[i], lwd=3))
            with(results[id==field.list[i] & quantile==0.500], points(1:length(unique(bin)), V1, col=col.list[i], pch=19))
        }
    }
    if(legend==TRUE) {with(results, legend("topright", legend=field.list, fill=b, bty="n"))}
}
