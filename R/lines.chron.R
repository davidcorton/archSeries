#' Plot all simulation runs as lines.
#'
#' Plots every single simulation run for each specified variable as a separate semi-transparent line.
#' @param results A list resembling the output from date.simulate or a related function, or a data table resembling the first component
#'      thereof - i.e. the full simulation results.
#' @param field.list A character vector of columns in 'results' to be plotted. Defaults to NULL, in which case all suitable
#'      columns in 'results' are used. Nb. variables entitled "catch" or "effort" will be ignored, due to their roles in output from
#'      cpue.
#' @param col.list Character vector: colours to be used for each column plotted. Defaults to c("darkred", "darkgreen", "blue", "grey",
#'      "goldenrod"). Nb. if more than five columns are plotted then the default colours will start to recycle.
#' @param opacity Numeric: opacity of each line. Defaults to 20.
#' @param ylim Numeric: an easy way to override the built-in scaling in plot - if a vector of length 1 is passed it will
#'      be converted into c(0, ylim) to be passed to the ylim argument in plot. Alternatively a vector of length 2 will be passed
#'      straight to plot as is. Defaults to NULL, in which case the built-in scaling in plot takes over.
#' @param small.n Character vector of colours to be used for boxes marking periods of low "effort", when plotting cpue results. Has no
#'      effect if 'results' doesn't have a third element called "small.n". Defaults to NULL, in which case no boxes are plotted.
#' @param small.n.op Numeric vector of length equal to 'small.n' (or otherwise recycled to that length) specifying opacity for small.n
#'      boxes. Defaults to 126 (i.e. ~50pc).
#' @param add Logical: should data be added to current plot, or should axis.setup be called to start a new plot? Defaults to FALSE.
#' @param legend Logical: should an automatic legend of column names and corresponding colours be plotted? Defaults to TRUE.
#' @param ... Other graphical arguments to be passed to plot. Nb. (a) includes special arguments for axis.setup (currently just 'lab.sp'),
#'      (b) 'ylab' will default to "Estimated frequency density", as per axis.setup, unless specified here.
#' @return None.
#' @export lines.chron
#' @examples
#' date.ranges <- data.table(Start=c(450, 450, 600), End=c(700, 800, 650), frags=c(3, 6, 25))
#' x <- freq.simulate(date.ranges, weight=date.ranges$frags, bin.width=50, reps=200, summ=FALSE)
#' lines.chron(x)

lines.chron <- function(results, field.list=NULL, col.list=c("darkred", "darkgreen", "blue", "grey", "goldenrod"), opacity=20,
                         ylim=NULL, small.n=NULL, small.n.op=126, add=FALSE, legend=TRUE, ...) {
    rep.no <- boxes <- NULL
    if(class(results)[1]=="list") {
        if(length(results)==3) {boxes <- results$small.n}
        results <- results[[1]]
    }
    if(is.null(field.list)==TRUE) {field.list <- colnames(results)[!colnames(results) %in% c("bin", "bin.no", "rep.no", "catch", "effort", "n")]}
    if(add==FALSE) {axis.setup(results, field.list=field.list, ylim=ylim, ...)}
    if(!is.null(small.n) & !is.null(boxes)) {grey.zones(boxes, small.n, small.n.op, ylim[length(ylim)])} #Sets up boxes to highlight small n
    plist <- data.frame(field.list, col.list[1:length(field.list)], opacity)
    a <- col2rgb(plist[, 2])
    b <- character()
    for(i in 1:nrow(plist)) {
        b[i] <- rgb(a[1, i], a[2, i], a[3, i], plist[i, 3], maxColorValue=255)
    }
    for(i in 1:max(results$rep.no)) {
        for(j in 1:length(field.list)) {
            lines(results[rep.no==i, get(field.list[j])], col=b[j])
        }
    }
    if(legend==TRUE) {with(results, legend("topright", legend=field.list, fill=col.list[1:length(b)], bty="n"))}
}
