#' Plots aorist output as a barplot.
#'
#' Just a wrapper for barplot with some tweaks added, e.g. to make bars line up with data plotted by other archSeries functions.
#' @param aorist The output from an aorist call, or any data table with a character or factor column called 'bin' and a numeric column
#'      called aorist.
#' @param col Character: colour of bars. Defined formally here rather than passed via ... simply so that it can be combined with 'opacity'.
#' @param opacity Numeric: opacity of each line. Defaults to 80.
#' @param lab.sp Integer: intervals at which to place bin labels. Defaults to 1, i.e. labelling every bin.
#' @param add Logical: should data be added to current plot, or should axis.setup be called to start a new plot? Defaults to FALSE.
#' @param ... Other graphical arguments to be passed to barplot.
#' @return None.
#' @export
#' @examples
#' date.ranges <- data.table(Start=c(450, 450, 600), End=c(700, 800, 650), frag.count=c(3, 6, 25))
#' x <- aorist(date.ranges, weight=date.ranges$frag.count, 500, 1500, bin.width=50)
#' aorist.plot(x, col="grey60", ylab="Total probability density")

#Function for plotting aorist output (this is really just barplot with some tweaks
#added, e.g. to make bars line up with points plotted by other functions)

aorist.plot <- function(aorist, col="grey", opacity=255, lab.sp=1, add=FALSE, ...) {
    col <- col2rgb(col)
    col <- rgb(col[1, ], col[2, ], col[3, ], opacity, maxColorValue=255)
    with(aorist, barplot(aorist, space=c(0.5, rep(0, length(bin) - 1)), add=add, ...))
    if(add==FALSE) {
        ticks <- seq(1, nrow(aorist), by=lab.sp)
        axis(1, at=ticks, labels=aorist$bin[ticks], las=2)
    }
}
