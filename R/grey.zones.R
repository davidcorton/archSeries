#' Plot boxes to indicate areas of low sample size.
#' 
#' A utility function sometimes called by main archSeries plotting functions when dealing with output from cpue. Can in theory also be used 
#'      alone after a call to axis.setup.
#' @param boxes A list resembling the third element of the output from a cpue call ("small.n").
#' @param colours Character vector: colours to be used for boxes at each level of small sample. No default.
#' @param opacity Numeric vector: opacity of boxes at each level of small sample. No default
#' @param ylim Numeric: an easy way to override the built-in scaling in plot - if a vector of length 1 is passed it will 
#'      be converted into c(0, ylim) to be passed to the ylim argument in plot. Alternatively a vector of length 2 will be passed 
#'      straight to plot as is. Defaults to NULL, in which case the built-in scaling in plot takes over.
#' @return None.
#' @export
#' @examples
#' dates <- data.table(unit=c(1, 2, 3, 4), Start=c(450, 450, 600, 1000), End=c(700, 800, 650, 1200), frags=c(3, 6, 2, 1), vol=c(40, 40, 40, 40))
#' x <- cpue(dates, dates, dates$frags, dates$vol, context.fields=NULL, small.n=1, reps=1000)
#' axis.setup(x)
#' grey.zones(x$small.n, "grey30", 200, ylim=130)


grey.zones <- function(boxes, colours, opacity, ylim=0) {
    if(is.null(ylim)) {ylim <- 0}
    id <- 1:length(boxes)
    cols <- data.frame(id, colours=colours, opacity=opacity)
    a <- col2rgb(cols[,2])
    b <- character()
    for(i in 1:nrow(cols)) {
        b[i] <- rgb(a[1,i], a[2,i], a[3,i], cols[i,3], maxColorValue=255) #Builds opacity into each colour
    }
    for(i in 1:length(boxes)) {
        for(j in 1:length(boxes[[i]])) {
            x <- boxes[[i]][[j]]
            x$weak.y[2:3] <- pmax(x$weak.y[2:3], ylim)
            with(x, polygon(weak.x, weak.y, col=b[i], border=NA))
        }
    }
}
