#' Calculate aoristic sum.
#'
#' Calculates aoristic sum from a table of entities with defined date ranges, based on assumption of uniform probability between limits.
#' @param data Data table with (minimally) two numeric columns called Start and End.
#' @param weight Numeric vector. The weight to be applied to each row in `data`, or a constant weight to be applied to all.
#'      Defaults to 1.
#' @param start.date Numeric. Start of time period to be considered. Defaults to lowest value in data$Start.
#' @param end.date Numeric. End of time period to be considered. Defaults to highest value in data$End.
#' @param bin.width Numeric. The resolution of the analysis, in units of time. Defaults to 100.
#' @return data table with two named columns: 'bin', a character vector specifying the date range represented by each chronological bin;
#'      'aorist', a numeric vector giving the total probability mass assigned to each bin.
#' @export
#' @examples
#' date.ranges <- data.table(Start=c(450, 450, 600, 1000, 1200), End=c(700, 800, 650, 1200, 1550), frag.count=c(3, 6, 25, 1, 8))
#' x <- aorist(date.ranges, weight=date.ranges$frag.count, 500, 1500, bin.width=50)

aorist <- function(data, weight=1, start.date=NULL, end.date=NULL, bin.width=100) {
    #Load required package
    require(data.table)

    #Tidies up input data
    data <- data.table(cbind(data, weight)) #appends weights to list of date ranges, recycling if necessary (e.g. for uniform weight)

    #Read start and end dates from input data if not specified
    if(is.null(start.date)) {
        start.date <- min(data$Start)
    }
    if(us.null(end.date)) {
        end.date <- max(data$End)
    }
    data <- data[End >= start.date & Start <= end.date] #excludes ranges that fall entirely outside the study period

    #Set up columns for duration and for weight per year
    data[, duration := End - Start]
    data[, weight.per.unit := weight / duration]

    #Set up breaks and labels
    breaks <<- seq(start.date, end.date, bin.width) #creates and saves vector of breaks
    labels <- numeric(0)
    for(i in 1:(length(breaks) - 1)) {
        labels[i] <- paste(breaks[i], breaks[i + 1], sep="-") #sets bin labels
    }

    #Set frame for results
    aorist <- data.table(bin = labels, bin.no = 1:length(labels), aorist = 0)

    #Cycle through bins, assigning probability mass to cases where appropriate

    for(i in 1:length(labels)) {
        bin.1 <- breaks[i] #Find start date of bin
        bin.2 <- breaks[i + 1] #Find end date of bin
        data[, assign("a", labels[i]) := 0]
        data[Start >= bin.1 & Start < bin.2, assign("a", labels[i]) := (bin.2 - Start) * weight.per.unit]
        data[End > bin.1 & End <= bin.2, assign("a", labels[i]) := (End - bin.1) * weight.per.unit]
        data[Start < bin.1 & End > bin.2, assign("a", labels[i]) := bin.width * weight.per.unit]
        data[Start >= bin.1 & End <= bin.2, assign("a", labels[i]) := as.double(weight)]
        aorist$aorist[i] <- sum(data[, get(labels[i])], na.rm=TRUE)
    }
    aorist
}
