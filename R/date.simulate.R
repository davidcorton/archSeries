#' Simulate date distributions.
#'
#' Simulates chronological distributions from a table of entities with defined date ranges, based on assumption of uniform probability
#'       between limits, then (optionally) simulates a dummy set of the same size drawing from a specified distribution.
#' @param data Data table (or object that can be coerced to one) with, minimally, two numeric columns called Start and End.
#' @param probs Numeric vector defining a null model from which to sample the dummy set. Will be recycled up to nrow(data), so passing a
#'      single value results in a uniform null model. If length > 1 then probs is used to set number of bins, overriding bin.width.
#'      Defaults to 1.
#' @param weight Numeric vector: the weight to be applied to each row in `data` (and to its counterpart in the dummy set), or a constant
#'      weight to be applied to all. Defaults to 1.
#' @param ds.fun Function: the summary function to be applied to the entities in each bin during each simulation run. Defaults to sum, i.e.
#'      calculates frequency distributions, but can be set to e.g. mean or median to deal with e.g. metrical data.
#' @param real Logical: should the date distribution of the empirical data be simulated? Defaults to TRUE.
#' @param dummy Logical: should a dummy set be simulated in addition to the empirical data? Defaults to FALSE.
#' @param comp.field Character: optional name of column to be used to subset 'data'. Defaults to NULL.
#' @param comp.values Optional vector specifying values of 'comp.field' to compare. Defaults to NULL, in which case all
#'      unique values (exlcuding blanks and NAs) are compared if comp.field is not NULL.
#' @param context.fields Character vector specifying the column(s) in data which define the minimal stratigraphic entities to analyse.
#'      Add more column names if you want to group the data by additional criteria prior to simulations - For example, should different
#'      taxa be treated separately rather than lumped together when analysing bone remains from a table of contexts? Defaults to "SITE_C".
#' @param quant.list Numeric vector of quantiles to be calculated in a summary table. Defaults to c(0.025,0.25,0.5,0.75,0.975).
#' @param start.date Numeric: the start of time period to be considered. Defaults to 0.
#' @param end.date Numeric: the end of time period to be considered. Defaults to 2000.
#' @param bin.width Numeric: the resolution of the analysis, in units of time. Defaults to 100.
#' @param reps Integer: the number of times the simulation will be run. Defaults to 100.
#' @param RoC Rate of Change. Character: how should rates of change between adjacent bins be calculated alongside the raw counts? In
#'      absolute terms ("a") or relative to the current bin ("r"). Defaults to NULL, in which case not calculated at all.
#' @param summ Logical: should a summary table be calculated (allowing plotting with poly.chron, for example)? Defaults to TRUE.
#' @return A list with two named elements:
#'      "full" is a long-format data table with at least five named columns: 'rep.no', integer specifying simulation run; 'bin',
#'      character specifying chronological bin in terms of date range; 'bin.no' integer specifying number of bin, counting from earliest;
#'      'count', numeric giving the number of entities (or total weight) assigned to the given bin in the given simulation run; 'dummy',
#'      giving the number of entities (or total weight) assigned to a bin in the dummy version of a given simulation run. If RoC=TRUE there
#'      will be two more columns: 'RoC.count' and 'RoC.dummy' give the rate of change between this bin and the next for 'count' and 'dummy'
#'      respectively.
#'      "summary" is a second long format data table with four named columns: 'bin', as above; 'V1', the relevant value for the given bin
#'      at a given quantile; 'quantile', the quantile at which V1 is calculated; 'id', character specifying which column from "full" V1 is
#'      based on: e.g. "count", "dummy", "RoC.count", "RoC.dummy".
#' @export
#' @examples
#' date.ranges <- data.table(unit=c(1, 2, 3, 4), Start=c(450, 450, 600, 1000), End=c(700, 800, 650, 1200), frag.count=c(3, 6, 25, 1))
#' x <- date.simulate(date.ranges, weight=date.ranges$frag.count, context.fields="unit", start.date=500, end.date=1500)

date.simulate <- function(data, probs=1, weight=1, ds.fun=sum, real=TRUE, dummy=FALSE, comp.field=NULL, comp.values=NULL, context.fields=c("ID"),
                          quant.list=c(0.025, 0.25, 0.5, 0.75, 0.975), start.date=NULL, end.date=NULL, bin.width=100, reps=100, RoC=NULL, summ=TRUE) {
    #Load required packages
    require(data.table)

    #Tidy up input data
    data <- data.table(cbind(data, weight)) #appends weights to list of date ranges, recycling if necessary (e.g. for uniform weight)

    #Read start and end dates from input data if not specified
    if(is.null(start.date)) {
        start.date <- min(data$Start)
    }
    if(us.null(end.date)) {
        end.date <- max(data$End)
    }
    data <- data[End >= start.date & Start <= end.date]  #drops records outside the date range FROM BOTH SIMULATION SETS

    #Aggregate data
    data <- data[, j=list(weight=sum(as.numeric(weight))), by=c(context.fields, comp.field, "Start", "End")]
    if(length(weight)==1) {data[, weight := 1]}

    # Separate comparison groups, if relevant
    if(!is.null(comp.field)) {
        if(is.null(comp.values)) {comp.values <- unique(data[!get(comp.field)=="" & !is.na(get(comp.field)), get(comp.field)])}
        data <- data[get(comp.field) %in% comp.values]
        x <- as.formula(paste(paste(context.fields, collapse="+"), "+Start+End~", comp.field, sep=""))
        data <- dcast.data.table(data, formula=x, value.var="weight", fill=0)
    }

    #Reset bin.width based on probs, if necessary
    if(is.vector(probs)==TRUE & length(probs)>1) {bin.width <- (end.date - start.date) / length(probs)}  #if probs supplied, use to set bin.widths
    if(sum(class(probs)=="data.frame")==1) {bin.width <- (end.date - start.date) / nrow(probs)}  #likewise if supplied as data.frame/data.table

    #Set up breaks and labels
    breaks <- seq(start.date, end.date, bin.width) #sets breaks
    labels <- numeric(0)
    for(i in 1:(length(breaks) - 1)) {
        labels[i] <- paste(breaks[i], breaks[i + 1], sep="-") #sets bin labels
    }
    probs <- cbind(1:length(labels), probs) #recycles probs to length of labels, if necessary

    #Set up for simulation
    rep.no <- rep(1:reps, each = nrow(data))
    data <- cbind(rep.no, data) #recycles input data 'reps' times to provide frame for simulation
    results <- data.table(rep.no=rep(1:reps, each=length(labels)), bin.no=rep(1:length(labels), reps), bin=rep(labels, reps))

    #Simulate real data
    if(real==TRUE) {

        #Perform simulation
        data[, sim := {x <- runif(nrow(data)); (x * (data$End - data$Start)) + data$Start}] #simulates a date for each row
        data[, bin := cut(sim, breaks, labels=labels)] #finds the relevant bin for each simulated date

        #Aggregate by bin
        if(is.null(comp.field)) {
            real <- data[is.na(bin)==FALSE, j=list(count=ds.fun(as.numeric(weight))), by=list(rep.no, bin)] #sums weights by bin and rep number
        } else {
            real <- data[is.na(bin)==FALSE, lapply(.SD, ds.fun), by=list(rep.no, bin), .SDcols=as.character(comp.values)] #sums weights by bin and rep number
            setnames(real, old=as.character(comp.values), new=paste(comp.values, ".count", sep=""))
        }

        #Merge with frame
        results <- merge(results, real, by=c("rep.no", "bin"), all=TRUE)
    }

    #Simulate dummy set, if required
    if(dummy==TRUE) {

        #Perform simulation
        data[, dummy.bin := sample(labels, size=nrow(data), replace=TRUE, prob=probs[, 2])]

        #Aggregate by bin
        if(is.null(comp.field)) {
            dummy <- data[, j=list(dummy=ds.fun(as.numeric(weight))), by=list(rep.no, dummy.bin)] #sums weights by bin and rep number
        } else {
            dummy <- data[is.na(bin)==FALSE, lapply(.SD, ds.fun), by=list(rep.no, dummy.bin), .SDcols=as.character(comp.values)] #sums weights by bin and rep number
            setnames(dummy, old=as.character(comp.values), new=paste(comp.values, ".dummy", sep=""))
        }

        #Merge with frame
        setnames(dummy, "dummy.bin", "bin")
        results <- merge(results, dummy, by=c("rep.no", "bin"), all=TRUE)
    }

    #Replace any missing values with 0
    results[is.na(results)] <- 0
    results <- results[order(rep.no, bin.no)]

    #Calculate rate of change variables
    if(!is.null(RoC)) {results <- RoC.fun(results, type=RoC)}

    #Create summary dataset, if required and return results
    if(summ==TRUE) {
        summary <- sim.summ(results)
        results <- list(full=results, summary=summary)
    }

    #Return results
    results
}
