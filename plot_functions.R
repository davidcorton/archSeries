
axis.setup <- function(results, field.list=NULL, lab.sp=1, main, ylab) {
    if(is.null(field.list)==TRUE) {field.list <- colnames(results)[4:ncol(results)]}
    minmaxer <- numeric(0)
    for(i in 1:length(field.list)) {minmaxer <- c(minmaxer, get(field.list))}
    minmaxer <- cbind(minmaxer, unique(bin.no))
    plot(minmaxer[,2], minmaxer[,1], xlab="", xaxt="n", main=main, ylab=ylab, type="n")
    names <- unique(bin)
    ticks <- seq(1, length(names),by=lab.sp)
    axis(1, at=ticks, labels=unique(bin)[ticks], las=2)
}

#Function to plot full dataset as semi-transparent lines. This is a stand-alone function - calls axis.setup.

lines.chron <- function(results, field.list=NULL, col.list=c("darkred", "darkgreen", "blue", "grey"), opacity=20, lab.sp=1, main="", ylab="Estimated frequency density") {
    if(class(results)[1]=="list") {results <- results[[1]]}
    attach(results)
    if(is.null(field.list)==TRUE) {field.list <- colnames(results)[4:ncol(results)]}
    axis.setup(results, field.list=field.list, lab.sp=lab.sp, main=main, ylab=ylab)
    col.list <- col.list[1:length(field.list)]
    x <- col2rgb(col.list)
    for(i in 1:length(col.list)) {
        col.list[i] <- rgb(x[1,i], x[2,i], x[3,i], opacity, maxColorValue=255)
    }
    for(i in 1:max(rep.no)) {
        for(j in 1:length(field.list)) {
            lines(get(field.list[j])[rep.no==i], col=col.list[j])
        }
    }
    detach(results)
}

#Function to plot polygons from summary data

poly.chron <- function(results, field.list=NULL, col.list=c("darkred", "darkgreen", "blue", "grey"), opacity=126, med.line=TRUE)
    if(class(results)[1]=="list") {results <- results[[2]]}
    attach(results)
    if(is.null(field.list)==TRUE) {field.list <- unique(id)}
    x <- c(1:length(unique(bin)), length(unique(bin)):1)
    for(i in 1:length(field.list)) {
        y <- c(results[id==field.list[i]&quantile==quant[1], V1], rev(results[id==field.list[i]&quantile==quant[2], V1]))
        polygon(x,y,col=col.list[i])
    }
    detach(results)
}
