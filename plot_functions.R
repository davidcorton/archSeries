
axis.setup <- function(results, field.list=NULL, lab.sp=1) {
    minmaxer <- numeric(0)
    for(i in 1:length(field.list)) {minmaxer <- c(minmaxer, get(field.list))}
    minmaxer <- cbind(minmaxer, unique(bin.no))
    plot(minmaxer[,2], minmaxer[,1], xlab="", xaxt="n", main=main, ylab=ylab, type="n")
    names <- unique(bin)
    ticks <- seq(1, length(names),by=lab.sp)
    axis(1, at=ticks, labels=unique(bin)[ticks], las=2)
}

lines.chron <- function(results, field.list=NULL, col.list=c("darkred", "darkgreen", "blue", "grey"), opacity=20, lab.sp=1, main="", ylab="Estimated frequency density") {
    if(class(results)=="list") {results <- results[[1]]}
    attach(results)
    if(is.null(field.list)==TRUE) {field.list <- colnames(results)[4:ncol(results)]}
        #try to move the above line into axis.setup
    axis.setup(results, field.list=field.list)
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

