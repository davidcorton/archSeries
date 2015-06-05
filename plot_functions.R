#Function to set up axes

axis.setup <- function(results, field.list=NULL, value.field="V1", lab.sp=1, main="", ylab="Estimated frequency density", type="full") {
    if(class(results)[1]=="list") {results <- results[[1]]}
    if(is.null(field.list)==TRUE) {field.list <- colnames(results)[4:ncol(results)]}
    minmaxer <- numeric(0)
    if(type=="full") {
        for(i in 1:length(field.list)) {minmaxer <- c(minmaxer, results[,get(field.list[i])])}
        minmaxer <- cbind(minmaxer, unique(results$bin.no))
    } else {
        minmaxer <- as.matrix(cbind(results[id%in%field.list, get(value.field)], 1:length(unique(results$bin))))
    }
    plot(minmaxer[,2], minmaxer[,1], xlab="", xaxt="n", main=main, ylab=ylab, type="n")
    names <- unique(results$bin)
    ticks <- seq(1, length(names),by=lab.sp)
    axis(1, at=ticks, labels=unique(results$bin)[ticks], las=2)
}

#Function to plot full dataset as semi-transparent lines. This is a stand-alone function - calls axis.setup.

lines.chron <- function(results, field.list=NULL, col.list=c("darkred", "darkgreen", "blue", "grey"), opacity=20, lab.sp=1, main="", ylab="Estimated frequency density", add=FALSE, legend=TRUE) {
    if(class(results)[1]=="list") {results <- results[[1]]}
    if(is.null(field.list)==TRUE) {field.list <- colnames(results)[4:ncol(results)]}
    if(add==FALSE) {axis.setup(results, field.list=field.list, lab.sp=lab.sp, main=main, ylab=ylab)}
    plist <- data.frame(field.list, col.list, opacity)[1:length(field.list),]
    a <- col2rgb(plist[,2])
    b <- character()
    for(i in 1:nrow(plist)) {
        b[i] <- rgb(a[1,i], a[2,i], a[3,i], plist[i,3], maxColorValue=255)
    }
    for(i in 1:max(results$rep.no)) {
        for(j in 1:length(field.list)) {
            lines(results[rep.no==i, get(field.list[j])], col=b[j])
        }
    }
    if(legend==TRUE) {with(results, legend("topright", legend=field.list, fill=col.list[1:length(b)], bty="n"))}
}

#Function to plot polygons from summary data

poly.chron <- function(results, field.list=NULL, quant=c(0.025, 0.975), col.list=c("darkred", "darkgreen", "blue", "grey"), opacity=126, value.field="V1", lab.sp=1, main="", ylab="Estimated frequency density", med.line=TRUE, add=FALSE, legend=TRUE) {
    if(class(results)[1]=="list") {results <- results[[2]]}
    if(is.null(field.list)==TRUE) {field.list <- unique(results$id)}
    if(add==FALSE) {axis.setup(results, field.list=field.list, lab.sp=lab.sp, main=main, ylab=ylab, value.field=value.field, type="summary")}
    plist <- data.frame(field.list, col.list, opacity)[1:length(field.list),]
    a <- col2rgb(plist[,2])
    b <- character()
    for(i in 1:nrow(plist)) {
        b[i] <- rgb(a[1,i], a[2,i], a[3,i], plist[i,3], maxColorValue=255)
    }
    x <- c(1:length(unique(results$bin)), length(unique(results$bin)):1)
    for(i in 1:length(field.list)) {
        y <- c(results[id==field.list[i]&quantile==quant[1], get(value.field)], rev(results[id==field.list[i]&quantile==quant[2], get(value.field)]))
        skip <- length(x)+1       
        if(substr(field.list[i],1,3)=="RoC") {skip <- length(unique(results$bin))}
        polygon(x[!x==skip],y[!x==skip],col=b[i])
        if(med.line==TRUE) {with(results[id==field.list[i]&quantile==0.500], lines(1:length(unique(bin)), get(value.field), col=col.list[i], lwd=2))}
    }
    if(legend==TRUE) {with(results, legend("topright", legend=field.list, fill=b, bty="n"))}
}

#Function to plot boxes from full data

box.chron <- function(results, field.list=NULL, col.list=c("darkred", "darkgreen", "blue", "grey"), opacity=255, lab.sp=1, main="", ylab="Estimated frequency density", add=FALSE) {
    if(class(results)[1]=="list") {results <- results[[1]]}
    if(is.null(field.list)==TRUE) {field.list <- colnames(results)[4:ncol(results)]}
    if(add==FALSE) {axis.setup(results, field.list=field.list, lab.sp=lab.sp, main=main, ylab=ylab)}
    plist <- data.frame(field.list, col.list, opacity)[1:length(field.list),]
    a <- col2rgb(plist[,2])
    b <- character()
    for(i in 1:nrow(plist)) {
        b[i] <- rgb(a[1,i], a[2,i], a[3,i], plist[i,3], maxColorValue=255)
    }
    for(i in 1:length(field.list)) {
        with(results, boxplot(get(field.list[i])~bin.no, outline=FALSE, xaxt="n", col=b[i], add=TRUE))
    }
}


