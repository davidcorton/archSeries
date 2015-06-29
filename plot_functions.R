#Function to set up axes

axis.setup <- function(results, field.list=NULL, value.field="V1", lab.sp=1, main="", ylab="Estimated frequency density", ylim=NULL, type=1) {
    if(class(results)[1]=="list") {results <- results[[type]]}
    minmaxer <- numeric(0)
    if(type==1) {
        if(is.null(field.list)) {field.list <- colnames(results)[!colnames(results)%in%c("bin","bin.no","rep.no","x","y")]}
        for(i in 1:length(field.list)) {minmaxer <- c(minmaxer, results[,get(field.list[i])])}
        minmaxer <- cbind(minmaxer, unique(results$bin.no))
    } else {
        if(is.null(field.list)) {field.list <- unique(results$id)}
        minmaxer <- as.matrix(cbind(results[id%in%field.list, get(value.field)], 1:length(unique(results$bin))))
    }
    if(!is.null(ylim)) {ylim <- c(0,ylim)}
    plot(minmaxer[,2], minmaxer[,1], xlab="", xaxt="n", main=main, ylab=ylab, type="n", ylim=ylim)
    names <- unique(results$bin)
    ticks <- seq(1, length(names),by=lab.sp)
    axis(1, at=ticks, labels=names[ticks], las=2)
}

#Function to plot full dataset as semi-transparent lines. This is a stand-alone function - calls axis.setup.

lines.chron <- function(results, field.list=NULL, col.list=c("darkred", "darkgreen", "blue", "grey", "goldenrod"), opacity=20, lab.sp=1,
                        main="", ylab="Estimated frequency density", ylim=NULL, small.n=NULL, small.n.op=126, add=FALSE, legend=TRUE) {
    boxes <- NULL
    if(class(results)[1]=="list") {
        if(length(results)==3) {boxes <- results$small.n}
        results <- results[[1]]
    }
    if(is.null(field.list)==TRUE) {field.list <- colnames(results)[!colnames(results)%in%c("bin","bin.no", "rep.no","x","y")]}
    if(add==FALSE) {axis.setup(results, field.list=field.list, lab.sp=lab.sp, main=main, ylab=ylab, ylim=ylim)}
    if(!is.null(small.n)&!is.null(boxes)) {grey.zones(boxes, small.n, small.n.op, ylim)} #Sets up boxes to highlight small n
    plist <- data.frame(field.list, col.list[1:length(field.list)], opacity)
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

poly.chron <- function(results, field.list=NULL, quant=c(0.025, 0.975), col.list=c("darkred", "darkgreen", "blue", "grey", "goldenrod"),
                       opacity=126, value.field="V1", lab.sp=1, main="", ylab="Estimated frequency density", med.line=TRUE, ylim=NULL, 
                       small.n=NULL, small.n.op=126, add=FALSE, legend=TRUE) {
    boxes <- NULL
    if(class(results)[1]=="list") {
        if(length(results)==3) {boxes <- results$small.n}
        results <- results[[2]]
    }
    if(is.null(field.list)==TRUE) {field.list <- unique(results$id)}
    if(add==FALSE) {axis.setup(results, field.list=field.list, lab.sp=lab.sp, main=main, ylab=ylab, ylim=ylim, value.field=value.field, type=2)}
    if(!is.null(small.n)&!is.null(boxes)) {grey.zones(boxes, small.n, small.n.op, ylim)} #Sets up boxes to highlight small n
    plist <- data.frame(field.list, col.list[1:length(field.list)], opacity)
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
        if(med.line==TRUE) {
            with(results[id==field.list[i]&quantile==0.500], lines(1:length(unique(bin)), get(value.field), col=col.list[i], lwd=3))
            with(results[id==field.list[i]&quantile==0.500], points(1:length(unique(bin)), get(value.field), col=col.list[i], pch=19))
        }
    }
    if(legend==TRUE) {with(results, legend("topright", legend=field.list, fill=b, bty="n"))}
}

#Function to plot boxes from full data

box.chron <- function(results, field.list=NULL, col.list=c("darkred", "darkgreen", "blue", "grey"), opacity=255, lab.sp=1, main="",
                     ylab="Estimated frequency density", ylim=NULL, small.n=NULL, small.n.op=126, add=FALSE) {
    boxes <- NULL
    if(class(results)[1]=="list") {
        if(length(results)==3) {boxes <- results$small.n}
        results <- results[[1]]
    }
    if(is.null(field.list)==TRUE) {field.list <- colnames(results)[!colnames(results)%in%c("bin","bin.no", "rep.no","x","y")]}
    if(add==FALSE) {axis.setup(results, field.list=field.list, lab.sp=lab.sp, main=main, ylab=ylab, ylim=ylim)}
    if(!is.null(small.n)&!is.null(boxes)) {grey.zones(boxes, small.n, small.n.op, ylim)} #Sets up boxes to highlight small n
    plist <- data.frame(field.list, col.list[1:length(field.list)], opacity)
    a <- col2rgb(plist[,2])
    b <- character()
    for(i in 1:nrow(plist)) {
        b[i] <- rgb(a[1,i], a[2,i], a[3,i], plist[i,3], maxColorValue=255)
    }
    for(i in 1:length(field.list)) {
        with(results, boxplot(get(field.list[i])~bin.no, outline=FALSE, xaxt="n", col=b[i], add=TRUE))
    }
}

#Function for plotting aorist output (this is really just barplot with some tweaks
#added, e.g. to make bars line up with points plotted by other functions)

aorist.plot <- function(aorist, col="grey", opacity=255, lab.sp=1, add=FALSE) {
    col <- col2rgb(col)
    col <- rgb(col[1,], col[2,], col[3,], opacity, maxColorValue=255)
    with(aorist, barplot(aorist, space=c(0.5,rep(0, length(bin)-1)), col=col, add=add))
    if(add==FALSE) {
        ticks <- seq(1, nrow(aorist),by=lab.sp)
        axis(1, at=ticks, labels=aorist$bin[ticks], las=2)
    }
}

#Function to plot small-n polygons as output by cpue

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

