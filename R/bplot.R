"bplot"<-
function(x, ..., xpos = NA, width, label, by, srt = 0, add = F, space = 0.25, 
        sort.names = T, xlab = "", ylab = "")
{
# Draws Horizontal Boxplots. Andreas Krause, Dec 1991.
# Doug Nychka April 28, 1992
        if(is.matrix(x)) data <- data.frame(x)
        if(data.class(x) == "numeric")
                data <- list(x, ...)
        if(is.list(x))
                data <- x
        if(!missing(by)) data <- cat.to.list(unlist(data), by)  
        # at this point the data should be in list form regradless of how the
# pieces were originally passed
        quant <- c(0.050000000000000003, 0.25, 0.5, 0.75, 0.94999999999999996)
        cols <- length(data)
        range.data <- range(as.numeric(unlist(data)), na.rm = T)
        if(is.na(xpos[1])) {
                xpos <- 1:cols
        }
        if(missing(width)) {
                width <- min(diff(sort(xpos))) * space
                if(cols == 1)
                        width <- space
        }
        if(length(width) == 1)
                width <- rep(width, cols)
        if(!add) {
                plot(range(c(xpos - (0.5 * width)/space, xpos + (0.5 * width)/
                        space)), range.data, type = "n", xaxt = "n", ylab = 
                        ylab, xlab = xlab)
        }
        for(i in 1:cols) {
                temp <- data[[i]]
                temp <- temp[!is.na(temp)]
                bb <- quantile(temp, quant)
                mid <- xpos[i]  #  i - 0.5
                low <- mid - width[i] * 0.5
                high <- mid + width[i] * 0.5
                if(length(temp) > 5) {
                        y <- c(bb[1], bb[1], NA, bb[1], bb[2], NA, bb[2], bb[2],
                                bb[4])
                        x <- c(high, low, NA, mid, mid, NA, high, low, low)
                        y <- c(y, bb[4], bb[2], bb[3], bb[3], NA, bb[4], bb[5], 
                                bb[5], bb[5])
                        x <- c(x, high, high, high, low, NA, mid, mid, high, 
                                low)
                        lines(x, y)
                }
                if(length(temp) > 5) {
                        outlier <- temp[(temp < bb[1]) | (temp > bb[5])]
                }
                else outlier <- temp
                olen <- length(outlier)
                if(olen > 0)
                        points(rep(mid, olen), outlier)
        }
        if(missing(label)) {
                if(is.null(names(data)))
                        label <- format(1:cols)
                else label <- names(data)
        }
        if(length(label) > 7)
                srt <- 90
        axis(1, xpos, label, tick = F, srt = srt, adj = 1)
        invisible()
}

