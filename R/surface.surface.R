"surface.surface"<-
function(obj, lab = NA, type = "b", zlab, xlab, ylab, graphics.reset = T, ...)
{
	old.par <- par()
        # error in R 0.61.3 ?
        old.par$fin<-NULL
	if(graphics.reset) {
		on.exit(par(old.par))
	}
	if(is.null(obj$xlab))
		obj$xlab <- "X"
	if(is.null(obj$ylab))
		obj$ylab <- "Y"
	if(missing(zlab)) {
		zlab <- "Z"
	}
	if(!missing(xlab)) {
		obj$xlab <- xlab
	}
	if(!missing(ylab)) {
		obj$xlab <- ylab
	}
	if(is.na(lab) & !is.null(obj$main)) {
		lab <- paste("Fixed Variables:  ", obj$main)
	}
#	if(type == "b")
#		set.panel(2, 1, T)
	if(type == "p" | type == "b") {
#	par(mar = c(3, 0, 0, 0))
#		persp(obj, xlab = obj$xlab, ylab = obj$ylab, zlab = zlab, ...)
          cat("Warning: persp not supported in R, using image instead.\n")
          image(obj$x,obj$y,obj$z, xlab = obj$xlab, ylab = obj$ylab, zlab = zlab, ...)
		if(!is.na(lab)) {
#	mtext(lab, 1, 2)
			title(lab)
		}
	}
	if(type == "c" | type == "b") {
#par(mar = c(3, 0, 0, 0))
          if(type == "b")
#            contour(obj$x,obj$y,obj$z, xlab = obj$xlab, ylab = obj$ylab,add=T, ...)
            contour(obj$x,obj$y,obj$z, add=T, ...)
          else
            contour(obj$x,obj$y,obj$z, ...)
		if(!is.na(lab) & type != "b") {
#		mtext(lab, 1, 1, outer = T)
			title(lab)
		}
	}
	invisible()
}
