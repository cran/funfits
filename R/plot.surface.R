"plot.surface"<-
function(obj, main = NULL, type = "b", zlab = NULL, xlab = NULL, ylab = NULL, 
	levels = NULL, zlim = NULL, graphics.reset = T, ...)
{
	old.par <- par()
        # error in R 0.61.3
        old.par$fin <- NULL
	if(graphics.reset)
		on.exit(par(old.par))
	if(is.null(xlab)) {
		if(is.null(obj$xlab))
			xlab <- "X"
		else xlab <- obj$xlab
	}
	if(is.null(ylab)) {
		if(is.null(obj$ylab))
			ylab <- "Y"
		else ylab <- obj$ylab
	}
	if(is.null(zlab)) {
		if(is.null(obj$zlab))
			zlab <- "Z"
		else zlab <- obj$zlab
	}
	if(is.null(main))
		if(!is.null(obj$main))
			main <- obj$main
#	if(type == "b") # makes no sense without persp
#		set.panel(2, 1, T)
	if(type == "p" | type == "b") {
##	par(mar = c(3, 0, 0, 0))
#		if(is.null(zlim)) persp(obj, xlab = xlab, ylab = ylab, zlab = 
#				zlab, ...)
#                else persp(obj, xlab = xlab, ylab = 
#				ylab, zlab = zlab, zlim = zlim, ...)
          cat("Warning: persp not supported in R, using image instead.\n")
		if(is.null(zlim))
                  image(obj$x,obj$y,obj$z, xlab = xlab, ylab = ylab,
                        zlab = zlab, ...)
                else
                  image(obj$x, obj$y, obj$z, xlab = xlab, ylab = ylab,
                        zlab = zlab, zlim = zlim, ...)
		if(!is.null(main))
                  title(main)
	}
	if(type == "c" | type == "b") {
##par(mar = c(3, 0, 0, 0))
		if(is.null(levels)) levels <- pretty(range(obj$z,na.rm=T), 5)
                if(type=="b")
#                  contour(obj$x,obj$y,obj$z, xlab = xlab, ylab = ylab, levels = levels,add=T, ...)
                  contour(obj$x,obj$y,obj$z, levels = levels,add=T, ...)
                else
                  contour(obj$x,obj$y,obj$z, xlab = xlab, ylab = ylab, levels = levels, ...)
		if((!is.null(main)) & type != "b")
			title(main)
	}
	invisible()
}
