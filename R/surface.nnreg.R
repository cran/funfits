"surface.nnreg"<-
function(obj, grid.list = NA, extrap = F, graphics.reset = T, xlab = NULL, ylab
	 = NULL, main = NULL, zlab = NULL, zlim = NULL, levels = NULL, model = 
	NA, ...)
{
## modified so that you can give main, and ylab as arguments
## in ... and have them passed correctly
	old.par <- par("mfrow", "oma")
	if(graphics.reset)
		on.exit(par(old.par))
	out.p <- predict.surface(obj, grid.list = grid.list, extrap = extrap, 
		model = model)
	if(!is.null(ylab))
		out.p$ylab <- ylab
	if(!is.null(xlab))
		out.p$xlab <- xlab
	if(!is.null(zlab))
		out.p$zlab <- zlab
	if(!is.null(main)) out.p$main <- main	##    else
##      out.p$main <- NULL
	plot(out.p, type = "b", graphics.reset = graphics.reset, levels = 
		levels, zlim = zlim, ...)
	invisible()
}
