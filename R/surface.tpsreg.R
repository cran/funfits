"surface.tpsreg"<-
function(tpsreg.object, vary, granularity = 25, transform, x.scale, ..., main)
{
	fnames <- dimnames(tpsreg.object$x)[[2]]	
	## Make sure the vary parameter is usable
	nfactors <- length(fnames)
	if(missing(vary)) {
		vary <- list()
		for(i in 1:nfactors)
			if(i <= 2) eval(parse(text = paste("vary$", fnames[i], 
				  " <- 'v'", sep = ""))) else eval(parse(text
				   = paste("vary$", fnames[i], " <- 'c'", sep
				   = "")))
	}
	if(length(vary) > nfactors)
		stop("Error: vary has too many components.")
	nmv <- names(vary)
	nvary <- 0
	for(nm in fnames) {
		if(is.na(charmatch(nm, nmv)))
			eval(parse(text = paste("vary$", nm, " <- 'c'", sep = 
				"")))
		eval(parse(text = paste("if(vary$", nm, 
			" == 'v') nvary <- nvary+1", sep = "")))
	}
	if(nvary != 2) stop("Error: Must have exactly 2 factors varying")	
	## ready to call grid routine 
	vnames <- names(vary)
	index <- charmatch(fnames, vnames)
	index.v <- rep(F, nfactors)
	mf <- tpsreg.object$x
	ranges <- t(apply(mf, 2, range))
	holdat <- rep(0, nfactors)
	for(i in 1:nfactors) {
		if(vary[[index[i]]] == "v")
			index.v[i] <- T
		else if(vary[[index[i]]] == "c") {
			if(missing(x.scale))
				holdat[i] <- mean(ranges[i,  ])
			else holdat[i] <- mean(x.scale[i])
		}
		else {
			if(missing(x.scale))
				holdat[i] <- vary[[index[i]]]
			else {
				a <- min(x.scale[, i])
				b <- max(x.scale[, i])
				holdat[i] <- (vary[[index[i]]] - a)/(b - a)
			}
		}
	}
	varies <- (1:nfactors)[index.v]
	if((length(varies) > nfactors) && !all(ii <- (ranges[ - varies, 1] <= 
		holdat[ - varies] && ranges[ - varies, 2] >= holdat[ - varies])
		))
		stop(paste("Fixed variables must be held at a value", 
			"in the range of the model data.", paste("Variable", 
			fnames[ - varies][!ii], "is outside its range", 
			collapse = "\n")))
	xx <- seq(ranges[varies[1], 1], ranges[varies[1], 2], length = 
		granularity)
	yy <- seq(ranges[varies[2], 1], ranges[varies[2], 2], length = 
		granularity)
	xm <- rep(xx, rep(granularity, granularity))
	ym <- rep(yy, granularity)
	x.pred <- matrix(NA, ncol = nfactors, nrow = granularity * granularity)
	x.pred[, varies[1]] <- xm
	x.pred[, varies[2]] <- ym
	nvaries <- (1:nfactors)[!index.v]
        # error in R !?
        if(length(nvaries)!=0)
          for(i in 1:length(nvaries))
            x.pred[, nvaries[i]] <- rep(holdat[nvaries[i]], granularity * 
                                        granularity)
	dimnames(x.pred) <- list(1:(granularity * granularity), fnames)
	ypred.mat <- predict(tpsreg.object, x.pred)
	z.mat <- matrix(as.numeric(ypred.mat), nrow = granularity, ncol = 
		granularity, byrow = T)
	if(!missing(transform))
		eval(parse(text = paste("z.mat <- ", paste(substitute(transform
			)), "(z.mat)", sep = "")))
	if(!missing(x.scale)) {
		a <- min(x.scale[, varies[1]])
		b <- max(x.scale[, varies[1]])
		xx <- xx * (b - a) + a
		a <- min(x.scale[, varies[2]])
		b <- max(x.scale[, varies[2]])
		yy <- yy * (b - a) + a
	}
	result <- surface.default(xx, yy, z.mat, ..., xlab = vnames[varies[1]], 
		ylab = vnames[varies[2]])
	if(missing(main)) {
		main <- paste("Surface plot for", tpsreg.object$yname)
		if(length(nvaries) > 0) {
			for(i in 1:length(nvaries)) {
				if(missing(x.scale))
				  holdx <- holdat[nvaries[i]]
				else {
				  a <- min(x.scale[, nvaries[i]])
				  b <- max(x.scale[, nvaries[i]])
				  holdx <- holdat[nvaries[i]] * (b - a) + a
				}
				if(i == 1)
				  main <- paste(main, "\nHolding ", fnames[
				    nvaries[i]], "=", holdx, sep = "")
				else main <- paste(main, ",", fnames[nvaries[i]
				    ], "=", holdx, sep = "")
			}
		}
	}
	title(main = main)
	invisible(result)
}
