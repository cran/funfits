"make.lags"<-
function(x, lags, cov = NA, nobs = 3500)
{
    if(is.null(x)) stop("x should exist!")
	x <- as.matrix(x)
	xd <- ncol(x)
	m <- length(lags)
	N <- min(nobs, nrow(x) - max(lags))
	n <- min(nobs, N)
	if(N > nobs)
		warning(" series lengh truncated to\ndefault length in make.lags"
			)
	start <- max(lags) + 1
	temp <- matrix(0, ncol = xd * (length(lags)), nrow = n)
	for(k in 1:length(lags)) {
		a <- start - lags[k]
		b <- a + n - 1
		temp[, (1:xd) + (k - 1) * xd] <- x[(a:b),  ]
	}
	a <- start
	b <- a + n - 1
	if(xd == 1)
		lab <- format(paste("lag", rep(lags, rep(xd, length(lags))), 
			sep = ""))
	else lab <- format(paste(, rep(1:xd, length(lags)), "lag", rep(lags, 
			rep(xd, length(lags))), sep = ""))
	dimnames(temp) <- list(NULL, lab)
	skip <- NA
	if(!is.na(cov[1])) {
		cov <- as.matrix(cov)
		temp <- cbind(temp, cov[a:b,  ])
		cat(a, b)
		skip <- (1:ncol(cov)) + m * xd
	}
	if(xd == 1)
		y <- c(x[a:b])
	else y <- x[a:b,  ]
	list(x = temp, y = y, nvar = m, lags = lags, skip = skip, start = a, 
		end = b)
}
