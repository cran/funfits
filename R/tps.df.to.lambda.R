"tps.df.to.lambda"<-
function(df, D, guess = 1)
{
	if(is.na(df))
		return(NA)
	if(df < sum(D == 0)) {
		warning("df too small to match with a lambda value")
		return(NA)
	}
	if(df > length(D)) {
		warning(" df too large to match a lambda value")
		return(NA)
	}
	l1 <- guess	########## find upper lambda
	for(k in 1:8) {
		tr <- sum(1/(1 + l1 * D))
		if(tr <= df)
			break
		l1 <- l1 * 2
	}
########## find lower lambda
##########
	l2 <- guess
	for(k in 1:8) {
		tr <- sum(1/(1 + l2 * D))
		if(tr >= df)
			break
		l2 <- l2/2
	}
	info <- list(D = D, df = df, N = length(D))
	out <- bisection.search(log(l1), log(l2), tps.fdf, tol = 0.0001, 
		f.extra = info)$x
 + exp(out)
}
