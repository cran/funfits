"lle.default"<-
function(jac, lags = NA, nprod = c(5, 10, 20, 40, 80), skip = NA, statevector
	 = F, clean = T, verbose = F)
{
	if(!is.na(skip))
		jac <- as.matrix(jac[, -1 * skip])
	n <- nrow(jac)
	m <- length(nprod)	
	# pad the jacobain marix with zero columns for lags not appearing  in
# the model. ( i.e. some lags have a partial deriviatve of zero)
	if(!is.na(lags[1])) {
		hold <- jac	# fill up a matirx with zeroes
		jac <- matrix(0, ncol = max(lags), nrow = n)	
	#replace nonzero columns with columns of the passed jacobian. 
		jac[, lags] <- hold
	}
	if(!is.na(nprod[1])) {
		temp1 <- matrix(nrow = n, ncol = m)
		temp2 <- matrix(nrow = n, ncol = m)
		temp3 <- matrix(nrow = n, ncol = m)
		dimnames(temp1) <- list(NULL, paste(format(nprod), "steps"))
		dimnames(temp2) <- list(NULL, paste(format(nprod), "steps"))
		for(k in 1:length(nprod)) {
			if(nprod[k] < n)
				temp <- make.lle(jac, nprod[k], statevector = 
				  statevector, clean = clean, verbose = verbose
				  )
			temp1[1:(n - nprod[k] + 1), k] <- temp[, 1]
			temp2[1:(n - nprod[k] + 1), k] <- temp[, 2]
			temp3[1:(n - nprod[k] + 1), k] <- temp[, 3]
		}
	}
	else {
		temp1 <- NA
		temp2 <- NA
		temp3 <- NA
	}
	glb <- make.lle(jac, -1, statevector = statevector, clean = clean, 
		verbose = verbose)[, 2]
	temp <- list(local.svd = temp1, local.qr = temp2, local.11 = temp3, 
		nprod = nprod, glb = glb)
	class(temp) <- "lle"
	temp
}
