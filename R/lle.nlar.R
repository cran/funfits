"lle.nlar"<-
function(obj, model = NA, nprod = c(5, 10, 20, 40, 80), verbose = F, clean = T, 
	...)
{
	statevector <- F
	if(is.na(model)) {
		model <- obj$fit$best.model
	}
	jac <- predict(obj$fit, derivative = 1, model = model, ...)	#
# omit columns  of Jacobian that are not state variables. 	
#
	lle.default(jac, lags = obj$lags, skip = obj$skip, nprod = nprod, clean
		 = clean, verbose = verbose)
}
