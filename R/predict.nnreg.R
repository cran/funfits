"predict.nnreg"<-
function(out, x = out$x, model = NA, derivative = 0, type = "full")
{
	if(is.na(model)) {
		model <- out$best.model
	}
	if(model == out$best.model & out$nfits > 1) {
		cat("Note: the nnreg output data object fit has more than one model",
			fill = T)
		cat(" just the best one is being evaluated.", fill = T)
	}
	fit <- out$model[[model]]
	predict(fit, x, derivative, type = type)
}
