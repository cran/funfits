"summary.nlar"<-
function(obj, ...)
{
	cat("Call:", fill = T)
	print(obj$call)
	cat("summary of fitted model(s)", fill = T)
	summary(obj$fit, ...)
}
