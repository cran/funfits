"print.nnreg"<-
function(fit)
{
#	cat(" Call: ", fit$call,sep="\n", fill = T)	
	cat(" Call: ")
        print(fit$call)
	#	cat("Dimension of surface:", fit$d, fill = T)
#out <- cbind(fit$k, fit$np)
#	dimnames(out) <- list(rep(" ", length(fit$k)), c("# units", 
#		"#  parameters   "))
#print(out)
	cat("Detailed output from FORTRAN program", " \n \n", fill = T)
	print.text(fit$summary)
	invisible()
}
