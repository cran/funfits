"netfit.reformat.greedy"<-
function(model)
{
	nm <- length(model)
	temp <- list(1:nm)
	temp[[1]] <- model[[1]]
	if(nm > 1) {
		for(jj in 2:nm) {
			temp[[jj]] <- add.netfit.model(temp[[jj - 1]], model[[
				jj]])
		}
	}
	temp
}
