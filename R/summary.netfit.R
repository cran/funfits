"summary.netfit"<-
function(fit, standardized = F)
{
	d <- fit$d
	k <- fit$k
	theta <- fit$theta
	beta <- theta[1:(k + 1)]
	mu <- theta[(1:k) + k + 1]
	gamma <- matrix(theta[(1:(d * k)) + 2 * k + 1], ncol = d, nrow = k, 
		byrow = T)
	if(standardized == F) {
		gamma <- gamma * matrix(1/fit$xsd, ncol = fit$d, nrow = fit$k, 
			byrow = T)
		temp <- (gamma * matrix(fit$xm, ncol = fit$d, nrow = fit$k, 
			byrow = T)) %*% rep(1, fit$d)
		mu <- mu - c(temp)
		beta <- beta * fit$ysd
		beta[1] <- beta[1] + fit$ym
	}
	list(beta = beta, mu = mu, gamma = gamma, standardized = standardized)
}
