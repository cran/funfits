"COR"<-
function(dat)
{
	m <- ncol(dat)
	hold <- c(1, 1)
	temp <- matrix(NA, nrow = m, ncol = m)
	temp2 <- matrix(NA, nrow = m, ncol = m)
	ntemp <- matrix(0, nrow = m, ncol = m)
	hold <- NULL
	hold.stat <- stats(dat)
	for(k in 1:m) {
		for(j in 1:k) {
			hold <- pair.na(dat[, j], dat[, k])
			n <- nrow(hold)
			ntemp[j, k] <- ntemp[k, j] <- n
			if(nrow(hold) > 0) {
				cv <- cor(hold[, 1], hold[, 2])
				temp[j, k] <- temp[k, j] <- cv
			}
		}
		cat(k, " ", fill = T)
	}
	list(cor = temp, N = ntemp, sd = c(hold.stat[3,  ]), mean = c(hold.stat[
		2,  ]))
}
