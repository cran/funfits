"ushell.design"<-
function(n.factors, factor.names.arg = std.factor.names(n.factors), n.shells = 
	1, logx = F, digits = 4)
{
## creates uniform shell designs as described in
## David D. Doehlert (1970). "Uniform Shell Designs," Applied Statistics
##(JRSS, Series C) Vol 19, 231-239
	if(length(logx) == 1) logx <- rep(logx, n.factors)
	if(length(logx) != n.factors) stop(paste("Length of", substitute(logx), 
			"can only be 1 or", substitute(n.factors)))	
	## make up the list of generates from p. 233/239 of the paper
##  gen.list <- list(d2=rbind(c(0,0),c(1,0),c(.5,.86602)))
##  gen.list$d3 <- rbind(cbind(gen.list[[1]],rep(0,3)), 
##		       c(.5,.28868,.81650))
##  gen.list$d4 <- rbind(cbind(gen.list[[2]],rep(0,4)), 
##		       c(.5,.28868,.20413,.79057))
	gen.list <- vector("list", length = 10)
	names(gen.list) <- paste("d", 2:11, sep = "")
	gen.list[[1]] <- rbind(c(0, 0), c(1, 0), c(1/2, sqrt(3/4)))
	for(j in 2:10) {
## j = d-1
		gen.mat <- gen.list[[j - 1]]
		gen.mat <- cbind(gen.mat, rep(0, j + 1))
		new.pt <- gen.mat[j + 1, 1:(j - 1)]
		new.pt <- c(new.pt, 1/sqrt(2 * (j + 1) * j), sqrt((j + 2)/(2 * (
			j + 1))))
		gen.mat <- rbind(gen.mat, new.pt)
		gen.list[[j]] <- gen.mat
	}
	max.factors <- length(gen.list) + 1
	if(n.factors > max.factors)
		stop(paste("\nError: Program will handle only up to", 
			max.factors, "factors right now.\n"))
	gen.mat <- gen.list[[n.factors - 1]]	
	## if there are d factors, there are d+1 rows in the generating
## matrix. When each point is subtracted from every other point, 
## then it creates d new points. So the total number of points
## is (d+1)*d. One point of the design is the origin, and the
## remaining d(d+1) points lie on a sphere of radius 1
	d <- n.factors
	n.gen <- d + 1
	n <- d * (d + 1) + 1
	des.mat <- matrix(NA, nrow = n, ncol = n.factors)
	des.mat[1:n.gen,  ] <- gen.mat
	i <- n.gen + 1
	for(j in 2:n.gen) {
		for(l in 1:n.gen) {
			if(l != j) {
				des.mat[i,  ] <- gen.mat[l,  ] - gen.mat[j,  ]
				i <- i + 1
			}
		}
	}
## if there are multiple shells create them here
	if(n.shells > 1)
		add.shells <- T
	else add.shells <- F
	while(add.shells) {
		n.gen <- nrow(des.mat)
		n <- n.gen + (n.gen - 1)^2
		new.des <- matrix(NA, nrow = n, ncol = d)
		new.des[1:n.gen,  ] <- des.mat
		j <- n.gen + 1	
	## add and subtract the generating points to the existing
## matrix
		for(i in 2:n.gen) {
			for(k in 2:(n.gen)) {
				new.des[j,  ] <- new.des[i,  ] - new.des[k,  ]
				j <- j + 1
			}
		}
		des.mat <- new.des
		des.mat <- des.mat[!dup.matrix(round(des.mat, digits)),  ]	
	## now check to see how many shells there are
		r.vals <- sqrt(apply(des.mat^2, 1, sum))
		r.vals[r.vals <= 10^-6] <- 0
		r.vals <- round(r.vals, digits)
		u.rvals <- sort(unique(r.vals))[-1]
		n.r <- length(u.rvals) - 1	## don't include 0
		if(n.r >= n.shells) add.shells <- F	#      browser()
	}
## now pull off the desired number of shells
	if(n.shells > 1) {
		r.cut <- u.rvals[n.shells]
		des.mat <- des.mat[r.vals <= r.cut,  ]
		r.vals <- u.rvals[1:n.shells]
	}
## now assign the factor names and levels
	if(!is.list(factor.names.arg)) {
		fnames <- as.list(rep(0, n.factors))
		names(fnames) <- factor.names.arg
		for(i in 1:n.factors)
			if(logx[i]) fnames[[i]] <- log(c(0.10000000000000001, 
				  10)) else fnames[[i]] <- c(-1, 1)
	}
	else {
		if(length(factor.names.arg) != n.factors)
			stop(paste(
				"Number of elements in factor.names.arg must equal",
				n.factors))
		if(!all(sapply(factor.names.arg, function(x)
		length(x) == 2)))
			stop(paste("If a list, each component of", 
				"factor.names.arg must have length 2"))
		fnames <- factor.names.arg
		for(i in 1:n.factors)
			if(logx[i]) fnames[[i]] <- log(fnames[[i]])
	}
	des.mat <- des.mat/max(abs(des.mat))
	des <- data.frame(des.mat)
	des <- eval(parse(text = paste("des[order(", paste("des[,", 1:n.factors,
		"]", sep = "", collapse = ","), "),]", sep = "")))
	for(i in 1:n.factors) {
		x <- fnames[[i]]
		y <- des[, i]
		a <- x[1]
		b <- x[2]
		center.x <- 0.5 * (a + b)
		scale.x <- (b - a)/2
		y <- y * scale.x + center.x
		attributes(y) <- list(log = logx[i], scale = scale.x, center = 
			center.x, class = "rsm.factor")
		des[, i] <- y
	}
	dimnames(des) <- list(1:nrow(des), names(fnames))
	class(des) <- c("rsm.design", "design", "data.frame")
	attr(des, "n.shells") <- n.shells
	if(n.shells > 1)
		attr(des, "r.vals") <- r.vals
	des
}
