"as.cov.isocor"<-
function(loc, sigma, grid.list = NA, nx = 50, ny = 50, lim.theta = NA, M = ncol(
        sigma), nugget = F, ...)
{
        if(!nugget) rho <- 1    #
        if(is.na(grid.list)) {
                grid.list <- list(seq(min(loc[, 1]), max(loc[, 1]),  , nx), seq(
                        min(loc[, 2]), max(loc[, 2]),  , ny))
        }
###
        out <- list()
        grid <- make.surface.grid(grid.list)
        temp <- (tps(loc, c(sqrt(diag(sigma))), return.matrices = F))
        tps.sum <- (temp)
        print(tps.sum)
        temp2 <- (predict(temp, grid))^2
        temp2 <- as.surface(grid, temp2)
        dst <- rdist.earth(loc, loc)
        ind <- col(dst) > row(dst)
        dst <- dst[ind]
        temp <- 1/sqrt(diag(sigma))
        sigma <- diag(temp) %*% sigma %*% diag(temp)
        print(sigma[1:4, 1:4])
        y <- sigma[ind]
        if(is.na(lim.theta[1])) {
                t1 <- min(dst)
                t2 <- max(dst)
        }
        else {
                t1 <- lim.theta[1]
                t2 <- lim.theta[2]
        }
        theta <- seq(t1, t2,  , 80)
        ss <- rep(NA, length(theta))
        for(k in 1:length(theta)) {
                temp <- exp( - c(dst)/theta[k])
                if(nugget) {
                        rho <- sum(temp * y)/sum(temp^2)
                }
                ss[k] <- sum((y - rho * temp)^2)
        }
        theta.est <- theta[min(ss) == ss]
        list(loc = loc, var = temp2, ss = ss, theta = theta, theta.est = 
                theta.est, rho = rho, summary = tps.sum, dst = dst, y = y, call
                 = match.call())
}

