"as.EOF.cov"<-
function(loc, sigma, grid.list = NA, nx = 50, ny = 50, M = ncol(sigma))
{
#
        if(is.na(grid.list)) {
                grid.list <- list(seq(min(loc[, 1]), max(loc[, 1]),  , nx), seq(
                        min(loc[, 2]), max(loc[, 2]),  , ny))
        }
###
        grid <- make.surface.grid(grid.list)
        hold <- eigen(sigma)
        print(hold$values)
        fits <- as.list(1:M)
        for(k in 1:M) {
                cat(k, " ")
                temp <- predict(tps(loc, c(hold$vectors[, k]), 0, 
                        return.matrices = F), grid)
                cat(k, " ")
                fits[[k]] <- as.surface(grid, temp)
        }
#        Now make up weights 
        temp.e <- hold$values[1:M]
        call <- match.call()
        list(loc = loc, delta = temp.e, fits = fits, M = M, call = match.call()
                )
}

