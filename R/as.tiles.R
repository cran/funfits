"as.tiles"<-
function(grid.list, FUN)
{
        grid <- make.surface.grid(grid.list)
        g.info <- make.surface.grid(grid.list, info.list = T)   
        # grid has the nu*nv in the right order for reformatting as a matrix
# persp plotting. ( X varies fastest Y held fixed)
#
        xyz <- FUN(grid)        
        # Ok now is the fun  part. Need to take the xyz values and repeat to
#make the tiles. 
#
#  reformat as a matrix to make the indexing easier. 
#
        nu <- g.info$nx
        nv <- g.info$ny #
#     order for constructing tile is     1 2
#                                        4 3
#
        temp <- matrix(xyz[, 1], ncol = nv, nrow = nu)  
        # Now X is a vector in groups of 4 coordinates representing the x
#coordinates of teh tiel corners. The NA is there to end the polygon
#
#
        X <- (rbind(c(temp[1:(nu - 1), 1:(nv - 1)]), c(temp[2:(nu), 1:(nv - 1)]
                ), c(temp[2:(nu), 2:(nv)]), c(temp[1:(nu - 1), 2:nv]))) 
        #  Now do the same thing for the y  and z
        temp <- matrix(xyz[, 2], ncol = nv, nrow = nu)
        Y <- (rbind(c(temp[1:(nu - 1), 1:(nv - 1)]), c(temp[2:(nu), 1:(nv - 1)]
                ), c(temp[2:(nu), 2:(nv)]), c(temp[1:(nu - 1), 2:nv]))) #
        temp <- matrix(xyz[, 3], ncol = nv, nrow = nu)
        Z <- (rbind(c(temp[1:(nu - 1), 1:(nv - 1)]), c(temp[2:(nu), 1:(nv - 1)]
                ), c(temp[2:(nu), 2:(nv)]), c(temp[1:(nu - 1), 2:(nv)])))       
#
        one <- rep(0.25, 4)
        list(X = t(X)[, c(1:4, 1)], 
             Y = t(Y)[, c(1:4, 1)], 
             Z = t(Z)[, c(1:4, 1)], 
             uv = grid, 
             center = cbind(t(X) %*% one, t(Y) %*% one, t(Z) %*% one))
}

