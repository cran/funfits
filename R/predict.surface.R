"predict.surface" <-
function (out, grid.list = NA, extrap = F, chull.mask, model = NA, 
        nx=30,ny=30) 
{
        if ((length(grid.list) == 1) | (is.na(grid.list)[1])) {
                grid.list <- as.list(rep("c", ncol(out$x)))
                grid.list[[1]] <- "x"
                grid.list[[2]] <- "y"
                ##print(grid.list)
                temp <- dimnames(out$x)[[2]]
                if (!(is.null(temp))) {
                        if (!(temp[1] == "")) 
                                names(grid.list) <- temp
                }
        }
        if (is.null(out$x)) 
                xg <- make.surface.grid(grid.list, nx=nx, ny=ny)
        else xg <- make.surface.grid(grid.list, X = out$x, nx=nx,ny=ny)
        out2 <- as.surface(xg, predict(out, xg, model = model))
        if (!extrap) {
                if (missing(chull.mask)) {
                        ind <- c(attr(xg, "format")[, 1])
                        chull.mask <- out$x[, ind]
                }
                chull.mask <- unique.matrix(chull.mask)
                mask.temp <- interp(chull.mask[, 1], chull.mask[, 
                        2], rep(1, length(chull.mask[, 1])), 
                        xo = out2$x, yo = out2$y)$z
                out2$z <- ifelse(mask.temp == 1, out2$z, NA)
                out2$zlab <- out$yname
        }
        out2
}
