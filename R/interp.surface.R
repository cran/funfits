"interp.surface" <-
function (obj, loc) 
{
        # if(missing(loc))
        # loc <- make.surface.grid(list(range(obj$x), range(obj$y)))
        xg <- unique(obj$x)
        yg <- unique(obj$y)
        nx <- length(xg)
        ny <- length(yg)
        xa <- min(xg)
        xb <- max(xg)
        xr <- xb - xa
        ya <- min(yg)
        yb <- max(yg)
        yr <- yb - ya
        lx <- ((nx - 1) * (loc[, 1] - xa))/xr + 1
        ly <- ((ny - 1) * (loc[, 2] - ya))/yr + 1
        lx1 <- ifelse(lx == nx, nx - 1, trunc(lx))
        ly1 <- ifelse(ly == ny, ny - 1, trunc(ly))
        lx1 <- ifelse(lx1 < 1 | lx1 > nx, NA, lx1)
        ly1 <- ifelse(ly1 < 1 | ly1 > ny, NA, ly1)
        ex <- lx - lx1
        # simple linear interpolation of the rom obj$z
        ey <- ly - ly1
        #
        #
        # if(class(loc) == "surface.grid")
        temp <- (obj$z[cbind(lx1, ly1)] * (1 - ex) * (1 - ey) + 
                obj$z[cbind(lx1 + 1, ly1)] * (ex) * (1 - ey) + 
                obj$z[cbind(lx1, ly1 + 1)] * (1 - ex) * (ey) + 
                obj$z[cbind(lx1 + 1, ly1 + 1)] * ex * ey)
        temp
}
