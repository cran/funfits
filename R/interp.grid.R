"interp.grid" <-
function (loc, grid) 
{
        xg <- unique(great.lakes.rom.grid[, 1])
        yg <- unique(great.lakes.rom.grid[, 2])
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
        ex <- lx - lx1
        # simple linear interpolation of the rom grid
        ey <- ly - ly1
        #
        #
        (grid[cbind(lx1, ly1)] * (1 - ex) * (1 - ey) + grid[cbind(lx1 + 
                1, ly1)] * (ex) * (1 - ey) + grid[cbind(lx1, 
                ly1 + 1)] * (1 - ex) * (ey) + grid[cbind(lx1 + 
                1, ly1 + 1)] * ex * ey)
}
