"surface" <-
function (x, grid.list = NA, extrap = F, xlab = NULL, ylab = NULL, 
        zlab = NULL, main = NULL, levels = NULL, zlim = NULL, nx=30, ny=30,
        ...) 
UseMethod("surface")
