"surface.addreg" <-
function (obj, grid.list = NA, extrap = F, graphics.reset = T, 
        ...) 
{
        old.par <- par("mfrow", "oma")
        if (graphics.reset) 
                on.exit(par(old.par))
        out.p <- predict.surface(obj, grid.list = grid.list, 
                extrap = extrap)
        # was plot:
        surface(out.p, type = "b", graphics.reset = graphics.reset, 
                ...)
        invisible()
}
