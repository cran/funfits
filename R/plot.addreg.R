"plot.addreg" <-
function (out, ...) 
{
        old.par <- par("mfrow", "oma")
        on.exit(par(old.par))
        m <- ncol(out$predicted.comp)
        set.panel(min(m/2 + 1.5, 3), 2)
        plot(out$fitted.values, out$residuals, ...)
        hist(out$residuals)
        for (k in (1:m)) {
                if (out$lam[k] != 0) 
                        ix <- order(out$x[, k])
                ytemp <- out$y - out$fitted.values + out$predicted.comp[, 
                        k]
                plot(out$x[ix, k], ytemp[ix], xlab = paste("variable", 
                        k), ylab = "estimated ridge function", 
                        pch = ".", ...)
                lines(out$x[ix, k], out$predicted.comp[ix, k])
        }
}
