"predict.addreg" <-
function (out, xnew = out$x, model = NA) 
{
        temp <- rep(0, nrow(xnew))
        m <- ncol(out$x)
        for (k in 1:m) {
                temp <- temp + splint(out$x[, k], out$predicted.comp[, 
                        k], xnew[, k])
        }
        temp
}
