"INTERP"<-
function(x, y, z, ...)
{
        ind <- !(dup(x) & dup(y))
        x <- x[ind]
        y <- y[ind]
        z <- z[ind]
        interp(x, y, z, ...)
}

