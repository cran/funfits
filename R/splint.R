"splint" <-
function (x, y, xgrid, derivative = 0) 
{
        # This is a lean and mean cubic spline interpolation routine
        # although this is the same subroutie that will also smooth 
        # data the correct set of job codes with provide a simple
        # and fast interpolation wiothout much overhead
        #
        if (!is.loaded(symbol.For("css"))) {
          stop("Compiled code has not been dynamically loaded")
                # 
#                temp <- dyn.load(paste(FUNFITS.BIN, "css.o", 
#                        sep = "/"), 2)
        }
        if (is.matrix(x)) {
                xgrid <- y
                y <- x[, 2]
                x <- x[, 1]
        }
        # remove duplicate X's
        if (is.list(x)) {
                xgrid <- y
                y <- x$y
                x <- x$x
        }
        ind <- !dup(x)
        x <- x[ind]
        y <- y[ind]
        if ((derivative > 2) | (derivative < 0)) 
                stop("derivative must be 0,1,2")
        if (length(x) != length(y)) 
                stop("Lengths of x and y must match")
        n <- length(x)
        #       subroutine css(h,npoint,x,y,wght,sy,trace,diag,vlam,  
        #    +                  ngrid,xg,yg,job,ideriv,ierr)  
        .Fortran("css", as.double(0), as.integer(n), as.double(x), 
                as.double(y), as.double(rep(0, n)), as.double(rep(1, 
                        n)), as.double(1), as.double(1), as.double(1), 
                as.integer(length(xgrid)), as.double(xgrid), 
                ygrid = as.double(rep(0, length(xgrid))), as.integer(c(2, 
                        3, 0)), as.integer(derivative), as.integer(0))$ygrid
}
