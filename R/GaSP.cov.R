"GaSP.cov"<-
function(x1, x2, theta = rep(1, ncol(x1)), p = rep(1, ncol(x1)), C = NA)
{
        if(!is.loaded(symbol.For("gaspbs"))) {
                stop("Compiled code has not been dynamically loaded")
                #temp <- dyn.load(paste(FUNFITS.BIN, "TPS.o", sep = "/"), 2)
                #temp2 <- dyn.load(paste(FUNFITS.BIN, "GASP.o", sep = "/"), 2)
        }
        if(!is.matrix(x1))
                x1 <- as.matrix(x1)
        if(missing(x2))
                x2 <- x1
        if(!is.matrix(x2))
                x2 <- as.matrix(x2)     #
        d <- ncol(x1)
        n1 <- nrow(x1)
        n2 <- nrow(x2)  # scale both X's to reflect theta parameter
        x1 <- t(t(x1)/theta)
        x2 <- t(t(x2)/theta)    #       make.gaspb(x1, x2, p)
        par <- p
        if(is.na(C[1])) {
#
#
# return the full covariance matrix
                temp <- .Fortran("gaspbs",
                        nd = as.integer(d),
                        x1 = as.double(x1),
                        n1 = as.integer(n1),
                        x2 = as.double(x2),
                        n2 = as.integer(n2),
                        par = as.double(par),
                        k = as.double(rep(0, n1 * n2)))
                matrix(temp$k, ncol = n2, nrow = n1)
        }
        else {
###
##
# return the covaraince matrix multiplied by the vector C
#
                .Fortran("multgb",
                        nd = as.integer(d),
                        x1 = as.double(x1),
                        n1 = as.integer(n1),
                        x2 = as.double(x2),
                        n2 = as.integer(n2),
                        par = as.double(par),
                        c = as.double(C),
                        h = as.double(rep(0, n1)),
                        work = as.double(rep(0, n2)))$h
        }
}

