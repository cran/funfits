"make.surface.grid" <-
function (grid.list, X, nx = 30, ny = 30, info.list = F, FUN = median) 
{
        #
        # what array tells where in the grid.list the x and y variable are. 
        #
        # initialize the what array
        #
        what <- rep(NA, 2)
        # If a X matrix is supplied fill out the grid.list with all the names
        #  of the variables
        if (!is.list(grid.list)) 
                stop("Must supply a list to describe grid limits")
        #
        #
        # at this point the grid.list has been fixed up so that all components
        # are numeric.
        #   The what array has been filled in if the x and y indicators have
        #   been found  
        #
        # ind variable keeps track of where x and y variables are .
        if (!missing(X)) {
                if (data.class(X) != "data.frame") {
                        names.of.X <- dimnames(X)[[2]]
                        if ((is.null(names.of.X))) {
                                names.of.X <- format(1:ncol(X))
                        }
                }
                else names.of.X <- names(X)
                #
                m <- length(names.of.X)
                #print(names.of.X)
                # add integer names to the grid.list if they are missing. 
                #
                #
                # make sure that any component of the grid.list also appears in the 
                # names of X
                #   default integer names have been added at this point so this should
                #work even for grid.lists and X matrices without names
                if (is.null(names(grid.list))) {
                        if (length(grid.list) < m) 
                                stop(" grid.list must be as long as the number of columns of X!")
                        names(grid.list) <- format(1:length(grid.list))
                }
                test <- match(names(grid.list), names.of.X)
                #
                #
                #   default is to center at vairbales that do not appear in the grid
                # specification
                #
                if (!(all(!is.na(test)))) {
                        print("names in grid.list")
                        print(names(grid.list))
                        print("names for columns of X matrix")
                        print(names.of.X)
                        stop(" some of the\ngrid.list names are not found in the names of the X columns")
                }
                # 
                temp <- as.list(rep("c", m))
                names(temp) <- names.of.X
                # now add the grid list info to this master list 
                #
                #
                temp[names(grid.list)] <- grid.list
                for (k in 1:length(temp)) {
                        test <- temp[[k]]
                        if (length(test) == 1) {
                                if (test == "c") 
                                 temp[[k]] <- FUN(X[, k])
                                if (test == "x") {
                                 temp[[k]] <- seq(min(X[, k]), 
                                  max(X[, k]), , nx)
                                 what[1] <- k
                                }
                                if (test == "y") {
                                 temp[[k]] <- seq(min(X[, k]), 
                                  max(X[, k]), , ny)
                                 what[2] <- k
                                }
                        }
                }
                #
                # Now update the orignal grid.list
                #
                grid.list <- temp
        }
        #
        ind <- unlist(lapply(grid.list, length))
        #  check to make sure that the grid list has only two components that
        #have more than 1 value ( i.e the X and Y grid info)
        #
        if (sum(ind > 1) > 2) {
                stop("Only two components can have more than one\nvalue in the grid list")
        }
        #
        nl <- length(grid.list)
        #  fill in what vector in the case when X has not been passed
        #
        if (is.na(what[1])) {
                what <- (1:nl)[ind > 1]
        }
        x1 <- grid.list[[what[1]]]
        x2 <- grid.list[[what[2]]]
        if (length(x1) == 2) {
                x1 <- seq(min(x1), max(x1), , nx)
        }
        if (length(x2) == 2) {
                x2 <- seq(min(x2), max(x2), , ny)
        }
        nx <- length(x1)
        ny <- length(x2)
        # OK at this point x1 and x2 are the real grids in the right order
        nr <- nx * ny
        # Now fill in the constant levels for the other variables
        # 
        if (!info.list) {
                xg <- matrix(NA, ncol = nl, nrow = nr)
                #
                dimnames(xg) <- list(NULL, names(grid.list))
                # attributes contain enough information to reformat the grid for
                # surface plotting or anything else you might think of
                #
                attr(xg, "format") <- cbind(what, ind[what])
                attr(xg, "surface.info") <- list(x = x1, y = x2, 
                        nx = nx, ny = ny, xlab = names(grid.list)[what[1]], 
                        ylab = names(grid.list)[what[2]], fixed.variables = grid.list[-what], 
                        nvar = nl)
                #
                attr(xg, "grid.list") <- grid.list
                # stuff the set of grid points into the right columns of the matrix 
                #
                #
                xg[, what] <- cbind(rep(x1, ny), rep(x2, rep(nx, 
                        ny)))
                for (k in 1:nl) {
                        if (ind[k] == 1) {
                                xg[, k] <- rep(grid.list[[k]], 
                                 nr)
                        }
                }
                class(xg) <- "surface.grid"
                return(xg)
        }
        else {
                #
                # return surface.info component of attributes
                #
                return(list(x = x1, y = x2, nx = nx, ny = ny, 
                        xlab = names(grid.list)[what[1]], ylab = names(grid.list)[what[2]], 
                        fixed.variables = grid.list[-what], grid.list = grid.list, 
                        nvar = nl))
        }
}
