"as.surface" <-
function (grid.list, z, order.variables = "xy") 
{
        ##
        # arg 1 can either be the grid.list or the actual surface.grid object
        #
        if (is.null(class(grid.list)) | (class(grid.list) != 
                "surface.grid")) {
                hold <- make.surface.grid(grid.list, info.list = T)
        }
        else {
                hold <- attributes(grid.list)$surface.info
        }
        if (hold$nx * hold$ny != length(z)) 
                stop("Problems\nmatching grid info with your z vector. Check your grid dimensions!")
        if (hold$nvar > 2) {
                main.title <- paste(names(hold$fixed.variables), 
                        " = ", unlist(hold$fixed.variables), 
                        sep = "")
                main.title <- paste(main.title, collapse = " ")
        }
        else {
                main.title <- NULL
        }
        if (order.variables == "xy") {
                out <- list(x = hold$x, y = hold$y, z = matrix(z, 
                        ncol = hold$ny, nrow = hold$nx), xlab = hold$xlab, 
                        ylab = hold$ylab, main = main.title)
                class(out) <- "surface"
        }
        else {
                out <- list(x = hold$y, y = hold$x, z = t(matrix(z, 
                        ncol = hold$ny, nrow = hold$nx)), xlab = hold$ylab, 
                        ylab = hold$xlab, main = main.title)
                class(out) <- "surface"
        }
        out
}
