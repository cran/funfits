"std.factor.names"<-
function(n)
{
	ll <- LETTERS[ - c(6, 20)]
	if((m <- ceiling(n/length(ll))) > 1)
		ll <- c(ll, outer(ll, ll[1:(m - 1)], function(x, y)
		paste(x, y, sep = "")))
	ll[1:n]
}
 
"unix" <-
function (command, input, output.to.S = T) 
{
        if (!missing(input)) {
                file <- tempfile("unix")
                on.exit(unlink(file))
                cat(input, file = file, sep = "\n")
                command <- paste("<", file, command)
        }
        system(command, output.to.S)
}
