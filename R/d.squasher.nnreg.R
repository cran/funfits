"d.squasher.nnreg"<-
function(u)
{
	au <- abs(u)
	su <- ifelse(u < 0, -1, 1)
	N <- (u * (1 + 0.5 * au))
	D <- (2 + au + 0.5 * au * au)
	((1 + au) * D - N * (su + u))/(D * D)
}
