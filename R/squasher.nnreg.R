"squasher.nnreg"<-
function(u)
{
	au <- abs(u)
	(u * (1 + 0.5 * au))/(2 + au + 0.5 * au * au)
}
