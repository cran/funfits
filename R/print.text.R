"print.text"<-
function(txt)
{
	tf <- tempfile()
	write(txt, file = tf, ncol = 1)
	unix(paste("cat ", tf), output = F)
	invisible()
}
