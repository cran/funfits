"remove.file"<-
function(file.name)
{
	touch.file(file.name)
	unix(paste("rm ", file.name))
	invisible()
}
