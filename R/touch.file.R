"touch.file"<-
function(file.name)
{
	unix(paste("touch ", file.name))
}
