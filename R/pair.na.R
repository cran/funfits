"pair.na"<-
function(temp, b)
{
	if(!missing(b)) {
		temp <- cbind(temp, b)
	}
	temp[!(is.na(temp[, 1]) | is.na(temp[, 2])),  ]
}
