"qr.q2ty"<-
function(qr, y)
{
	qrqr <- qr$qr
	if(is.null(qrqr))
		stop("First argument should be a qr object")
	qra <- qr$qraux
	rank <- qr$rank
	dq <- dim(qrqr)
	if(is.matrix(y)) {
		y <- as.matrix(y)
		dy <- dim(y)
	}
	else {
		dy <- c(length(y), 1)
		y <- as.matrix(y)
	}
	if(dy[1] != dq[1])
		stop("y and qr$qr should have same number of rows")
	if(!(cmplx <- mode(qra) == "complex"))
		storage.mode(y) <- "double"
	else mode(y) <- "complex"
        # orig:
	#.Fortran(if(!cmplx) "dqrsl1" else "zqrsl1",
	#	qrqr,
	#	as.integer(dq),
	#	qra,
	#	as.integer(rank), 
	#	y,
	#	as.integer(dy[2]),
	#	qy = y,
	#	if(!cmplx) 0 else 0i,
	#	as.integer(1000),
	#	as.integer(1))$qy[(rank + 1):dy[1],  ]
        matrix(.Fortran("dqrqty",
                 as.double(qrqr),
                 as.integer(dq[1]),
                 as.integer(rank),
                 as.double(qra),
                 as.double(y),
                 as.integer(dy[2]),
                 qy = as.double(y)
                 )$qy,dy[1],dy[2])[(rank + 1):dy[1], ]
}
