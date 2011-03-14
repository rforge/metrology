# Mandel-Paule algorithm
# Based on code by S Cowen
#

mpaule <- function(x, ..., tol=.Machine$double.eps^0.25, maxiter=25) {
	mandel.paule(x, ..., tol=.Machine$double.eps^0.25, maxiter=25)
}

mandel.paule <- function(x, ..., tol=.Machine$double.eps^0.25, maxiter=25) {
	UseMethod("mandel.paule")
}

mandel.paule.default <- function(x, u=NULL, n=NULL, groups=NULL, tol=.Machine$double.eps^0.25, maxiter=25)		
{
#If n present, u is interpreted as sd of n observations
#	x <- input.data[,1]
#	u <- input.data[,2]
#	n <- input.data[,3] + 1
	
	
	count<-function(x) sum(!is.na(x))
	
	if(!is.null(groups)) {
		groups <- factor(groups)
		x.i <- as.vector(tapply(x, groups, mean, na.rm=TRUE))
		u.i <- as.vector(tapply(x, groups, sd, na.rm=TRUE))
		u.i <- u.i/sqrt(as.vector(tapply(x, groups, count)))
	} else {
		x.i <- x
		if(is.null(n)) 
			u.i <- u
		else
			u.i <- u/sqrt(rep(n, length(x.i)))  #Recycles n

	}
	v <- var(x.i)

	iter <- 0
	dv <- v
	while(iter < maxiter && abs(dv) > tol )	
	{
		iter<- iter + 1
		wt <- 1 / ( u.i^2 + v )
		cons.mean <- sum(wt * x.i) / sum(wt)
		F <- sum( wt * (x.i - cons.mean)^2 ) - (length(x.i) - 1)	
		dv <- F / sum( wt^2 * (x.i - cons.mean)^2 )
		v <- v + dv
		if(v < 0) {
			v <- 0
			dv<-0
		}
		
	}
	
	if(abs(dv) >= tol) warning("Maximum iterations reached; M-P may not have converged", call.=TRUE)

	rv <- .construct.loc.est( x=cons.mean, u=1/sqrt(sum(wt)), df=NA, 
		xi=x.i, ui=u.i, u.eff=sqrt(v+u.i^2), w=rep(1, length(x.i)), method="Mandel-Paule", method.details=list(iter=iter, var.between=v))
	
	return(rv)
}
