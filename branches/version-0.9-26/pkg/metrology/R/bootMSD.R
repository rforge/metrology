#Calculates parametric bootstrapped median scaled difference for observations x given sd's s
#if s is a scalar function, it is applied to x to obtain an estimate of s

#Define generic
bootMSD <- function(x, ...) {
	UseMethod("bootMSD")
}

#Define variant for class mtR.msd (metRology-msd)

bootMSD.MSD <- function(x, B=3000, probs=c(0.95, 0.99), 
	method=c("rnorm", "lhs"), keep=FALSE, labels=names(x), ...) {
	cat("Using MSD variant\n")
	bootMSD.default(attr(x, "x"), attr(x, "s"), B=B, probs=probs, 
		method=method, keep=keep, labels=names(x), ...)
}

bootMSD.default<-function(x, s=mad , B=3000, probs=c(0.95, 0.99), 
	method=c("rnorm", "lhs"), keep=FALSE, labels=names(x), ...) {

	cat("Using default\n")
 	cat(sprintf("with object of class %s\n", class(x)))
 	
 	method <- match.arg(method)

        #Get standard deviations if not a vector
        ss <- if(is.function(s)) {
                rep(s(x, ...), length(x))
        } else {
                if(length(s) == 1) {
                        rep(s, length(x))
                } else {
                        s
                }
        }
        
        #Get the raw values of msd from the data
        t0 <- c(msd(x, ss)) #c() strips unwanted attributes from msd object
        
	N <- length(x)
	
	#Generate simulated data
	if(method=="rnorm") {
		m <- matrix(rnorm(B*N, mean=0, sd=ss), byrow=TRUE, ncol=N)
	} else if(method=="lhs") {
		require(lhs)
		m <- sweep(qnorm(randomLHS(B, N)), MARGIN=2, STATS=ss, FUN='*')
	} else {
		stop(paste("Method", method, "not implemented"))
	}
	
	#Generate the msd simulation
	t <- t(apply(m, 1, function(x, s) as.vector(msd(x, s)), s=ss))

	#Quantiles (critical values)
	ncrit <- length(na.omit(probs))
	if(ncrit > 0) q <- apply(t, 2, quantile, probs=probs)
	if(ncrit==1) q <- matrix(q, nrow=1)
	
	p <- apply(sweep(t, MARGIN=2, STATS=t0, FUN="-"), 2, function(x) sum(x > 0) )
	p <- p/B
	
	structure(list(msd=t0, labels=labels, probs=probs, critical.values=q, 
			pvals=p, B=B, method=method, t=if(keep) t else NA),
			class="bootMSD")
        
}

print.bootMSD <- function(x, ...) {
	print(c(x$msd), ...)
}

summary.bootMSD <- function(x, p.adjust="none", digits=NULL, ...) {
	p.adjust <- match.arg(p.adjust, p.adjust.methods)
	p.adj <- p.adjust(x$pvals, method=p.adjust)
	structure(c(x[c("msd", "labels", "probs", "critical.values", "pvals")], 
		list(p.adjust=p.adjust, p.adj=p.adj),
		x[c('B', "method")]),
		class="summary.bootMSD")
}

print.summary.bootMSD <- function(object, signif.stars = getOption("show.signif.stars"), 
		signif.legend=signif.stars, ...) 
{
	cat("Median Scaled Difference parametric bootstrap\n")
	cat(sprintf("%s replicates\n", format(object$B)))
	cat(sprintf("Sampling method: %s\n", object$method))
	cat(sprintf("P-value adjustment: %s\n", object$p.adjust))
	df.crit <- as.data.frame(t(object$critical.values), check.names=FALSE)
	names(df.crit) <- paste("Upper", names(df.crit))
	m <- format( cbind(data.frame(MSD=object$msd, row.names=object$labels), 
		     df.crit, "P(>MSD)"=object$p.adj) )
	which.0 <- which(object$pvals < 1/object$B)
	if(length(which.0) > 0) {
		#Recalculate adjusted p '+1' for zeros
		p.plus <- p.adjust( pmax(1 + object$B * object$pvals)/object$B,
					method=object$p.adjust )
		m[["P(>MSD)"]][which.0] <- sprintf(" < %7.1e", p.plus[which.0] )
	}
	if (signif.stars && any( object$p.adj < 0.1 )) {
		Signif <- symnum(object$p.adj, corr = FALSE, na = FALSE, 
		cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
			symbols = c("***", "**", "*", ".", " "))
                m <- cbind(m, " "=format(Signif))
	} else {
	       #Nothing significant so no legend
	       signif.legend <- FALSE
	}
 	print(m)
	if (signif.legend) {
		#Using code borrowed from print.Coefmat ...
		if ((w <- getOption("width")) < nchar(sleg <- attr(Signif, 
		    "legend"))) 
		    sleg <- strwrap(sleg, width = w - 2, prefix = "  ")
		cat("---\nSignif. codes:  ", sleg, sep = "", fill = w + 
		    4 + max(nchar(sleg, "bytes") - nchar(sleg)))
	}
 	invisible(object)
}

plot.bootMSD <- function(x, ...) {
	UseMethod("barplot")
}

barplot.bootMSD <- function(height, ylab="MSD", names.arg=height$labels, 
	crit.vals=FALSE, lty.crit=c(2,1), col.crit=2, lwd.crit=c(1,2), ... ) {
	
	mids <- barplot(height$msd, ylab=ylab, names.arg=names.arg, ...)
	
	if(crit.vals && length(na.omit(height$probs)) > 0) {
		ncrit <- length(height$probs)
		dmid <- diff(mids)[1]
		cw <- 0.98*dmid/2
		lty.crit <- rep(lty.crit, ncrit)
		lwd.crit <- rep(lwd.crit, ncrit)
		col.crit <- rep(col.crit, ncrit)
		for(i in 1:ncrit) segments(mids-cw, height$critical.values[i,], mids+cw, height$critical.values[i,],
			lty=lty.crit[i], lwd=lwd.crit[i], col=col.crit[i], lend=2)
	}
	return(invisible(mids))
}

