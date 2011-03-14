#uncertMC object
print.uncertMC<-function(uncertMC, digits=NULL, right=FALSE, ..., simplify=TRUE, minimise=FALSE){
	maxwidth<-12L
	cat("\nUncertainty evaluation\n\n")

	cat("Call:\n  ",deparse(uncertMC$call), sep="")
	cat("\n\n")
	cat("Expression: ")
	if(class(uncertMC$expr)=="formula" ) {
		cat(paste(uncertMC$expr, collapse=""))
	} else if(is.function(uncertMC$expr)) {
		cat( deparse(uncertMC$expr)[1] )
	} else if(is.expression(uncertMC$expr)) {
		cat( deparse(uncertMC$expr[[1]]) )
	} else if(is.na(uncertMC$expr)) {
		cat("NA")
	}

	cat("\n\n")
	cat(paste("Evaluation method: ",uncertMC$method, "\n\n"))

	cat("Budget:\n")
	dp<-uncertMC$budget[sapply(uncertMC$budget, function(x) !all(is.na(x)))]
	if(!is.null(uncertMC[["distrib", exact=TRUE]]) ) {
		distrib.labels<- as.vector(
				sapply(uncertMC$distrib, function(x) if(is.function(x)) deparse(x)[1] else paste(x)) 
			)
		dp$distrib<-sub(paste("(.{",maxwidth,",",maxwidth,"})(.+)", sep=""), "\\1...",distrib.labels)
	}
	if(!is.null(uncertMC$distrib.pars)) {
		dp$distrib.pars <- vector("character", length=nrow(uncertMC$budget) )
		for(nn in row.names(dp) ) {
			dp[nn,"distrib.pars"]<-
				paste(names(uncertMC$distrib.pars[[nn]]), format(uncertMC$distrib.pars[[nn]], digits=digits), sep="=", collapse=", ")
		}
		
	}
	print.data.frame(dp,digits=digits, right=right, ...)
	if(!is.null(uncertMC$additional) ) print(as.data.frame(uncertMC$additional), ...)

	cat("\n   y: ", format(uncertMC$y))
	cat("\nu(y): ", format(uncertMC$u.y), "\n")

	if(!simplify) {
		cat("\nCovariance matrix:\n")
		print(uncertMC$cov)
		cat("\nCorrelation matrix:\n")
		print(uncertMC$cor)
		cat("\n")
		if(!is.null(uncertMC$cov.xy)) {
			cat("\nX-Y Covariances:\n")
			print(uncertMC$cov.xy)
		}
		if(!is.null(uncertMC$cor.xy)) {
			cat("\nX-Y Correlations:\n")
			print(uncertMC$cor.xy)
		}
	}
	
	cat(sprintf("\nMonte Carlo evaluation using %d replicates:\n", uncertMC$B))
	cat("\n   y:\n")
	if(simplify) {
		print(summary(uncertMC$MC$y))
	} else {
		print(uncertMC$MC$y)
		if(!is.null(uncertMC$MC$x) ) {
			cat("\n   x:\n")
			print(summary(uncertMC$MC$x))
		}
	}
	invisible(uncertMC)
}

summary.uncertMC<-function(uncertMC, digits=NULL, right=FALSE, ..., simplify=TRUE, minimise=FALSE){
	print.uncertMC(uncertMC, digits=digits, right=right, ..., simplify=simplify, minimise=minimise)
}

plot.uncertMC<-function(uncertMC, which=1:2, main=paste("Monte Carlo evaluation -",deparse(substitute(uncertMC))),
		ask = prod(par("mfcol")) < length(which) && dev.interactive(), 
		caption=list("Histogram", "Q-Q plot", "Density", "Correlation x-y", "Covariance x-y"), 
		xlab=paste(deparse(substitute(uncertMC)), "$y", sep=""),
		..., cex.caption=1.0, cex.main=1.25, lwd.y=2, col.y=2, lty.y=1,
		col.qqline=NULL, lty.qqline=NULL, lwd.qqline=NULL ) {
	
	arglist <- list(...)
	gpars<-arglist[names(arglist) %in% names( par() )]
	one.fig <- prod(par("mfcol")) == 1
	if(prod(par("mfcol"))>2) cex.caption <- cex.caption*0.8
	
	if (ask) {
		oask <- devAskNewPage(TRUE)
		on.exit(devAskNewPage(oask))
	}
	
	if( class(uncertMC)[1]=="uncertMC" ) {
		at<-NULL
		x.names<-row.names(uncertMC$budget)
		if(1 %in% which) {
			histpars<-arglist[names(arglist) %in% names(c(formals(hist.default),par()))]
			do.call(hist.default, c(list(x=uncertMC$MC$y, main="", xlab=xlab), histpars))
			if(lwd.y >= 1) abline(v=uncertMC$y, col=col.y, lwd=lwd.y, lty=lty.y)
			mtext(caption[[1]], side = 3, line=0.25, cex=cex.caption)
			if(one.fig) title(main=main)
		}

		if(2 %in% which) {
			qqpars<-arglist[names(arglist) %in% names(c(formals(qqnorm),par()))]
			if(is.null(qqpars$datax)) qqpars$datax = TRUE
			do.call(qqnorm.default, c(list(y=uncertMC$MC$y,  main=""), qqpars))
			
			qqlpars<-arglist[names(arglist) %in% names(c(formals(qqline),par()))]
			if(is.null(qqlpars$datax)) qqlpars$datax = TRUE
			if(!is.null(col.qqline)) qqlpars$col<-col.qqline
			if(!is.null(lty.qqline)) qqlpars$lty<-lty.qqline
			qqlpars$y<-uncertMC$MC$y
			do.call(qqline, qqlpars)
			mtext(caption[[2]], side = 3, line=0.25, cex=cex.caption)
			if(one.fig) title(main=main)
		}
		
		if(3 %in% which) {
			dpars<-arglist[names(arglist) %in% names(formals(density.default))]
			dpars$x<-uncertMC$MC$y
			d<-do.call(density.default, dpars)
			
			dppars<-arglist[names(arglist) %in% names(c(formals(plot.default), par()))]
			do.call(plot.density, c(list(x=d, main=""), dppars))
			if(lwd.y >= 1) abline(v=uncertMC$y, col=col.y, lwd=lwd.y, lty=lty.y)
			mtext(caption[[3]], side = 3, line=0.25, cex=cex.caption)
			if(one.fig) title(main=main)
		}

		if(4 %in% which) {
			corpars<-arglist[names(arglist) %in% names(formals(cor))]
			cor.xy<-function(x, y, p) {
			            do.call(cor, c(list(x=x, y=y), p))
			}
			
			if(!is.null(corpars$method)) {
				c.methods<- eval(formals(cor)$method)
				c.m<-c.methods[pmatch(corpars$method,c.methods)][1]
				if(is.na(c.m)) {
					warning(sprintf("Correlation method %s not found: using pearson",corpars$method) )
					corpars$method<-"pearson"
				} else corpars$method <-c.m
			} else {
				corpars$method<-"pearson"
			}
			
			if(corpars$method %in% row.names(uncertMC$cor.xy)) {
				xycor<- unlist(uncertMC$cor.xy[corpars$method, ])
			}  else {
				if(!is.null(uncertMC$MC$x)) 
				    xycor<-sapply(uncertMC$MC$x, cor.xy, y=uncertMC$MC$y, p=corpars)
				else xycor<-NULL
			} 
			
			if( !is.null(xycor) ) {
				names.arg<-if(is.null(names(xycor))) x.names else names(xycor)

				barpars<-arglist[names(arglist) %in% names(c(formals(barplot.default), par()))]
				at<-do.call(barplot, c(list(height=as.vector(xycor), names.arg=names.arg), barpars))
				corMethod<-paste(toupper(substring(corpars$method, 1,1)), 
					substring(corpars$method, 2), sep="", collapse=" ")

				c4<-if( caption[[4]] == eval(formals(plot.uncertMC)$caption)[[4]] )
					paste(corMethod,caption[[4]]) 
				      else caption[[4]]
				mtext(c4, side = 3, line=0.25, cex=cex.caption)
				if(one.fig) title(main=main)
			} else {
				warning(sprintf("Missing %s$MC$x and $cor.xy; correlation plot not made", deparse(substitute(uncertMC))), call.=TRUE)
			}
		}

		if(5 %in% which) {
			covpars<-arglist[names(arglist) %in% names(formals(cov))]
			cov.xy<-function(x, y, p) {
			            do.call(cov, c(list(x=x, y=y), p))
			}
			
			if(!is.null(covpars$method)) {
				c.methods<- eval(formals(cov)$method)
				c.m<-c.methods[pmatch(covpars$method,c.methods)][1]
				if(is.na(c.m)) {
					warning(sprintf("Covariance method %s not found: using pearson",covpars$method) )
					covpars$method<-"pearson"
				} else covpars$method <-c.m
			} else {
				covpars$method<-"pearson"
			}
			
			if(covpars$method %in% row.names(uncertMC$cov.xy)) {
				xycov<- unlist(uncertMC$cov.xy[covpars$method, ])
			}  else {
				if(!is.null(uncertMC$MC$x)) 
				    xycov<-sapply(uncertMC$MC$x, cov.xy, y=uncertMC$MC$y, p=covpars)
				else xycov<-NULL
			} 
			
			if( !is.null(xycov) ) {
				names.arg<-if(is.null(names(xycov))) x.names else names(xycov)

				barpars<-arglist[names(arglist) %in% names(c(formals(barplot.default), par()))]
				at<-do.call(barplot, c(list(height=as.vector(xycov), names.arg=names.arg), barpars))
				covMethod<-paste(toupper(substring(covpars$method, 1,1)), 
					substring(covpars$method, 2), sep="", collapse=" ")

				c5<-if( caption[[5]] == eval(formals(plot.uncertMC)$caption)[[5]] )
					paste(covMethod,caption[[5]]) 
				      else caption[[5]]
				mtext(c5, side = 3, line=0.25, cex=cex.caption)
				if(one.fig) title(main=main)
			} else {
				warning(sprintf("Missing %s$MC$x and $cov.xy; covariance plot not made", deparse(substitute(uncertMC))), call.=TRUE)
			}
		}

		if (!one.fig ) 	{	  
			mtext(main, outer = TRUE, cex = cex.main, 
				line=if(par("oma")[3L] >= 1) 0 else -1.5)
		}
		
		invisible(NULL)
		
	} else {	
		stop(paste(deparse(substitute(uncertMC))), "is not an 'uncertMC' object")
	}	
	
}

