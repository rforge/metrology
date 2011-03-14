#Provide a simple constructor for uncert class objects, largely to 
#centralise code maintenance and avoid unnecessary inter-function calls.

.construct.uncert<-function(call.uncert=NULL, expr=NULL, y=NULL,u.y=NULL, method=NULL, 
		x=NULL, df=NULL, u=NULL, ci=NULL, u.c=NULL, x.names=NULL, 
		cov=NULL, cor=NULL, deriv=NULL, distrib=NULL, distrib.pars=NULL, 
		cov.xy=NULL, cor.xy=NULL, delta=NULL, MC=NULL, B=NULL, keep.x=NULL, ...)
	{
	rv<-list()
	rv$call <- if(!is.null(call.uncert)) call.uncert else NA
	rv$expr <- if(!is.null(expr)) expr else NA
	rv$y <- if(!is.null(y)) y else NA
	rv$u.y <- if(!is.null(u.y)) u.y else NA
	rv$method <- if(!is.null(method)) method else NA
	
	L<-max(length(x), length(u), length(ci), length(df), length(distrib))
	if(L>0) {	#Budget omitted if none of the above are included
		rv$budget<-data.frame(
			x=if(!is.null(x)) unlist(x) else rep(NA, L),
			u=if(!is.null(u)) unlist(u) else rep(NA, L),
			df=if(!is.null(df)) unlist(df) else rep(NA, L),
			c=if(!is.null(ci)) unlist(ci) else rep(NA, L)
			)
		rv$budget$u.c<-if(!is.null(u.c)) unlist(u.c) else rv$budget$u*rv$budget$c
		if(max(abs(rv$budget$u*rv$budget$c-rv$budget$u.c), na.rm = TRUE) > 
			2*.Machine$double.eps ) {
			warning("Supplied value for u*c is inconsistent with supplied u and c") 
		}
		if( is.null( x.names ) ) {
			if(!is.null(names(x)) ) {
				x.names <- names(x)
			} else if(!is.null(names(u)) ) {
				x.names <- names(u)
			} else if( !is.null(names(c)) ) {
				x.names <- names(c)
			} else x.names <- paste("X",1:L, sep="")
		}
		row.names(rv$budget) <- x.names
		
	}
	if(!is.null(distrib) ) { #Omitted if not provided
		#Add distrib (can be vector)
		if(is.list(distrib)) {
			if(is.null(names(distrib))) 
				names(distrib)<-x.names
			
			if(!all(names(distrib) %in% x.names)) {
				warning("names(distrib) do not match x.names: omitting distrib")
				rv$distrib<-NULL
			} else {
				#fill in any missing names
				missing.names<-x.names[!(x.names %in% names(distrib))]
				if(length(missing.names)>0) distrib[[missing.names]] <- NA

				#order as x
				rv$distrib<-distrib[x.names]

			}
		} else {
			rv$distrib<-as.list(distrib)
			if( is.null(names(rv$distrib)) ) names(rv$distrib) <- x.names
		}

	}
	#Add distrib.pars
	if(!is.null(rv$distrib) && !is.null(distrib.pars)) { #Omitted if not provided
		if(is.null(names(distrib.pars))) 
			names(distrib.pars<-x.names)
			
		if(!all(names(distrib.pars) %in% x.names)) {
			warning("names(distrib.pars) do not match x.names: omitting distrib.pars")
			rv$distrib.pars<-NULL
		} else {
			#fill in any missing names
			missing.names<-x.names[!(x.names %in% names(distrib.pars))]
			if(length(missing.names)>0) distrib.pars[[missing.names]] <- NA
			
			#order as x
			rv$distrib.pars<-distrib.pars[x.names]
			
		}
	}
	
	if(!missing(...)) rv$additional <- list(...)
	
	if(!is.null(cov)) rv$cov <- cov
	if(!is.null(cov) && is.null(dimnames(rv$cov))) dimnames(rv$cov)<-list(x.names, x.names)
	
	if(!is.null(cor)) rv$cor <- cor
	if(!is.null(cor) && is.null(dimnames(rv$cor))) dimnames(rv$cor)<-list(x.names, x.names)

	if(!is.null(cov.xy)) {
		if(is.null(names(cov.xy)))  names(cov.xy)<-x.names
		if(!(is.data.frame(cov.xy)))  cov.xy<-as.data.frame( t(cov.xy) )
		rv$cov.xy <- cov.xy
	}

	if(!is.null(cor.xy)) {
		if(is.null(names(cor.xy)))  names(cor.xy)<-x.names
		if(!(is.data.frame(cor.xy)))  cor.xy<-as.data.frame( t(cor.xy) )
		rv$cor.xy <- cor.xy
	}

	if(method %in% c("NUM", "kragten", "k2") ) {
		rv$delta=delta
	}
	
	rv$deriv <- if(!is.null(deriv)) deriv
	rv$B <- if(!is.null(B)) B
	rv$keep.x <- if(!is.null(keep.x)) keep.x
	rv$MC <- if(!is.null(MC)) MC
	
	class(rv) <- c("uncert", "list")
	
	return(rv)

}

print.uncert<-function(uncert, digits=NULL, right=FALSE, ..., simplify=TRUE){
	maxwidth<-12L
	cat("\nUncertainty evaluation\n\n")
	cat("Call:\n  ",deparse(uncert$call), sep="")
	cat("\n\n")
	cat("Expression: ")
	if(class(uncert$expr)=="formula" ) {
		cat(paste(uncert$expr, collapse=""))
	} else if(is.function(uncert$expr)) {
		cat( deparse(uncert$expr)[1] )
	} else if(is.expression(uncert$expr)) {
		cat( deparse(uncert$expr[[1]]) )
	} else if(is.na(uncert$expr)) {
		cat("NA")
	}
	cat("\n\n")
	cat(paste("Evaluation method: ",uncert$method, "\n\n"))
	cat("Uncertainty budget:\n")
	dp<-uncert$budget[sapply(uncert$budget, function(x) !all(is.na(x)))]
	if(!is.null(uncert[["distrib", exact=TRUE]]) ) {
		fdp<-function(x) { if(is.function(x)) deparse(x)[1] else paste(x) }
		distrib.labels<- as.vector( sapply(uncert$distrib, fdp ) )
		dp$distrib<-sub(paste("(.{",maxwidth,",",maxwidth,"})(.+)", sep=""), 
			"\\1...",distrib.labels)
	}
	if(!is.null(uncert$distrib.pars)) {
		dp$distrib.pars <- vector("character", length=nrow(uncert$budget) )
		for(nn in row.names(dp) ) {
			dp[nn,"distrib.pars"]<-
				paste(names(uncert$distrib.pars[[nn]]), 
					format(uncert$distrib.pars[[nn]], digits=digits),
					sep="=", collapse=", ")
		}
		
	}
	print.data.frame(dp,digits=digits, right=right, ...)
	if(!is.null(uncert$additional) ) {
		cat("Additional parameters:\n")
		print(as.data.frame(uncert$additional), ...)
	}
	cat("\n   y: ", format(uncert$y))
	cat("\nu(y): ", format(uncert$u.y), "\n")
	
	if(!simplify) {
		cat("\nCovariance matrix:\n")
		print(uncert$cov)
		cat("\nCorrelation matrix:\n")
		print(uncert$cor)
		cat("\n")
		if(!is.null(uncert$cov.xy)) {
			cat("\nX-Y Covariances:\n")
			print(uncert$cov.xy)
		}
		if(!is.null(uncert$cor.xy)) {
			cat("\nX-Y Correlations:\n")
			print(uncert$cor.xy)
		}
	
	}
	if(!is.null(uncert$MC) ) {
		cat(sprintf("\nMonte Carlo evaluation using %d replicates:\n", uncert$B))
		cat("\n   y:\n")
		if(simplify) {
			print(summary(uncert$MC$y))
		} else {
			print(uncert$MC$y)
			if(!is.null(uncert$MC$x) ) {
				cat("\nMC x:\n")
				print(uncert$MC$x)
			}
		}
	}
	invisible(uncert)
}


summary.uncert<-function(uncert, ... ,  simplify=TRUE){

	print.uncert(uncert, ..., simplify=simplify)	
	invisible(uncert)
}

plot.uncert<-function(uncert, which=c(1,2,4,5), main=paste(deparse(substitute(uncert))),
		ask = prod(par("mfcol")) < length(which) && dev.interactive(), 
		caption=list("Variance and covariance contributions",
		expression(sqrt(group("|","Variance and covariance contributions","|"))),
		expression("Contribution "*u[i](y)==c[i]*u[i]), "Combined contribution",
		"Correlation (x,y)", "Covariances (x,y)"),
		cex.caption=1.0, ...) {
	
	one.fig <- prod(par("mfcol")) == 1
	if(prod(par("mfcol"))>2) cex.caption <- cex.caption*0.8
	
	if (ask) {
		oask <- devAskNewPage(TRUE)
		on.exit(devAskNewPage(oask))
	}
	
	if(is.null(uncert$MC)) {
		at<-NULL
		x.names<-row.names(uncert$budget)
		contrib<-with(uncert$budget, outer(c,c,"*"))*uncert$cov
		covars<- as.vector(uncert$cov %*% uncert$budget$c)
		cors<- covars/(uncert$budget$u*uncert$u.y)
		if(1 %in% which) {
			indices<-which(abs(uncert$cor)>2*.Machine$double.eps & 
						upper.tri(uncert$cor), arr.ind=TRUE)
			x<-c(diag(contrib), 2*contrib[indices])
			names.arg<-c(x.names, paste(x.names[indices[,1]], 
					x.names[indices[,2]], sep=":"))
			at<-barplot(x, names.arg=names.arg, ...)
			mtext(caption[[1]], side = 3, line=0.25, cex=cex.caption)
			if(one.fig) title(main=main)
		}

		if(2 %in% which) {
			indices<-which(abs(uncert$cor)>2*.Machine$double.eps & 
					upper.tri(uncert$cor), arr.ind=TRUE)
			x<-sqrt(abs(c(diag(contrib), 2*contrib[indices])))
			names.arg<-c(x.names, paste(x.names[indices[,1]], 
					x.names[indices[,2]], sep=":"))
			at<-barplot(x, names.arg=names.arg, ...)
			mtext(caption[[2]], side = 3, line=0.25, cex=cex.caption)
			if(one.fig) title(main=main)
		}
		if(3 %in% which) {
			x<-uncert$budget$u.c
			names.arg<-x.names
			at<-barplot(x, names.arg=names.arg, ...)
			mtext(caption[[3]], side = 3, line=0.25, cex=cex.caption)
			if(one.fig) title(main=main)
		}

		if(4 %in% which) {
			x<-rowSums(contrib)+colSums(contrib)-diag(contrib)
			names.arg<-x.names
			at<-barplot(x, names.arg=names.arg, ...)
			mtext(caption[[4]], side = 3, line=0.25, cex=cex.caption)
			if(one.fig) title(main=main)
		}

		if(5 %in% which) {
			names.arg<-x.names
			at<-barplot(cors, names.arg=names.arg, ...)
			mtext(caption[[5]], side = 3, line=0.25, cex=cex.caption)
			if(one.fig) title(main=main)
		}
		if(6 %in% which) {
			names.arg<-x.names
			at<-barplot(covars, names.arg=names.arg, ...)
			mtext(caption[[6]], side = 3, line=0.25, cex=cex.caption)
			if(one.fig) title(main=main)
		}

		if (!one.fig ) 		  
			mtext(main, outer = TRUE, cex = par("cex.main"), 
				line=if(par("oma")[3L] >= 1) 0 else -1.5)

		invisible(at)
		
	} else {	
		plot.uncertMC(uncert$MC, main=main, ask=ask, cex.caption=cex.caption, ...)
		
	}	
	
}


update.uncert<-function(uncert, expr=NULL, method=NULL, x=NULL, u=NULL, c=NULL, 
		df=NULL, cov=NULL, cor=NULL, distrib=NULL, distrib.pars=NULL, 
		delta=NULL, B=NULL, keep.x=NULL, ...)
	{
	#Note restriction of update to _input_ parameters.
	#Note also that x.names cannot be updated
	
	#uncert<-uncert
	
	x.names<-row.names(uncert$budget)
	
	
	#get supplied params:
	nf<-names(formals())
	nf<-nf[!(nf %in% c('uncert', '...'))]
	
	bnames<-list(x="x", u="u", df="df", c="c", u.c="u.c")
	
	found<-list()
	
	for(n in nf) {
		g<-get(n)
		if(!is.null(g)) {
			found[[n]] <- g
			ulg<-unlist(g)
			if(!is.null(nb <- bnames[[n]])) {
				#A $budget parameter
				if(!is.null(names(ulg))) 
					uncert$budget[names(ulg),nb] <- ulg
				else
					uncert$budget[,nb] <- ulg
			
			} else {
				if(is.null(uncert[[n]])) {
					uncert[[n]]<-g
				} else {
					if(!is.null(ng<-names(g))){
						for(nng in ng) uncert[n][[nng]] <- g[[nng]]
					} else {
						uncert[[n]]<-g
					}
				}
			}
		
		}
	}
	#Update $additional if required
	l.<-list(...)
	for(na in names(l.)) uncert$additional[[na]] <- l.[[na]]
	
	#Extract $budget components and apply names...
	blist<-lapply(uncert$budget, function(x,nn) {names(x)<-nn; as.list(x)}, 
					nn=row.names(uncert$budget))
	
	#Extract parameters that match original function call formals
	name.fun <-deparse(uncert$call[[1]])
	argnames.fun<-names(formals(name.fun))
	allnames <- c(names(uncert$call)[-1], names(found))	#excludes '...'. uncert
	passnames <- unique(allnames[ allnames %in% argnames.fun ])
	passlist<-c(unclass(uncert), blist)
	passlist$FUN <- passlist$formula <- uncert$expr
	passlist<-c(passlist[passnames], uncert$additional)
	
	#call the original constructor with updated parameters:
	rv <- do.call(name.fun, passlist)
	
	rv$call <- match.call()
	
	return(rv)
}
