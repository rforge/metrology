duewer.plot<-function(x, ...) UseMethod("duewer.plot")

duewer.plot.default<-function(x,s,mu=median(x),sigma=mad(x), s0=median(s), labels=NA,
        radius=1:3, basis=c("radius","prob"), probs=1-(2*(1-pnorm(radius))), df=NA, 
        units=c("z","x"), scale.x=TRUE,
        main, xlab, ylab, 
        xlim, ylim,  at.xax=NULL, at.yax=NULL, aspect, 
        col.boundary="lightgrey", lty.boundary=par("lty"), lwd.boundary=par("lwd"),
        label.contours=T, 
        cex=par("cex"), cex.label=0.7, pos=3, adj=NULL, 
        pos.clab="bottomright", col.clab=col.boundary,
        cex.axis=par("cex.axis"), pch=par("pch"), las=par("las"), 
        col=par("col"), bg=par("bg"), ...) {
        #Produces a Duewer concordance/apparent precision plot 
        #of f=s/s0 against z=(x-mu)/sigma, centred on zero, with
        #(semicircular) boundary lines at (f^2+s^2)=radius
        #If basis=="prob" the boundaries are drawn at constant likelihood ratio
        #
        
        #units: if units=="z", plot z-score and s/s0; if "x" plot x and s with boundaries scaled to 
        #       radius*sigma
        
        #pos.clab:      Contour labels for prob basis are placed approximately at the location indicated
        #               and adjusted outward appropriately. Options are ‘"top"’, ‘"topright"’, ‘"right"’,
        #               ‘"bottomright"’, ‘"bottom"’, ‘"bottomleft"’, ‘"left"’, ‘"topleft"’.
        #               For basis="prob" the positions are taken as provided. For basis="radius", 
        #               "bottomright" and "bottomleft" are as for "right" and "left" but just below the x-axis,
        #               and "bottom" is replaced with c("bottomright", "bottomleft").
        
        
        pos.clab <- match.arg(pos.clab, choices=c("top", "topright", "right",
                        "bottomright", "bottom", "bottomleft", "left", "topleft"), several.ok=TRUE)
                        
        units <- match.arg(units)
        
        if(units=="x") {
                z<-(x-mu)
                f<-s
                x.scale<-sigma
        } else {
                z<-(x-mu)/sigma
                f<-s/s0
                x.scale<-1
        }       

        basis<-match.arg(basis)
        
        ###The following may need to be basis-specific...
        if(missing(xlim)) {
                xlim=c(-1,1)*max(c(radius*x.scale,pretty(abs(z))), na.rm=T)
        }

        if(missing(ylim)) {
                ylim=c(0, max(c(radius*x.scale,pretty(abs(f))), na.rm=T) )
        
        }
        
        
        if(missing(aspect)){
                if(basis=="radius")
                        aspect<-1.0
                else
                        aspect<-NA
        }
        
        plot.new()
        
        plot.window(xlim=xlim, ylim=ylim, asp=aspect)
        
        ###
        
        #Calculate boundaries
        t<-seq(0,pi, length.out=200)
        max.sx<-0
        if(is.na(basis) || basis=="radius") {
                xpd <- par(xpd=TRUE)
                for(r in radius) {
                        lines(x.x <- r*x.scale*cos(t), s.x <- r*x.scale*sin(t), col=col.boundary, lty=lty.boundary, lwd=lwd.boundary)
                                max.sx <- max(max.sx, max(s.x))
                        if(label.contours) {
                                if("bottom" %in% pos.clab) {
                                        pos.clab <- unique(c(pos.clab[-match("bottom", pos.clab)], "bottomleft", "bottomright"))
                                }
                                for(i in 1:length(pos.clab)) {
                                        l.pos <- switch(pos.clab[i], top=100, topright=50, right=1,
                                                        bottomright=1, bottom=1, bottomleft=200, 
                                                        left=200, topleft=150)
                                        l.adj <- switch(pos.clab[i], top=c(-0.1,-0.2), topright=c(-0.1,-0.2), right=c(-0.1,-1),
                                                        bottomright=c(0,1.1), bottom=c(-0.1,-1), bottomleft=c(1,1.1), 
                                                        left=c(1.1,-1), topleft=c(1,-0.2))
                                        text( x.x[l.pos], s.x[l.pos], labels=paste("p=",format(probs[match(r, radius)], digits=3), sep=""), 
                                                        adj=l.adj, cex=cex.label, col=col.clab)
                                }
                        }
                }
                par(xpd=xpd)

        } else {
                if(basis=="prob") {
                        #df is _required_ for probabilistic calculation, so warn if not given

                        if(is.na(df)) stop("df is required if a probabilistic basis is chosen")

                        #will need n:
                        
                        n<-df+1
                        
                        #Calculate probabilistic boundaries
                        #NB: probs defaults to 1-2*pnorm(radius),
                        #on the assumption that radius is set for 2-sided 95% intervals on z.
                        #Boundaries b.x, b.s, for s are calculated as likelihood ratio
                        #boundaries for Helmert's distribution* based on qchisq(probs, 2)
                        #See Y Pawitan, (2001) In all likelihood: Statistical 
                        #modelling and Inference Using Likelihood, 
                        #Clarendon Press, Oxford, pp258-9 for the reason that this
                        #corresponds to a density ratio equal to probs.
                        #
                        #*Helmert's distribution: See
                        #W Kruskal, American Mathematical Monthly 53, 435-438, (1946)
                        
                        #Helmert's distribution has density
                        #K1*( s^(n-2))*exp( (-n/(2*sigma^2))*(s^2+u^2) )
                        #where u is the error x-mu and
                        #K1 = ( n^(n/2) ) / ( 2^((n-2)/2) *sqrt(pi) * Gamma((n-1)/2) * sigma^2 )
                        #
                        #K1 cancels out of the likelihood ratio, so we ignore it and solve for 
                        #s given u at required density 
                        #
                        
                        #Need max density, which is at sqrt((n-2)/n). Ignoring K1 and with sigma=1:
                        smax<-sqrt((n-2)/(n))
                        
                        dhmax<-.dhelmert(n=n, u=0, s=smax, sigma=1)
                        
                        cat(paste("Calculating for probs=c(",paste(round(probs, 3), collapse=", "),")\n"))

                        for(p in probs) {
                                #Generate smin...smax
                                smin<-.qdhelmert(n=n, d=dhmax*(1-p), u=0, lower.tail=TRUE)
                                smax<-.qdhelmert(n=n, d=dhmax*(1-p), u=0, lower.tail=FALSE)
                                #Generate sine-spaced sequence (produces a smoother contour)
                                b.s <- smin+(smax-smin)*sin(seq(0, pi/2, length.out=100))
                                #solve for x
                                b.x <- .qdhelmert(n=n, d=dhmax*(1-p), s=b.s[c(-1,-length(b.s))] )
                                if(scale.x) b.x<-sqrt(n)*b.x
                                
                                b.s<-c(b.s, rev(b.s)[c(-1,-length(b.s))])
                                b.x<-c(0,b.x, 0,-rev(b.x))
                                polygon(b.x*x.scale,b.s*x.scale, border=col.boundary, lty=lty.boundary, col=NA)
                                if(label.contours) {
                                        for(i in 1:length(pos.clab)) {
                                                l.pos <- switch(pos.clab[i], top=100, topright=70, right=which.max(b.x),
                                                                bottomright=10, bottom=1, bottomleft=190, 
                                                                left=which.min(b.x), topleft=130)
                                                l.adj <- switch(pos.clab[i], top=c(-0.1,-0.2), topright=c(-0.1,-0.2), right=c(-0.1,0.5),
                                                                bottomright=c(0,1), bottom=c(-0.1,1.1), bottomleft=c(1,1), 
                                                                left=c(1.1,0.5), topleft=c(1,-0.2))
                                                text( b.x[l.pos], b.s[l.pos], labels=paste("p=",format(p, digits=3), sep=""), 
                                                        adj=l.adj, cex=cex.label, col=col.clab)
                                        }
                                }
                        }
                }
                
        }

        usr<-par("usr")
        
        segments(c(usr[1],0),c(0,0 ), c(usr[2],0), c(0,usr[4]))
        
        axis(1, pos=0, las=las, cex.axis=cex.axis, at=at.xax)
        
        if(is.null(at.yax)) {
                at.yax <- pretty(c(0, usr[4]), n=4)
                at.yax <- at.yax[at.yax > 0]
        } 
        axis(2, pos=0, las=las, cex.axis=cex.axis, at=at.yax)
        
        points(z,f,pch=pch, col=col, bg=bg, cex=cex)
        
        if(!all(is.na(labels))) text(z, f, as.character(labels), pos=pos, adj=adj, cex=cex.label)
        
        if(!missing(main)) title(main=main)
        if(!missing(ylab)) title(ylab=ylab)
        if(!missing(xlab)) {
                h <- strheight(xlab, units = "user")
                text(0, -h*4, labels=xlab)
        }
        
        

}

.dhelmert <- function(n, u=0, s=1, sigma=1) {
        K1<-( n^(n/2) ) / ( 2^((n-2)/2) *sqrt(pi) * gamma((n-1)/2) * sigma^2 )
        K1*(s^(n-2))*exp( (-n/(2*sigma^2))*(s^2+u^2) )
}

.qdhelmert<-function(n,d,u,s,sigma=1, umax=6, smax=6*sigma, lower.tail=TRUE) {
        if(missing(u) && missing(s)) stop("One of u or s must be given")
        find.u <- function(u, n,d,s,sigma) .dhelmert(n=n,u=u,s=s,sigma=sigma )-d
        find.s <- function(s, n,d,u,sigma) .dhelmert(n=n,u=u,s=s,sigma=sigma )-d
        if(missing(u)) {
                #find u
                objective <- find.u
                u <- vector(length=length(s))
                for(i in 1:length(u)) {
                        u[i] <- uniroot(objective, c(0, umax), n=n, s=s[i], d=d, sigma=sigma)$root
                }
                return(u)
        } else if(missing(s)){
                #find s
                s.max<-sigma*sqrt((n-2)/(n))
                #find u
                objective <- find.s
                s <- vector(length=length(u))
                if(lower.tail) interval <- c(0,s.max) 
                        else interval <- c(s.max, smax)
                for(i in 1:length(s)) {
                        s[i] <- uniroot(objective, interval, n=n, u=u[i], d=d, sigma=sigma)$root
                }
                return(s)
        } else {
                stop("Only one of u and s should be given")
        }
}


