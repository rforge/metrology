\name{metRology-package}
\alias{metRology-package}
\alias{metRology}
\docType{package}
\title{
\packageTitle{metRology}
}
\description{
\packageDescription{metRology}
}
\details{
The metRology package includes functions for:
  \itemize{ 
	\item Plotting for Key Comparisons (dot-and-bar, consistency)
	\item Uncertainty evaluation using algebraic or numeric differentiation, with support for correlation
	\item Monte Carlo evaluation of uncertainty (including correlation for normally distributed variables)
	\item Classes and functions for location estimates for metrology comparisons
	\item Mandel's h and k statistics and plots for interlaboratory studies
	\item Support functions for an excel interface
  } 
This version (0.9-23) 
  \itemize{
	\item includes a REML location estimate in addition to vr.mle
	\item adds a Youden plot (see \code{\link{yplot}}), with various options for confidence ellipses
	\item adds an improved version of \code{\link{msd}} which is faster and less memory-intensive 
	\item enables log and log.p in dt.scaled
	\item amends plot.drop1.uncert to give a plot for each measure of change specified in \code{which}
	\item corrects a grep warning appearing in drop1.uncert
	\item corrects an unnecessary 'missing \code{u}' error message in \code{uncert()} when \code{cov} 
	      was specified and \code{u} was not*
	}
*This is the sole change from version 09.22.
}
\author{
\packageAuthor{metRology}

Maintainer: \packageMaintainer{metRology}
}
%\references{
%~~ Literature or other references for background information ~~
%}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%~~ the R documentation directory ~~
\keyword{ package }
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
%\examples{
%~~ simple examples of the most important functions ~~
%}
