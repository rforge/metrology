\name{duewer.plot}
\alias{duewer.plot}
\alias{duewer.plot.default}

\title{
Duewer concordance/apparent precision plot 
}
\description{
Produces a Duewer concordance/apparent precision plot, showing relative precision or
uncertainty plotted against (relative) deviation from assigned value.
}
\usage{

duewer.plot(x, \dots)

\method{duewer.plot}{default}(x, s, mu = median(x), sigma = mad(x), s0 = median(s), 
			labels = NA, radius = 1:3, basis = c("radius", "prob"), 
			probs = 1 - (2 * (1 - pnorm(radius))), df = NA, 
			units = c("z", "x"), scale.x = TRUE, main, xlab, ylab, 
			xlim, ylim, at.xax = NULL, at.yax = NULL, aspect, 
			col.boundary = "lightgrey", lty.boundary = par("lty"), 
			lwd.boundary = par("lwd"), label.contours = T, 
			cex = par("cex"), cex.label = 0.7, pos = 3, adj = NULL, 
			pos.clab = "bottomright", col.clab = col.boundary, 
			cex.axis = par("cex.axis"), pch = par("pch"), las = par("las"), 
			col = par("col"), bg = par("bg"), ...)
}

\arguments{
  \item{x}{Numeric vector of values to be plotted.
}
  \item{s}{Numeric vector of standard deviations, standard errors or uncertainties
  of length \code{length(x)} associated with x.
}
  \item{mu}{A single location against which to compare x.
 }
  \item{sigma}{A measure of dispersion against which deviations x-mu can be compared.
}
  \item{s0}{A typical, expected or reference value for the standard uncertainties s
}
  \item{labels}{An optional vector of point labels of length \code{length(x)}. 
  \code{labels} are coerced to character on use, so may be character, factor etc..
}
  \item{radius}{A vector of radii for reference lines in the classic Duewer plot.
}
  \item{basis}{Controls the nature of the reference or contour lines drawn on
  the plot. See Details for specification.
}
  \item{probs}{A vector of probabilities for countours when \code{basis=="prob"}. 
}
  \item{df}{A single value for degrees of freedom, appropriate to the standard deviations 
  \code{s}. Ignored if \code{basis=="radius"}
}
  \item{units}{Controls scaling of the plot. If set to \code{"z"}, a classic Duewer plot 
  of \code{s/s0} vs. \code{(x-mu)/sigma} is produced. If \code{units=="x"}, the plot is 
  drawn without scaling by \code{sigma} coordinates.
}
  \item{scale.x}{Controls scaling for contours when \code{basis="prob"}. If \code{TRUE},
  	contours are scaled up in x by a factor \code{sqrt(df+1)}. This is appropriate when \code{sigma}
  	is taken from the dispersion of \code{x} when \code{x} represents mean values. If 
  	\code{sigma} represents the sd of the population from which samples of size n are averaged to 
  	give \code{x}, \code{scale.x} should be set to \code{FALSE}.
 }
  \item{main}{Main title for the plot, passed to \code{title()}.
}
  \item{xlab, ylab}{x- and y-axis labels, passed to \code{title()}.
}
  \item{xlim, ylim}{x- and y-limits for the plot. 
}
  \item{at.xax, at.yax}{Locations for x- and yaxis tick marks, passed to \code{axis()}
}
  \item{aspect}{The aspect ratio for the plot, passed to \code{plot.window}. This defaults to 1.0 for
  	\code{basis=="radius"}, giving semicircular countours, and NA otherwise.
}
  \item{col.boundary, lty.boundary, lwd.boundary}{Colour, line type and line width for contour lines.
}
  \item{label.contours}{Logical, controlling whether countour lines are labelled with approximate probabilities.
}
  \item{cex}{Expansion factor for plotted symbols.
}
  \item{cex.label}{Expansion factor for point labels.
}
  \item{pos, adj}{Specifies position/adjustment of point lables. Passed to \code{text}; see \code{?text} for details.
}
  \item{pos.clab}{Specification for location of contour labels. Options are '"top"', '"topright"', '"right"',
		'"bottomright"', '"bottom"', '"bottomleft"', '"left"', '"topleft"'. A vector can be provided
		to give multiple labels. Contour labels for \code{basis="prob"} are placed approximately 
		at the location indicated and adjusted outward appropriately. For basis="radius", 
     		"bottomright" and "bottomleft" are as for "right" and "left" but just below the x-axis,
     		and "bottom" is replaced with c("bottomright", "bottomleft").
}
  \item{col.clab}{Colour for contour labels.
}
  \item{cex.axis}{Expansion factor for axis labels.
}
  \item{las}{Axis label orientation, passed to \code{axis}.
}
  \item{pch, col, bg}{Graphical parameters passed to \code{points}.
}
  \item{\dots}{Other parameters passed to plotting functions. Currently unused.
}
}
\details{
A Duewer plot is a plot of dispersion against location. Classically, this has been applied to
multiple observations from laboratories. Locations \code{x} are mean results of the form 
\code{(x-mu)/s} and and dispersions \code{s} are the associated sd. The principle has also 
been applied to multiple results for different measurands per laboratory, by calculating 
z-scores for all observations relative to the assigned value and dispersion for each measurand
and then plotting mean and sd of the scores. More recently the plot has been used to summarise
reported values and (usually) standard uncertainties in metrology comparisons to allow
quick assessment of anomalies within data sets.

The traditional plot included visual guides in teh form of semicircular 
contours at multiples of (x-mu)/sigma for the x-axis and s/s0 for the y-axis, 
s0 being a median or other estimate of the typical sd. This corresponds to
\code{basis="radius"}. 

In this implementation, a probabilistic basis for contours is also available. This uses 
Helmert's distribution to provide constant likelihood contours for a known mean 
and standard deviation.

Because Helmert's distribution is for means and sd's of samples of a particular size
from a single population, the contours show at best what might be expected if all labs
performed similarly and when similar degrees of freedom for each sd are involved. This is
reasonably meaningful when taking averages and sd's of the same number of observations, 
but the contour probabilities are at best indicative in other circumstances.
The major advantage of the Helmert-based contours is that they reflect the low 
likelihood of very small dispersions rather better than semicircular guides around the origin.

In using the probabilistic basis, note that \code{scale.x} should be set to TRUE if \code{sigma} is a simple 
dispersion estimate taken from the observations \code{x}. For comparing with means from a population 
of known sigma, for example when summarising z-scores with expected sd of 1.0, use  \code{scale.x=FALSE}.

}
\value{
This function is called for its side effect, which is the production of a plot.
}
\references{
Duewer, D, Probably in Anal. Chem. in about 1990
%%  [### Find ref - SLRE ###] 
}
\author{S Ellison \email{s.ellison@lgc.co.uk}
}

\seealso{
\code{\link{axis}} for axis control, \code{\link{points}},  \code{\link{text}} for 
plotting parameters
}
\examples{
	require(metRology)
	data(Pb)
	Pb
	
	duewer.plot(Pb$value, Pb$u)

	duewer.plot(Pb$value, Pb$u, basis="prob", df=5)

	#Illustrate contour labelling
	duewer.plot(Pb$value, Pb$u, pos.clab="bottom")
}

\keyword{ hplot }
