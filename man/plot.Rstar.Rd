\name{plot.RStar}
\title{Plot method for a 'RStar' object}
\usage{
 \method{plot}{RStar}(x, what = c("RStar", "corrTerm", "logQonR", "QStatistic"), single = FALSE,
	compare = FALSE, legend = TRUE, \dots)
}
\alias{plot.RStar}
\arguments{
 \item{x}{ an object of class \code{RStar}.}
 \item{what}{ optional character string indicating what should be plotted. \code{"RStar"} plots
   the modified ratio, \code{"corrTerm"} plots the correction term, \code{"logQonR"} plots
	the logarithm of the Q statistic over the signed root of the (log) likelihood ratio, and
	\code{"QStatistic"} plots the Q statistic. Partial matching is used so only the inital character is needed.
	Defaults to "RStar".}
 \item{single}{ optional logical value indicating whether a single graph per window should be plotted.
	Defaults to \code{FALSE}.}
 \item{compare}{ optional logical value indicating whether only the modified ratio (if selected in \code{what})
   should be plotted, or also the signed root should be in the plot. Only the modified ratio is plotted
	by default.}
 \item{legend}{ optional logical value indicating whether a legend should be included in the plot, identifiying the lines.
  Defaults to \code{TRUE}.}
\item{\dots}{ some methods for the generic require additional
     arguments. None are used in this method.}
}
\description{
  Plots the modified ratio test and/or other quantities of interest. 
}
\value{
  Plots of the modified ratio test in the context of a linear mixed effects model fit.
  For every parameter present in
  the \code{RStar} object, a plot is produced, with lines indicating the 95% and 99%
  critical levels based on the standard normal approximation.
  Alternatively, other quantities are plotted
  depending on selection.
}
\author{Sigfrido Iglesias-Gonzalez \email{sigfrido@alumni.utoronto.ca}}
\seealso{
 \code{\link{intervals.RStar}},
 \code{\link{QStatistic}},
 \code{\link{plot.QStatistic}},
 \code{\link{plot.RStar}},
 \code{\link{plot.signedRoot}},
 \code{\link{RStar}},
 \code{\link{signedRoot}}
}
\examples{
# 'lmeObject' fit
Orth.fit <- lmeObject(distance ~ 1 + Sex + I(age - 11), random = pdSymm(~ 1 + I(age - 11)),
	data = Orthodont) # 'lmeObject'
hyp.scales <- list(fixed = list(c(3,3), c(3,3) , c(3,3)), reStruct = list(list(c(2.1,2.1), c(1.4,2), c(2,1.5)))) # scales
cons.Orth.fit<- conslmeObjects(Orth.fit, precision = 10, scale = hyp.scales) # 'conslmeObjects'
rs.Orth.fit <- RStar(Orth.fit, cons.Orth.fit) # computing the likelihood ratio
plot(rs.Orth.fit) # plots
plot(rs.Orth.fit, compare = TRUE, single = TRUE) # plots (comparing with the standard ratio)
}
