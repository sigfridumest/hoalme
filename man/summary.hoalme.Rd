
\name{summary.hoalme}
\title{Summary method for objects of class \code{hoalme}}
\usage{
\method{summary}{hoalme}(object,  level = .95, add.tests = c("all", "none", "r", "r~"), sig = 5, \dots)
}
\alias{summary.hoalme}
\arguments{
 \item{object}{an object of class \code{hoalme}.}
 \item{level}{a positive numeric value representing the confidence
 level for computing confidence intervals. Defaults to 0.95.}
 \item{add.tests}{a character string.  If \code{"all"}, the standard likelihood
 ratio and every modified ratio are computed to produce the summary.
 Choose \code{"none"} for computing only the modified ratio using the Fraser and Reid approach.
 Optionally, a particular modified ratio can be
 chosen by specifying \code{"r*"} or \code{"r~"}. Defaults to \code{"all"}.}
 \item{sig}{a positive integer value representing the number of significant decimals to show. Defaults to 5.}
\item{\dots}{Some methods for the generic require additional arguments. None are used in this method.}
}
\description{
  Tests for the fixed effects parameters in a mixed linear model fit are computed, together with confidence intervals.
  For the variance parameters only confidence intervals
  are computed. The returned object is printed by the \code{print.summary.hoalme} method.
}
\value{
  An object of class \code{summary.hoalme}. Components of this object
  are \code{tests.RStar}, \code{tests.signedRoot}, and 
  \code{intervals.RStar} and \code{intervals.signedRoot}.
  Later versions will include some other results based on other procedures.
}

\author{Sigfrido Iglesias-Gonzalez \email{sigfrido@alumni.utoronto.ca}}
\seealso{
 \code{\link{hoalme}},
 \code{\link{intervals.RStar}},
 \code{\link{intervals.signedRoot}},
 \code{\link{RStar}},
 \code{\link{signedRoot}}
}
\examples{
#lmeObject
Orth.fit <- lmeObject(distance ~ 1 + Sex + I(age - 11), random = pdSymm(~ 1 + I(age - 11)), data = Orthodont)
# specyfing hypotheses
hyp.scales <- list(fixed = 3, reStruct = list( list(c(2.1, 2.1), c(1.4, 2), c(2, 1.5)) ) )
#constrained optimization
cons.Orth.fit <- conslmeObjects(Orth.fit, scale = hyp.scales)
# 'hoalme' object
hl.Orth.fit <- hoalme(Orth.fit, cons.Orth.fit)
# test and confidence intervals
summary(hl.Orth.fit) 
}
