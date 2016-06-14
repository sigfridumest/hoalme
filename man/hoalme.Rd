
\name{hoalme}
\title{'hoalme' methods for conducting Higher Order Inference in the Linear Mixed-Effects Models}
\usage{
hoalme(object, \dots)
 \method{hoalme}{lmeObject}(object, conslmeObjs, \dots)
 \method{hoalme}{RStar}(object, \dots)
}
\alias{hoalme}
\alias{hoalme.lmeObject}
\alias{hoalme.RStar}
\arguments{
 \item{object}{an object for which higher order asymptotic methods are desired to be computed.}
 \item{conslmeObjs}{an object of class \code{conslmeObjects}, such as the one resulting from using the
  \code{conslmeObjects} constructor.}
 \item{\dots}{some methods for the generic may require additional
     arguments. None are used in these methods.} 
}
\description{
  This generic function generates an object of class \code{hoalme}, which provides Higher Order Inference
  on the parameters of a fitted model. The methods here included are only for
  linear mixed effects models. Use \code{summary} method for hypothesis tests
  and confidence intervals.
}
\value{
  An object of class \code{hoalme}, containing components necessary to compute hypothesis tests and confidence
  intervals based on Higher Order Asymptotic theory. Techniques implemented in this version include only the modified ratio
  using results by  Fraser and Reid (2001).
}

\author{Sigfrido Iglesias-Gonzalez \email{sigfrido@alumni.utoronto.ca}}
\seealso{
  \code{\link{lmeObject}},
  \code{\link{conslmeObjects}},
  \code{\link{RStar}},
  \code{\link{summary.hoalme}}
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
}
