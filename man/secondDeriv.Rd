\name{secondDeriv}
\title{Second Derivatives of the log likelihood}
\usage{
 secondDeriv(object, \dots)
 \method{secondDeriv}{lmeObject}(object, \dots)
}
\alias{secondDeriv}
\alias{secondDeriv.lmeObject}
\arguments{
 \item{object}{ an object from which to compute or extract second derivatives of the (log)
  likelihood evaluated at the maximum likelihood estimator.}
 \item{\dots}{ some methods require additional arguments. None are used in the method
  here described.}
}
\description{
  This method extracts the matrix of second derivatives of the log likelihoof evaluated at
  the maximum likelihood estimator. The method here provided gives the
  matrix of second derivatives of the (log) likelihood of a linear mixed model fit.
}
\value{
  An object of class \code{secondDeriv}, which is a matrix.
}
\author{Sigfrido Iglesias-Gonzalez \email{sigfrido@alumni.utoronto.ca}}
\seealso{
  \code{\link{lmeObject}},
}
\examples{
# second derivatives
Orth.fit <- lmeObject(distance ~ 1 + Sex + I(age - 11), random = pdSymm(~ 1 + I(age - 11)),
									data = Orthodont) # 'lmeObject'
secondDeriv(Orth.fit) # second derivatives matrix
}
