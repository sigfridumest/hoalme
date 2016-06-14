\name{V}
\title{V vectors indicating ancillary directions in a linear mixed model fit.}
\usage{
 \method{V}{lmeObject}(object, \dots)
}
\alias{V}
\alias{V.lmeObject}
\arguments{
 \item{object}{ an object of class \code{lmeObject}.}
 \item{\dots}{ some methods may require additional arguments. None are used in this method.}
}
\description{
  This method extracts the V vectors indicating the direction in which the sample space derivatives
  should be computed in a linear mixed effects model. These derivatives are the basic ingredient
  in the computation of the sample space derivatives statistic "Q" by Fraser and Reid.
}
\value{
  An object of class \code{V}, which is a list, containing the V vectors for each subject present in the
  \code{lmeObject} object. Each component in the list contains the V vectors of each subject.
}
\author{Sigfrido Iglesias-Gonzalez \email{sigfrido@alumni.utoronto.ca}}
\seealso{
  \code{\link{lmeObject}},
  \code{\link{RStar}}
}
\examples{
Orth.fit <- lmeObject(distance ~ 1 + Sex + I(age - 11), random = pdSymm(~ 1 + I(age - 11)),
		data = Orthodont) # 'lmeObject'
V(Orth.fit) #V vectors
}
