\name{hypotheses.conslmeObjects}
\title{Extract the set of hypotheses used in the constrained optimization}
\usage{
 \method{hypotheses}{conslmeObjects}(object, natural = TRUE, ...)
}
\alias{hypotheses.conslmeObjects}
\arguments{
 \item{object}{ a \code{conslmeObjects} object.}
 \item{natural}{a logical value indicating whether a set of hypotheses for
  mixed effects linear model in a natural parameterization should be
  returned, or in its unconstrained form.}
 \item{\dots}{some methods for the generic may require additional arguments. None are used in this method.}
}
\description{
  Constrained optimization of the (log) likelihood of a linear mixed model is performed on a
  set of hypothetical values of each of the parameters present in the model.
  These values are recovered through this method.
}
\value{
  A matrix with the parameter values used in the constrained optimization.}

\author{Sigfrido Iglesias-Gonzalez \email{sigfrido@alumni.utoronto.ca}}
\seealso{
 \code{\link{conslmeObjects}},
 \code{\link{hypotheses}}
}
\examples{
#lmeObject
Orth.fit <- lmeObject(distance ~ 1 + Sex + I(age - 11), random = pdSymm(~ 1 + I(age - 11)), data = Orthodont)
# specyfing hypotheses
hyp.scales <- list(fixed = 3, reStruct = list( list(c(2.1, 2.1), c(1.4, 2), c(2, 1.5)) ) )
#constrained optimization
cons.Orth.fit <- conslmeObjects(Orth.fit, scale = hyp.scales)
hypotheses(cons.Orth.fit, FALSE) # hypothetical parameter values actually used in the constrained optimization
hypotheses(cons.Orth.fit, TRUE) # hypothetical parameter values in natural parameterization
}
