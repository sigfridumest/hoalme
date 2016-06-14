\name{update.RStar}
\title{update method for an 'RStar' object}
\usage{
 \method{update}{RStar}(object, delta, ...)
}
\alias{update}
\alias{update.RStar}
\arguments{
 \item{object}{ an object of class \code{RStar}.}
 \item{delta}{ a positive number specifying the absolute values of the modified ratio less than \code{delta}
  that to be disregarded in the smoothing.}
 \item{\dots}{ some methods for the generic require additional arguments. None are used in this method.}
}
\description{
  Updates the smoothing of an \code{RStar} object.
}
\value{
  An \code{RStar} object.
}
\author{Sigfrido Iglesias-Gonzalez \email{sigfrido@alumni.utoronto.ca}}
\examples{
# 'lmeObject' fit
Orth.fit <- lmeObject(distance ~ 1 + Sex + I(age - 11), random = pdSymm(~ 1 + I(age - 11)),
	data = Orthodont) # 'lmeObject'
hyp.scales <- list(fixed = list(c(3,3), c(3,3) , c(3,3)), reStruct = list(list(c(2.1,2.1), c(1.4,2), c(2,1.5)))) # scales
cons.Orth.fit<- conslmeObjects(Orth.fit, precision = 10, scale = hyp.scales) # 'conslmeObjects'
rs.Orth.fit <- RStar(Orth.fit, cons.Orth.fit) # computing the likelihood ratio
plot(rs.Orth.fit) # plots
#updating r* computation setting a smaller neighborhood of the hypothesis
rs.Orth.d05 <- update(rs.Orth.fit, .5)
#exhibits some odd points
plot(rs.Orth.d05, compare = TRUE)
}
