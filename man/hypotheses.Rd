
\name{hypotheses}
\title{Extract or compute a set of hypotheses or theoretical parameter values}
\usage{
 hypotheses(object, ...)
}
\alias{hypotheses}
\arguments{
 \item{object}{an object which is expected to provide information about hypothetical
  values of a model.}
 \item{\dots}{some methods for this generic may require additional arguments.}
}
\description{
 Hypothetical parameter values present in \code{object}. These are usually hypotheses to be tested
 or parameter values at which a profile likelihood is computed.
}
\value{
  A set of hypothetical values (usually a matrix) of the parameters present in a model.
}

\author{Sigfrido Iglesias-Gonzalez \email{sigfrido@alumni.utoronto.ca}}
\seealso{
 \code{\link{hypotheses.conslmeObjects}}
}
