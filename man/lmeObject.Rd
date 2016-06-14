\name{lmeObject}
\title{lmeObject methods for fitting Linear Mixed-Effects Models}
\usage{
lmeObject(fixed, \dots)
 \method{lmeObject}{formula}(fixed, data, random, correlation, weights, subset, method,
    na.action, control, contrasts = NULL, keep.data = TRUE, start = NULL, \dots)
 \method{lmeObject}{lme}(fixed, \dots)
}
\alias{lmeObject}
\alias{lmeObject.formula}
\alias{lmeObject.hoalme}
\alias{lmeObject.lme}
\arguments{
 \item{fixed}{either an object of class \code{lme}, representing
  a fitted linear mixed-effects model given as by the \code{lme}
  method from the \code{nlme} package,
  or a two-sided linear formula object describing the
   fixed-effects part of the model, with the response on the left of a
   \code{~} operator and the terms, separated by \code{+} operators.}
  \item{data}{an optional data frame containing the variables named in
   \code{fixed}, \code{random}, \code{correlation}, \code{weights}, and
   \code{subset}.  By default the variables are taken from the
   environment from which \code{lme} is called.}
  \item{random}{optionally, any of the following: (i) a one-sided formula
   of the form \code{~x1+...+xn | g1/.../gm}, with \code{x1+...+xn}
   specifying the model for the random effects and \code{g1/.../gm} the
   grouping structure (\code{m} may be equal to 1, in which case no
   \code{/} is required). The random effects formula will be repeated
   for all levels of grouping, in the case of multiple levels of
   grouping; (ii) a list of one-sided formulas of the form
   \code{~x1+...+xn | g}, with possibly different random effects models
   for each grouping level. The order of nesting will be assumed the
   same as the order of the elements in the list; (iii) a one-sided
   formula of the form \code{~x1+...+xn}, or a \code{pdMat} object with
   a formula (i.e. a non-\code{NULL} value for \code{formula(object)}),
   or a list of such formulas or \code{pdMat} objects. In this case, the
   grouping structure formula will be derived from the data used to
   fit the linear mixed-effects model, which should inherit from class
   \code{groupedData}; (iv) a named list of formulas or \code{pdMat}
   objects as in (iii), with the grouping factors as names. The order of
   nesting will be assumed the same as the order of the order of the
   elements in the list; (v) an \code{reStruct} object. See the
   documentation on \code{pdClasses} for a description of the available
   \code{pdMat} classes. Defaults to a formula consisting of the right
   hand side of \code{fixed}.}  
 \item{correlation}{an optional \code{corStruct} object describing the
   within-group correlation structure. See the documentation of
   \code{corClasses} for a description of the available \code{corStruct}
   classes. Defaults to \code{NULL},
   corresponding to no within-group correlations.}
 \item{weights}{an optional \code{varFunc} object or one-sided formula
   describing the within-group heteroscedasticity structure. If given as
   a formula, it is used as the argument to \code{varFixed},
   corresponding to fixed variance weights. See the documentation on
   \code{varClasses} for a description of the available \code{varFunc}
   classes. Defaults to \code{NULL}, corresponding to homocesdatic
   within-group errors.} 
 \item{subset}{an optional expression indicating the subset of the rows of
   \code{data} that should be used in the fit. This can be a logical
   vector, or a numeric vector indicating which observation numbers are
   to be included, or a  character  vector of the row names to be
   included.  All observations are included by default.}
 \item{method}{a character string.  If \code{"REML"} the model is fit by
   maximizing the restricted log-likelihood.  If \code{"ML"} the
   log-likelihood is maximized.  Defaults to \code{"ML"} and in this version is restricted to this default.}
 \item{na.action}{a function that indicates what should happen when the
   data contain \code{NA}s.  The default action (\code{na.fail}) causes
   \code{lme} to print an error message and terminate if there are any
   incomplete observations.}
 \item{control}{a list of control values for the estimation algorithm to
   replace the default values returned by the function \code{lmeControl}.
   Defaults to an empty list.}
 \item{contrasts}{an optional list. See the \code{contrasts.arg}
   of \code{model.matrix.default}.}
 \item{keep.data}{logical: should the \code{data} argument (if supplied
   and a data frame) be saved as part of the model object?}
 \item{start}{an optional numeric vector specifying a set of starting
  variance parameter values for the optimization procedure
  at the unconstrained level.}
 \item{\dots}{ some methods for the generic may require additional arguments.
  None are used in this method.}
}

\description{
  This generic function fits a linear mixed-effects model using an adaptation of the \code{lme.formula}
  from the \code{nlme} package by Pinheiro and Bates.
}

\value{
  An object of class \code{lmeObject} (different from the \code{lme} object generated by
  the \code{lme} function from the \code{nlme} package) representing the linear mixed-effects
  model fit. The components of the fit are the same as the ones given by \code{lme} from
  \code{nlme} package, plus some additional components useful for the implementation of
  methods in the \code{hoalme} package (not intended to be visible to users).
  The \code{print} method generates the same output 
  as with a \code{lme} object, plus information on the optimization used to perform the fit.
  Generic functions \code{summary} and \code{intervals} have methods to show the results of the fit.
}

\references{
 See the \code{lme} documentation from the \code{nlme} package.
}
\author{Sigfrido Iglesias-Gonzalez \email{sigfrido@alumni.utoronto.ca}}
\seealso{
  \code{\link{conslmeObjects}}
  \code{\link{lme}},
  \code{\link{lmeControl}},
  \code{\link{lmeObject}},
}
\examples{
# using 'lmeObject' directly
fit1 <- lmeObject(distance ~ age, data = Orthodont)
summary(fit1)
intervals(fit1)

# fit using an 'lme' object
fit.lme <- lme(distance ~ age, data = Orthodont)
fit2 <- lmeObject(fit.lme) #same as fit1

}
