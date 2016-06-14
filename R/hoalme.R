## file hoalme/R/hoalme.R, v 0.0-1 2009/11/18
##
##  Copyright (C) 2009 Sigfrido Iglesias-Gonzalez
##
##  This file is part of the "hoalme" package for R.  This program is 
##  free software; you can redistribute it and/or modify it under the 
##  terms of the GNU General Public License as published by the Free 
##  Software Foundation; either version 2 of the License, or (at your 
##  option) any later version.
##
##  This library is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with this program; if not, write to the Free Software
##  Foundation, Inc., 59 Temple Place, Suite 330, Boston, 
##  MA 02111-1307 USA or look up the web page 
##  http://www.gnu.org/copyleft/gpl.html.
##
##  Please send any comments, suggestions or errors found to:
##  Sigfrido Iglesias-Gonzalez
##  Email: sigfrido@alumni.utoronto.ca.  



##hoalme methods

hoalme <-
  function(object, ...)
	UseMethod("hoalme")

hoalme.lmeObject <-
  function(object, conslmeObjs, ...){
	result <- list()
	RStar <- RStar(object, conslmeObjs, ...)
	sigRoot <- signedRoot(RStar)
	result$signedRoot <- sigRoot
	result$RStar <- RStar
	attr(result, "lmeObject") <- object
# 	attr(result, "nhyp") <- ncol(sigRoot)
	attr(result, "npar") <- nrow(sigRoot)
	attr(result, "npar.fixed") <- attr(conslmeObjs, "npar.fixed")
	class(result) <- "hoalme"
	result
}

hoalme.RStar <-
  function(object, ...){
	lmeObj <- attr(object, "lmeObject")
	conslmeObjs <- lmeObj$conslmeObjects
	result <- list()
	sigRoot <- signedRoot(object)
	result$signedRoot <- sigRoot
	result$RStar <- object
	attr(result, "lmeObject") <- lmeObj
	attr(result, "nhyp") <- ncol(sigRoot)
	attr(result, "npar") <- nrow(sigRoot)
	attr(result, "npar.fixed") <- attr(conslmeObjs, "npar.fixed")
	class(result) <- "hoalme"
	result
}

print.hoalme <-
  function(x, ...){
	lmeObj <- attr(x, "lmeObject")
	fixedEffects.fomula <- lmeObj$call$fixed
 	npar.fixed <- attr(x, "npar.fixed")
	reSt <- lmeObj$modelStruct$reStruct
		
	cat("\n   ** Higher Order Asymptotics ** \n")
	cat( "           applied to the\n")
	cat("     Linear Mixed Effects Model\n")
	cat("              using\n")
	cat(" ** Maximum Likelihood Estimation **\n\n" )
	##Model Information
	cat("Model Information:\n\n")
	cat(" Fixed Effect Formula: \n")
	cat(" ", deparse(fixedEffects.fomula), "\n")
	cat(" Random Effects Formula:\n")
	printInf(reSt)
	dd <- lmeObj$dims	
	cat("\n Total Number of Observations:", dd[["N"]])
	cat("\n Number of Groups (Independent Vectors): ")
	Ngrps <- dd$ngrps[1:dd$Q]
	if ((lNgrps <- length(Ngrps)) == 1) {	# single nesting
		cat(Ngrps,"\n")
	} else {				# multiple nesting
		sNgrps <- 1:lNgrps
		aux <- rep(names(Ngrps), sNgrps)
		aux <- split(aux, array(rep(sNgrps, lNgrps),
		c(lNgrps, lNgrps))[!lower.tri(diag(lNgrps))])
		names(Ngrps) <- unlist(lapply(aux, paste, collapse = " %in% "))
		cat("\n")
		print(rev(Ngrps))
	}
	param.names <- dimnames(hypotheses(lmeObj$conslmeObjects))[[ 1 ]]
	nhyp <- attr(x, "nhyp")
#  print(nhyp)
	npar <- attr(x, "npar")
	cat("\n\nHigher Order Asysmptotic Inference based on", nhyp, "points is available on the following parameters\n")
	if(npar.fixed > 0){
		cat(" Fixed Effects:\n  ", param.names[ 1:npar.fixed], " \n")
		cat(" Random Effects Covariance parameters:\n  ", param.names[ -(1:npar.fixed)], "\n")
	}
	else
		cat(" Random Effects Covariance parameters:\n  ", param.names, "\n")
	cat("\n(use 'summary' method)\n\n")

##	information from lmeObject
##	Number of points being tested (put it on the constructor)
##	Higher order inference methods available: RStar so far
##
}

print.summary.hoalme <-
  function(x, ...){
  RStar <- RStar(x)
  lmeObj <- attr(RStar, "lmeObject")
  add.tests <- attr(x, "add.tests")
  reSt <- lmeObj$modelStruct$reStruct
  level <- attr(x, "level")
  fixedEffects.fomula <- lmeObj$call$fixed
  npar.fixed <- attr(x, "npar.fixed")
  cat("** Higher Order Inference for the Mixed Linear Model **\n")
	if(npar.fixed > 0 ){
		##fixed effects : Tests
		cat("\nFixed Effects two-sided hypothesis tests:\n\n")
		cat(" Formula: ", deparse(fixedEffects.fomula), "\n")
		cat("\n                 Modified Ratio (r*) based\n" )
		print(x$tests.RStar)
		if(add.tests == "all"){
		cat("\n                   Signed Ratio (r) based\n" )
			print(x$tests.signedRoot)
		}
		##Confidence Intervals
		cat("\n\nFixed Effects ", level*100, "% Confidence Intervals :\n", sep = "")
		cat("\n                 Modified Ratio (r*) based\n" )
		print(x$intervals.RStar[1:npar.fixed , ])
		if(add.tests == "all"){
			cat("\n                   Signed Ratio (r) based\n" )
			print(x$intervals.signedRoot[1:npar.fixed , ])
		}
		
		## Random Effects covariance paramaters
		cat("\n\nRandom Effects Covariance Parameters:\n\n")
		Levels <- names(reSt)
		for(i in seq_along(reSt)){
			St <- reSt[[ i ]]
			if (!is.list(St)) {
				if (!(is.null(form <- attr(St, "formula")))) {
				cat(paste(" Formula: "))
				if (inherits(form, "formula")) {
					cat(deparse(form))
					if (!is.null(Levels)) { cat( paste( " |", Levels[ i ] ) ) }
				} else {
					if (length(form) == 1) {
					cat(deparse(form[[1]]))
					if (!is.null(Levels)) { cat( paste( " |", Levels[ i ] ) ) }
					} else {
					cat(deparse(lapply(form,
										function(el) as.name(deparse(el)))))
					cat("\n Level:", Levels[ i ])
					}
				}
				cat( "\n" )
				}
				cat(paste("  Structure: ", attr(St, "class")[1], "\n", sep = ""))
			}
		}
		cat("\n ", level*100, "% Confidence Intervals:\n", sep = "")
		cat("\n                 Modified Ratio (r*) based\n" )
		print(x$intervals.RStar[ -(1:npar.fixed) , ])
		if(add.tests == "all"){
			cat("\n                   Signed Ratio (r) based\n" )
			print(x$intervals.signedRoot[ -(1:npar.fixed) , ])
		}
	} #there are fixed effects to be tested
	else{ #there are no fixed effects to be tested
		## Random Effects covariance paramaters
		cat("\n\n Random Effects Covariance Parameters\n\n", sep = "")
		Levels <- names(reSt)
		for(i in seq_along(reSt)){
			x <- reSt[[ i ]]
			if (!is.list(x)) {
				if (!(is.null(form <- attr(x, "formula")))) {
				cat(paste(" Formula: "))
				if (inherits(form, "formula")) {
					cat(deparse(form))
					if (!is.null(Levels)) { cat( paste( " |", Levels[ i ] ) ) }
				} else {
					if (length(form) == 1) {
					cat(deparse(form[[1]]))
					if (!is.null(Levels)) { cat( paste( " |", Levels[ i ] ) ) }
					} else {
					cat(deparse(lapply(form,
										function(el) as.name(deparse(el)))))
					cat("\n Level:", Levels[ i ])
					}
				}
				cat( "\n" )
				}
				cat(paste("  Structure: ", attr(x, "class")[1], "\n", sep = ""))
			}
		}
		cat("\n                 Modified Ratio (r*) based\n" )
		print(x$intervals.RStar)
		if(add.tests == "all"){
			cat("\n                   Signed Ratio (r) based\n" )
			print(x$intervals.signedRoot)
		}
	}
}

summary.hoalme <-
  function(object, level = .95, add.tests = c("all", "none", "r", "r~"), sig = 5, ...) {
  ## fixed effects r* & r tests
	lmeObj <- attr(object, "lmeObject")
	RStar <- object$RStar
	sigRoot <- object$signedRoot
	npar.fixed <- attr(object, "npar.fixed")
	hyp <- hypotheses(RStar)
	if(npar.fixed > 0){
		tests.RStar <- tests.sigRoot <- matrix(0, ncol = 3, nrow= npar.fixed)
		for(i in 1:npar.fixed){
			spline.obj <- smooth.spline(y = RStar[ i , ], x = hyp[ i , ], all.knots = TRUE, spar = 1e-7)
			tests.RStar[ i , 2] <- rs <- predict(spline.obj, x = 0)$y
			tests.RStar[ i , 3] <- round( (1 - pnorm(abs(rs)))*2 , sig)
			spline.obj <- smooth.spline(y = sigRoot[ i , ], x = hyp[ i , ], all.knots = TRUE, spar = 1e-7)
			tests.sigRoot[ i , 2 ] <- sr <- predict(spline.obj, x = 0)$y
			tests.sigRoot[ i , 3 ] <- round( (1 - pnorm(abs(sr)))*2 , sig)
		}
		
		dimnames(tests.sigRoot) <- dimnames(tests.RStar) <- list(dimnames(hyp)[[1 ]][ 1:npar.fixed ], c("estimate", "r*", "p-value"))
		dimnames(tests.sigRoot)[[ 2 ]][ 2 ] <- "r"
		intervals.RStar <- intervals(RStar, level)
		intervals.sigRoot <- intervals(sigRoot, level)
		tests.sigRoot[ , 1] <- tests.RStar[ , 1] <- intervals.RStar[1:npar.fixed, 2]
		result <- list(intervals.RStar = intervals.RStar, intervals.signedRoot = intervals.sigRoot, tests.RStar = tests.RStar, tests.signedRoot = tests.sigRoot)
	}
	else
		result <- list(intervals.RStar = intervals(RStar), intervals.signedRoot = intervals(sigRoot))
	add.tests <- match.arg(add.tests)
	attr(result, "RStar") <- RStar
	attr(result, "npar.fixed") <- npar.fixed
# print(npar.fixed)
	attr(result, "level") <- level
	attr(result, "add.tests") <- add.tests
	class(result) <- "summary.hoalme"
	result
}

printInf <-
  function(object, Levels = names(object)){
	for(i in seq_along(object)){
	x <- object[[ i ]]
	if (!is.list(x)) {
		if (!(is.null(form <- attr(x, "formula")))) {
		if (inherits(form, "formula")) {
			cat(" ", deparse(form))
			if (!is.null(Levels)) { cat( paste( " |", Levels[ i ] ) ) }
		} else {
			if (length(form) == 1) {
			cat(" ", deparse(form[[1]]))
			if (!is.null(Levels)) { cat( paste( " |", Levels[ i ] ) ) }
			} else {
			cat(" ", deparse(lapply(form,
								function(el) as.name(deparse(el)))))
			cat("\n Level:", Levels[ i ])
			}
		}
		cat( "\n" )
		}
		cat(paste("  Structure: ", attr(x, "class")[1], "\n", sep = ""))
	}
	else {	# composite structure
		cat(paste(" Composite Structure: ", class(x)[ 1 ], "\n", sep =""))
		for(j in seq_along(x)){
			x. <- x[[ j ]]
			if (!(is.null(form <- attr(x., "formula")))) {
				if (inherits(form, "formula")) {
					cat(" ", deparse(form))
					if (!is.null(Levels)) { cat( paste( " |", Levels[ i ] ) ) }
				} else {
					if (length(form) == 1) {
					cat(" ", deparse(form[[1]]))
					if (!is.null(Levels)) { cat( paste( " |", Levels[ i ] ) ) }
					} else {
					cat(" ", deparse(lapply(form,
										function(el) as.name(deparse(el)))))
					cat("\n Level:", Levels[ i ])
					}
				}
			cat( "\n" )
			}
			cat(paste("  Structure: ", attr(x., "class")[1], "\n", sep = ""))
		}
	}
	}#for
}

