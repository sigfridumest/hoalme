## file hoalme/R/signedRoot.R, v 0.0-1 2009/11/18
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



##signedRoot methods

signedRoot <-
  function(object, ...)
	UseMethod("signedRoot")	

signedRoot.conslmeObjects<-
  function(object, ...){
  	lmeObj <- attr(object, "lmeObject")
	hypotheses <- object$nat.hypotheses
	nhyp <- object$nhyp
	npar <- object$npar
	npar.fixed <- object$npar.fixed
	logliks <- object$logProfLiks

	signedroot <- sign(hypotheses[ , 1] - hypotheses [ , -1])*sqrt(2*(logliks[ , 1] - logliks[ , -1]))
 	result <- matrix(signedroot, nrow = npar, ncol = nhyp, byrow = F)
	attr(result, "hypotheses") <- hypotheses[ , -1]
	attr(result, "uncons.hypotheses") <- object$hypotheses[ , -1]

	hyp.nams <- dimnames(hypotheses)
	hyp.nams[[ 2 ]] <- hyp.nams[[ 2 ]][ -1 ]
	dimnames(result) <- list(hyp.nams[[ 1 ]], hyp.nams[[ 2 ]])

	attr(result, "nhyp") <- nhyp
	attr(result, "npar") <- npar
	attr(result, "npar.fixed") <- npar.fixed
  	attr(result, "lmeObject") <- lmeObj
	class(result) <- "signedRoot"
	result
}

signedRoot.RStar <-
  function(object, ...){
	attr(object, "signedRoot")
}

signedRoot.lmeObject <-
  function(object, conslmeObjs, ...){
  	if(is.null(conslmeObjs)){
  		if(is.null(object$conslmeObjects)) conslmeObjs <- conslmeObjects(object) #provides a 'conslmeObjects' with default parameters
			else conslmeObjs <- object$conslmeObjects #avoids computations
	}
	signedRoot(conslmeObjs)
}


plot.signedRoot <-
  function(x, natural = TRUE, ...){
	signedroots <- x[ , ]
	hypotheses <- hypotheses(x, natural)
	nams <- dimnames(hypotheses)[[ 1 ]]
	nhyp <- attr(x, "nhyp")
	npar <- attr(x, "npar")
	
	if(npar >= 3) par(mfrow = c(3, 1))
		else par(mfrow = c(2, 1))

	if(npar > 3) pages <- floor(npar/3)
		else pages <- 1
 	for( i in 1:npar){
		if( any(i == ((1:pages)*3 + 1) ) )
			readline("Press ENTER to continue...")
		shyp <- order(hypotheses[i , ])
		plot(x = hypotheses[ i , shyp], y = signedroots[ i , shyp], type = "l", xlab = nams[ i ], ylab = "signed root")
		abline(coef = c(qnorm(.95), 0), col = "orange"); abline(coef = c(-qnorm(.95), 0), col = "orange")
		abline(coef = c(qnorm(.99), 0), col = "red"); abline(coef = c(-qnorm(.99), 0), col = "red")
	}
}

print.signedRoot <-
  function(x, ...){
	print(x[ , ])
}

print.intervals.signedRoot <-
  function(x, ...){
  	lmeObj <- attr(x, "lmeObject")
	fixedEffects.formula <- lmeObj$call$fixed
	reSt <- lmeObj$modelStruct$reStruct
	npar.fixed <- attr(x, "npar.fixed")
	level <- attr(x, "level")
	if(npar.fixed > 0){
		fixed.intervals <- x[1:npar.fixed, ]
		var.intervals <- x[-(1:npar.fixed), ]
	}
	else
		var.intervals <- x
	
	cat("\nApproximate ", 100*level, "% Confidence Intervals using the Likelihood Ratio (r)\n", sep ="")
  	if(npar.fixed > 0){
		cat(" Fixed Effects Parameters:\n")
		cat("  Formula: ", deparse(fixedEffects.formula), "\n")
		print(fixed.intervals)
  	}
   	cat("\n Random Effects Covariance Parameters:\n")
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
  	print(var.intervals[ , ])
}

