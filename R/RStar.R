## file hoalme/R/RStar.R, v 0.0-1 2009/11/18
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

##RStar methods

RStar <-
  function(object, ...)
  	UseMethod("RStar")

RStar.default <-
  function(object, ...){
	attr(object, "RStar")
}

RStar.lmeObject <-
  function(object, conslmeObjs, delta = 0.7, ...){
  	sigRoot <- signedRoot(conslmeObjs)
 	npar <- conslmeObjs$npar
 	npar.fixed <- conslmeObjs$npar.fixed
 	conslmelogLiks <- logLik(conslmeObjs, object)
  	QStat <- QStatistic(conslmelogLiks, object)
  	hypotheses <- hypotheses(QStat)
  	n.hyp <- ncol(hypotheses)
  	Q.on.R <- QStat/sigRoot
  	log.Q.on.R <- log(Q.on.R)
  	corrTerm <- log.Q.on.R/sigRoot
  	RStar <- sigRoot + corrTerm
  	RStar.orig <- RStar
  	waldStat <- attr(QStat, "WaldStatistic")
   	for(i in 1:npar){
	  	points.in <- abs(sigRoot[ i , ]) > delta
	  	points.out <- abs(sigRoot[ i , ]) <= delta
	  	if(sum(ifelse(points.in, 1, 0)) < 4){
	  		points.in <- rep(F, n.hyp)
	  		points.in[c(1, 2, n.hyp - 1, n.hyp)] <- T
	  	}
		spline.obj <- smooth.spline(y = RStar[ i , points.in ], x = sigRoot[ i , points.in ], all.knots = TRUE)
		RStar[ i, points.out] <- predict(spline.obj, x = sigRoot[ i , points.out])$y
	}
		
	class(RStar) <- NULL
	attr(RStar, "Q.on.R") <- Q.on.R
	attr(RStar, "logQonR") <- log.Q.on.R	
	attr(RStar, "corrTerm") <- corrTerm
	attr(RStar, "signedRoot") <- sigRoot
	attr(RStar, "QStatistic") <- QStat
	attr(RStar, "WaldStatistic") <- waldStat
	attr(RStar, "RStar.orig") <- RStar.orig
	attr(RStar, "hypotheses") <- hypotheses
	attr(RStar, "npar.fixed") <- attr(conslmeObjs, "npar.fixed")
	attr(RStar, "lmeObject") <- object
	class(RStar) <- "RStar"
	RStar
}

RStar.hoalme <-
  function(object, ...){
	object$RStar
}

print.RStar <-
  function(x, ...){
 	print(x[ , ])
}

plot.RStar <-
  function(x, what = c("RStar", "corrTerm", "logQonR", "QStatistic"), single = FALSE, compare = FALSE, legend = TRUE, ...){
	quantity <- match.arg(what)
	hypotheses <- attr(x, "hypotheses")
	npar <- attr(x, "npar")
	nams <- dimnames(hypotheses)[[ 1 ]]
	if(quantity == "RStar")
		Qty <- x[ , ]
		else
			Qty <- attr(x, quantity)

 	if(quantity == "RStar" && compare){
 		par(mfrow = c(1,1))
		sigRoots <- signedRoot(x)
		for(i in 1:npar){
			maximo <- max(c(sigRoots[ i , ], Qty[ i , ]))
			minimo <- min(c(sigRoots[ i , ], Qty[ i , ]))	
			plot(y = Qty[ i , ], x = hypotheses[ i, ], type = "l", xlab = nams[ i ], ylab = quantity, ylim = c(minimo, maximo))
				abline(a = 1.96, b = 0, col = "orange");abline(a = -1.96, b = 0, col = "orange");abline(a = 2.32, b = 0, col = "red");abline(a = -2.32, b = 0, col = "red");
			lines(y = sigRoots[ i , ], x = hypotheses[ i, ], type = "l", col = "blue")
			if(legend)
				legend(x = "r", c("modified ratio r*", "likelihood ratio r"), col = c("black", "blue"), lty = c(1,1), bty = "n")
			if(i < npar)
				readline("Press ENTER to continue...")
		}	
	}
	else{
		if(!single){
			if(npar >= 3)
				par(mfrow = c(3, 1))
				else
					par(mfrow = c(npar, 1))
					
			if(npar > 3) pages <- floor(npar/3)
				else pages <- 1
			for(i in 1:npar){
				if( any(i == ((1:pages)*3 + 1) ) )
					readline("Press ENTER to continue...")
				plot(y = Qty[ i , ], x = hypotheses[ i, ], type = "l", xlab = nams[ i ], ylab = quantity)
				if(quantity == "RStar"){
					abline(a = 1.96, b = 0, col = "orange");abline(a = -1.96, b = 0, col = "orange");abline(a = 2.32, b = 0, col = "red");abline(a = -2.32, b = 0, col = "red")
				}
			}
		}
		else{
			par(mfrow = c(1,1))
			for(i in 1:npar){
				plot(y = Qty[ i , ], x = hypotheses[ i, ], type = "l", xlab = nams[ i ], ylab = quantity)
				if(quantity == "RStar"){
					abline(a = 1.96, b = 0, col = "orange");abline(a = -1.96, b = 0, col = "orange");abline(a = 2.32, b = 0, col = "red");abline(a = -2.32, b = 0, col = "red");
				}
				if(i < npar)
					readline("Press ENTER to continue...")
			}
		}
	}
}

print.intervals.RStar <-
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
	cat("\nApproximate ", 100*level, "% Confidence Intervals using the Modified Ratio (r*)\n", sep ="")
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

update.RStar <-
  function(object, delta, ...){
	sigRoot  <- attr(object, "signedRoot")
	RStar.orig <- attr(object, "RStar.orig")
   	npar <- nrow(object)
   	nhyp <- ncol(object)
   	for(i in 1:npar){
	  	points.in <- abs(sigRoot[ i , ]) > delta
	  	points.out <- abs(sigRoot[ i , ]) <= delta
	  	if(sum(ifelse(points.in, 1, 0)) < 4){
	  		points.in <- rep(F, nhyp)
	  		points.in[c(1, 2, nhyp - 1, nhyp)] <- T
	  	}
		spline.obj <- smooth.spline(y = RStar.orig[ i , points.in ], x = sigRoot[ i , points.in ], all.knots = TRUE)
		RStar.orig[ i , points.out] <- predict(spline.obj, x = sigRoot[ i , points.out])$y
	}
	attributes(RStar.orig) <- attributes(object)
	RStar.orig
}


#RStar.summary.hoalme <-
#   function(object, ...){
# 	attr(object, "RStar")
# }
