## file hoalme/R/QStatistic.R, v 0.0-1 2009/11/18
##
##  Copyright (C) 2009 Sigfrido Iglesias-Gonzlez
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


## QStatistic methods

QStatistic <- 
 function(object, ...)
	UseMethod("QStatistic")

QStatistic.lmeLogLik <-
  function(object, consLik, hyp.names, ...){
 	secondDer.full <- object$secondDeriv
	secondDer.det <- object$secondDeriv.det
	nuissecondDer.dets <- object$nuissecondDeriv.dets
	ssDer.full <- object$SSDeriv
	dirDer.full <- ssDer.full$dirDeriv
	fact <- ssDer.full$fact
	secondDer.cons <- consLik$secondDeriv
	ssDer.cons <- consLik$SSDeriv
	hypothesis <- attr(consLik, "hypothesis")
	dirDer.cons <- ssDer.cons$dirDeriv
	MixedDer.cons <- ssDer.cons$MixedDeriv
	param.names <- dimnames(secondDer.full)[[ 1 ]]
	hypothesis.name <- hypothesis
	hypothesis <- match(hypothesis, hyp.names)
	if(length(param.names) == 2 + length(hyp.names))
		hypothesis <- hypothesis + 1
	profObsInfo <- -secondDer.cons[ -hypothesis, -hypothesis]
	ObsInfo.det <- secondDer.det #comes in absolute value
	eigen.profObsInfo <- eigen(profObsInfo, T, T)$values
	profObsInfo.det <- abs(prod(eigen.profObsInfo))

 	if(any(eigen(profObsInfo)$values < 0)){
		warning("Non positive-definite Profile Observed Information matrix.")
		cat("\n hypothesis:", hypothesis.name, "\n")
		pd <- F
	}
	else pd <- T
	nuisObsInfo.det <- nuissecondDer.dets[hypothesis] #comes in absolute value
	j.p <- ObsInfo.det/nuisObsInfo.det

	deviation <- dirDer.full - dirDer.cons
	MixedDer.cons <- MixedDer.cons[ , -hypothesis]
	conspar <- c(consLik$fixedPar, consLik$varPar)
	par <- c(object$fixedPar, object$varPar)
	psi <- conspar[attr(consLik, "hypothesis")]
	MLE <- par[attr(consLik, "hypothesis")]
	waldStat <- (MLE - psi)*sqrt(j.p)

	##Decomposition
	deviation.mle <- fact%*%deviation
	MixedDer.cons.mle <- fact%*%MixedDer.cons
	Corr <- abs(det(MixedDer.cons.mle[ -hypothesis,]))/sqrt(nuisObsInfo.det*profObsInfo.det)
#	Corr <- sqrt(nuisObsInfo.det*profObsInfo.det)/abs(det(MixedDer.cons.mle[ -hypothesis,]))

 	l.p <- deviation.mle[ hypothesis , ] -  MixedDer.cons.mle[ hypothesis , ]%*%solve(MixedDer.cons.mle[ -hypothesis , ])%*%deviation.mle[-hypothesis , ]
 	l.p <- l.p/sqrt(j.p)
	QStat <- Corr*l.p
# 	QStat <- l.p/Corr

 	attr(QStat, "NP") <- Corr
 	attr(QStat, "INF") <- l.p
	attr(QStat, "hypothesis") <- conspar[attr(consLik, "hypothesis")]
	attr(QStat, "WaldStatistic") <- waldStat
	attr(QStat, "PD") <- pd
	class(QStat) <- "QStatistic"
	QStat
}

QStatistic.conslmeLogLik <-
  function(object, logLik.full, hyp.names, ...){
	QStatistic(logLik.full, object, hyp.names)
}

QStatistic.conslmeLogLiks <-
  function(object, lmeObj, ...){
	hyp.names <- dimnames(hypotheses(object))[[1]]
	logLik.full <- logLik(lmeObj)
	result <- lapply(object, QStatistic, logLik.full, hyp.names)
	waldStat <- unlist(lapply(result, function(x){
			if(is.na(x)) return(NA)
				else attr(x, "WaldStatistic")
		}))
	pd <- unlist(lapply(result, function(x){
					 attr(x, "PD")
				}))
	NP <- unlist(lapply(result, function(x){
			if(is.na(x)) return(NA)
				else attr(x, "NP")
		}))
	INF <- unlist(lapply(result, function(x){
			if(is.na(x)) return(NA)
				else attr(x, "INF")
		}))

	result <- unlist(result)

	npar <- attr(object, "npar")
	hypotheses <- attr(object, "hypotheses")
	nhyp <- ncol(hypotheses)
	result <- matrix(result, nrow = npar, ncol = nhyp, byrow = T)
	waldStat <- matrix(waldStat, nrow = npar, ncol = nhyp, byrow = T)
	pd <- matrix(pd, nrow = npar, ncol = nhyp, byrow = T)
	NP <- matrix(NP, nrow = npar, ncol = nhyp, byrow = T)
	INF <- matrix(INF, nrow = npar, ncol = nhyp, byrow = T)

	for(i in 1:npar){
		pd.indicator <- pd[ i , ]
		if( any(!pd.indicator) ){
			Q.fine <- result[ i , pd.indicator]
			smooth.obj <- smooth.spline(x = hypotheses[ i , pd.indicator], y = Q.fine, all.knots = T)
			Q.fixed <- predict( smooth.obj, x = hypotheses[ i , !pd.indicator] )$y
			result[ i , !pd.indicator] <- Q.fixed
		}
	}

	dimnames(result) <- dimnames(hypotheses)
	attr(result, "npar") <- attr(object, "npar")
	attr(result, "MLE") <- attr(object, "MLE")
	attr(result, "hypotheses") <- hypotheses
	attr(result, "NP") <- NP
	attr(result, "INF") <- INF
	attr(result, "WaldStatistic") <- waldStat
	attr(result, "lmeObject") <- attr(object, "lmeObject")

	class(result) <- "QStatistic"
	result
}

QStatistic.lmeObject <-
  function(object, conslmeObjs, ...){
	QStatistic(logLik(conslmeObjs), object)
}

QStatistic.RStar <-
  function(object, ...){
	attr(object, "QStatistic")
}

plot.QStatistic <-
  function(x, what = c( "QStatistic", "NP", "INF"), log = TRUE, ...){
	quantity <- match.arg(what)
 	npar <- attr(x, "npar")
	if(npar >= 3) par(mfrow = c(3, 1))
		else par(mfrow = c(2, 1))

	hypotheses <- attr(x, "hypotheses")
	nams <- dimnames(hypotheses)
	n.hyp <- ncol(hypotheses)

	if(quantity == "QStatistic")
		Qty <- x[ , ]
		else
			Qty <- attr(x, quantity)

	if(npar > 3) pages <- floor(npar/3)
		else pages <- 1
 	for( i in 1:npar){
		if( any(i == ((1:pages)*3 + 1) ) )
			readline("Press ENTER to continue...")
		if( quantity == "NP" && log )
			plot(x = hypotheses[ i, ], y = log(Qty[ i , ]), type = "l", ylab = paste("(log)", quantity, sep = ""), xlab = nams[[ 1 ]][ i ], ...)
			else
				plot(x = hypotheses[ i, ], y = Qty[ i , ], type = "l", ylab = quantity, xlab = nams[[ 1 ]][ i ], ...)
	}
}

print.QStatistic <-
  function(x, ...){
	print(x[ , ])
}

