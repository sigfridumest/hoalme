## file hoalme/R/ConsFixedEffects.R, v 0.0-1 2009/11/18
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



##consFixedEffects methods

consFixedEffects <- 
  function(object, ...)
	UseMethod("consFixedEffects")

consFixedEffects.lmeObject <-
  function(object, scale = 2, precision = 10, control = list(returnObject = T)){
 	##Extracting relevant objects
 	data.tmp <- object$evalEnv$data #data.tmp recovers the original data set
	idf <- is.data.frame(data.tmp)
 	call <- object$call
 	call$control <- control
 	orig.formula <- as.formula(call$fixed)
	resp <- as.character(orig.formula[[ 2 ]])
 	X <- model.matrix(orig.formula, data.tmp)  #recovers original design matrix, in original ordering,just the way it entered into the lme routine
 	coeffFixed <- object$coefficients$fixed
 	varFix <- object$varFix
 	nams <-  dimnames(object$evalEnv$X)[[ 2 ]] ##Need these, and no other names##
 	npar <- length(nams)
 	std.errors <- sqrt(diag(varFix))
	lims <- matrix(0, nrow = npar, ncol = 2)
	hypotheses <- list()

	##modifying the call directly, to provide a data 'environment' where lmeobject can evaluate all the arguments
	updated.formula <- as.formula(paste("Response ~ pred - 1")) #no intecept, as it is included in the model.matrix (or not)
	call$fixed <- updated.formula
  	if(!is.null(call[[ "data" ]])) call[[ "data" ]] <- NULL
	call$data <- quote(data.tmp)
	
	##constructing the hypotheses
	precision2 <- ( precision/2 +1)
 	for( i in 1:npar){
 		std.err <- std.errors[ i ]
 		f.e <- coeffFixed[ i ]
 		if(!is.list(scale) && length(scale) == 1){
	 		lims[ i , 2] <- f.e + scale*std.err
 			lims[ i , 1] <- f.e - scale*std.err
 		}
 		else{
	 		lims[ i , 2] <- f.e + scale[[ i ]][ 2 ]*std.err
 			lims[ i , 1] <- f.e - scale[[ i ]][ 1 ]*std.err
 		}
		hyp_seq1 <- seq(lims[ i,1 ], f.e, length = precision2)[ -precision2 ]
		hyp_seq2 <- seq(f.e, lims[ i , 2 ], , length = precision2)[ -1 ]
		hyp_seq <- c(hyp_seq1, hyp_seq2)
		hypotheses[[ 2*(i - 1) + 1 ]] <- sort(hyp_seq[ hyp_seq < f.e ], T) #decreasing
		attr(hypotheses[[ 2*(i - 1) + 1 ]], "name") <- nams[ i ]
		hypotheses[[ 2*(i - 1) + 2 ]] <- hyp_seq[ hyp_seq > f.e ]
		attr(hypotheses[[ 2*(i - 1) + 2 ]], "name") <- nams[ i ]
 	}

	##
	##Constrained optimization
	##
	results <- list()
	l <- 1
	for(i in 1:npar){
		pred <- X[ , - i]
		for(j in 1:2){
			hyp <- hypotheses[[ (i -1)*2 + j ]]
			for( k in seq_along(hyp)){
				constraint <- X[ , i ]*hyp[ k ]
				data.tmp$Response <- data.tmp[ , resp ] - constraint
				data.tmp$pred <- pred
				results[[ l ]] <- eval(call) #couldn't use 'update.formula', because some character strings don't convert to terms...!!
				full <- list(X = object$evalEnv$X, y = object$evalEnv$y)
				fixedPars <- coeffFixed
				fixedPars[ i ] <- hyp[k]
				fixedPars[ -i ] <- results[[ l ]]$coefficients$fixed
				results[[ l ]]$fixedPars <- fixedPars
				hypothesis <- list(constraint = constraint, hypothesis = hyp[ k ])
				results[[ l ]]$full <- full
				results[[ l ]]$hypothesis <- hypothesis
				attr(results[[ l ]], "hypothesis") <- nams[ i ]
				class(results[[ l ]]) <- c("consFixedEffects", "conslmeObject", "lmeObject", "Object")
				l <- l +1
			}
		}
	}
	
	##ordering the objects generated
	hypotheses <- matrix(unlist(hypotheses), ncol = precision, byrow = T)
	ord.hyp <- apply(hypotheses, 1, sort.int, index.return = T)
	ix <- matrix(unlist(lapply(ord.hyp, function(x) x[[ 2 ]])), ncol = precision, byrow = T)
	hypotheses <- matrix(unlist(lapply(ord.hyp, function(x) x[[ 1 ]])), ncol = precision, byrow = T)	
	if(npar == 1)
		index <- c(ix)
		else
			index <- c( t(ix + precision*c( 0,1:(npar-1) )) )
	results <- results[index]
	attr(results, "hypotheses") <- hypotheses
	results
}

model.matrix.consFixedEffects <-
  function(object, ...){
  	##Objects actually generated during the constrained optimization; they are (y & X) the constrained versions (Response and pred)
	Z <- object$evalEnv$Z
	y <- object$evalEnv$y
	X <- object$evalEnv$X
# 	constrain <- object$hypothesis$constrain
	parFix <- matrix(object$coefficients$fixed, ncol = 1)
	r <- y - X%*%parFix
	
	dims <- object$dims
	ncols <- dims$ncol #columns: Z 1st level ... highest level, matrix X, response
	N <- dims$N #total number of observations
	Q <- dims$Q #number of groups

	G <- rev(object$groups) #factors for grouping; come froem highest to lowest: reversing
# 	GZXy <- data.frame(G, Z, X, y)
	X <- object$full$X
	y <- object$full$y
	GZXyr <- cbind(G, Z, X, y, r)

	##indexing components	
	ends <- cumsum(ncols[1:Q]) 
	starts <- c(0, ends[Q-1]) + 1
	attr(GZXyr, "Zcols") <- list(cols = Q + 1:ends[Q], starts = starts, ends = ends)
	attr(GZXyr, "gcols") <- 1:Q
	attr(GZXyr, "Xcols") <- (Q + ends[Q] + 1):(ncol(GZXyr) - 2)
	attr(GZXyr, "ycol") <- ncol(GZXyr) - 1
	attr(GZXyr, "rescol") <- ncol(GZXyr)
	attr(GZXyr, "Q") <- Q
	class(GZXyr) <- c("grps.model.frame", "data.frame")
	GZXyr
}
