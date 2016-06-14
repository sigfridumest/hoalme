## file hoalme/R/hypotheses.R, v 0.0-1 2009/11/18
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


##hypotheses methods

hypotheses <-
  function(object, ...)
	UseMethod("hypotheses")

hypotheses.default <-
  function(object, ...){
	attr(object, "hypotheses")
}

hypotheses.pdMat <-
  function(object, lmeObj, nam, scale = 2, precision = 10, ...){
	##getting unconstrained parameter values at optimization level
	uns.Struct <- solve(object)
	pars <- mle <- coef(uns.Struct)
	npar <- length(pars)

	## getting approximate covariance matrix of estimates
	apVar <- unconsVar(lmeObj)
	##selecting entries from asymptotic cov. matrix
	apVar.names<- dimnames(apVar)[[ 1 ]]
	which.are <- match(nam, apVar.names) 
	diag.par.Var <- diag(apVar)

	diag.par.Var <- diag.par.Var[which.are] #selecting the relevant entries
	if(any(diag.par.Var < 0)){
		diag.par.Var[diag.par.Var <= 0] <- 1e-6
 		warning("Non-positive semidefinite approximate covariance matrix")
	}

	lims <- matrix(0, ncol = 2, nrow = npar)

	##constructing limits
	if(!is.list(scale)) scale <- rep(list(rep(scale, 2)), npar)
	for(i in seq_along(scale)){
		uplim <- pars[ i ] + scale[[ i ]][ 2 ]*sqrt(diag.par.Var[ i ])
		lowlim <- pars[ i ] - scale[[ i ]][ 1 ]*sqrt(diag.par.Var[ i ])
		lims[ i , ] <- c(lowlim, uplim)
	}

	## hypotheses in unconstrained form, departuring from the mle
	hypotheses <- list()
	scales <- list()

	precision2 <- ( precision/2 +1)
	for( i in 1:npar){
		hyp_seq1 <- seq(lims[ i, 1 ], mle[ i ], length = precision2)[ -precision2 ]
		hyp_seq2 <- seq(mle[ i ], lims[ i , 2 ], length = precision2)[ -1 ]
		hyp_seq <- c(hyp_seq1, hyp_seq2)
		hypotheses[[ 2*(i - 1) + 1 ]] <- sort(hyp_seq[ hyp_seq < mle[ i ] ], T) #decreasing
		attr(hypotheses[[ 2*(i - 1) + 1 ]], "name") <- nam[ i ]
		attr(hypotheses[[ 2*(i - 1) + 1 ]], "scales") <- (c(hypotheses[[ 2*(i - 1) + 1 ]]) - pars[ i ] )/sqrt(diag.par.Var[ i ])
		hypotheses[[ 2*(i - 1) + 2 ]] <- sort(hyp_seq[ hyp_seq > mle[ i ] ], F)
		attr(hypotheses[[ 2*(i - 1) + 2 ]], "name") <- nam[i]
		attr(hypotheses[[ 2*(i - 1) + 2 ]], "scales") <- (c(hypotheses[[ 2*(i - 1) + 2 ]]) - pars[ i ] )/sqrt(diag.par.Var[ i ])
	}
	attr(hypotheses, "mle") <- mle
	hypotheses
}

hypotheses.reStruct <-
  function(object, lmeObj, scale = 2, precision = 10, ...){
	## Constructing lists of hypotheses for each 'Positive definite matrix structure' present at each group in 'reStruct'
	namStruct <- names(object) #levels names
	result <- list()


	for( i in seq_along(object) ){ #for each level

		if(inherits(object[[ i ]], "pdBlocked")){
			result[[ i ]] <- list()
			nblocks <- length(object[[ i ]]) #number of blocks at the i-th level
			for(j in seq_along(object[[ i ]])){ #blocks
 				npar <- length(coef(object[[ i ]][[ j ]])) #parameters in j-th block
# 				if(npar > 1)
# 					nam <- paste("reStruct.", namStruct[ i ], 1:nblocks, sep = "") #reStruct.ith-levelname.blocknumber
# 					paste(nam, 1
# 				else
					nam <- paste("reStruct.", namStruct[ i ], 1:nblocks, sep = "") #reStruct.ith-levelname.blocknumber

				St <- object[[ i ]][[ j ]] #structure per block
				if( !is.list(scale) && length(scale) == 1 )
					result[[ i ]][[ j ]] <- hypotheses(St, lmeObj, nam[ j ], scale, precision)
					else
						result[[ i ]][[ j ]] <- hypotheses(St, lmeObj, nam[ j ], scale = scale[[ i ]][[ j ]], precision)
			}
		}
		else{
			npar <- length(coef(object[[ i ]]))
			if(npar == 1) nam <- paste("reStruct.", namStruct[ i ], sep = "")
				else	nam <- paste("reStruct.", namStruct[ i ], 1:npar, sep = "")
			## getting hypotheses for each level
			St <- object[[ i ]] #ith level
# 			for(j in 1:npar){
				if( !is.list(scale) && length(scale) == 1 ) 
					result[[ i ]] <- hypotheses(St, lmeObj, nam, scale, precision)
					else{
# 						print(scale[[ i ]])
						result[[ i ]] <- hypotheses(St, lmeObj, nam, scale = scale[[ i ]], precision)
					}
# 				}#
			}
	} #for
	result
}

hypotheses.conslmeObjects <-
 function(object, natural = TRUE, ...){
	if(natural) return(object$nat.hypotheses)
		else return(object$hypotheses)
}

hypotheses.signedRoot <-
  function(object, natural = TRUE, ...){
	if(natural)
		result <- attr(object, "hypotheses")
		else
			result <- attr(object, "uncons.hypotheses")
	result
}
