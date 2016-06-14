## file hoalme/R/pdDeriv.R, v 0.0-1 2009/11/18
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



##pdDeriv methods

pdDeriv <-
  function(object, ...)
	UseMethod("pdDeriv")

pdDeriv.pdIdent <-
  function(object, Z, ...){
	result <- tcrossprod(as.matrix(Z))
	result
}

pdDeriv.pdSymm <-
  function(object, Z, ...){
	Z <- as.matrix(Z)
	m <- ncol(Z)
	d <- m*(m+1)*.5
	n <- nrow(Z)
	l <- 1
	result <- array(0, c(n,n,d))
	for(j in 1:m){
		for(i in 1:j){
			if(i==j) result[ , , l ]<- tcrossprod(Z[ , i ])
				else{
					result[ , , l ] <- tcrossprod(Z[ , j], Z[ , i])
					result[ , , l ] <- result[ , , l ] + t(result[ , , l ])
				}
			l <- l + 1
		}
	}
	result
}

pdDeriv.pdLogChol <-
  function(object, Z, ...){
	Z <- as.matrix(Z)
	m <- ncol(Z)
	d <- m*(m+1)*.5
	n <- nrow(Z)
	l <- 1
	result <- array(0, c(n,n,d))
	for(j in 1:m){
		for(i in 1:j){
			if(i==j) result[ , , l ]<- tcrossprod(Z[, i ])
				else{
					result[ , , l ] <- tcrossprod(Z[ , j], Z[ , i])
					result[ , , l ] <- result[ , , l ] + t(result[ , , l ])
				}
			l <- l + 1
		}
	}
	result
}

pdDeriv.pdDiag <-
  function(object, Z, ...){
	Z <- as.matrix(Z)
	n <- ncol(Z)
	m <- nrow(Z)
	result <- array(0, c(m,m,n))
# print(result)
	for(i in 1:n){
		result[ , , i ] <- tcrossprod(Z[ , i ])
	}
	result
}

pdDeriv.pdCompSymm <-
  function(object, Z, ...){
	Z <- as.matrix(Z)
	n <- nrow(Z)
	dp <- tcrossprod(Z)
	offdp <- tcrossprod(matrix(apply(Z, 1, sum))) - dp
	result <- array(c(dp, offdp), c(n, n, 2))
	result
}

pdDeriv.pdBlocked <-
  function(object, Z, ...){
	Z <- as.matrix(Z)
	n <- nrow(Z)
	varderiv <- lapply(object, function(object, Z){ #per Block
							block.names <- dimnames(as.matrix(object))[[ 1 ]]
							newZ<- Z[ , block.names]
							deriv <- pdDeriv(object, newZ)
							deriv
					}, Z)
	result <- array( unlist(varderiv), c(n, n, length(varderiv)) )
	result
}

pdDeriv.reStruct <-
 function(object, meanEff, ...){
	lapply(meanEff, function(x, object){
			Z <- as.data.frame(x[ , attr(x, "Zcols")$cols])
			grps <- data.frame(x[ , attr(x, "gcols")])
			starts <- attr(x, "Zcols")$starts
			ends <- attr(x, "Zcols")$ends
			Q <- attr(x, "Q")
			result <- vector(mode = "list", length = Q)
			names(result) <- names(object)
			result[[ Q ]] <- pdDeriv(object[[ Q ]], as.data.frame(Z[ , starts[ Q ]:ends[ Q ] ]))
			if(Q > 1){
				for( i in (Q-1):1){
					grps.factor <- factor(grps[ , i ], levels = unique(grps[ , i ]))
					newZ <- split(as.data.frame(Z[ , starts[ i ]:ends[ i ] ]), grps.factor)
					result[[ i ]] <- lapply(newZ, function(x, st){
									res <- pdDeriv(st, x)
									res
									}, st = object[[ i ]])
				}
			}
			result
		}, object)
}

pdDeriv.lmeStruct <-
 function(object, meanEff, ...){
	result <- list()
	result[[ "reStruct" ]] <- pdDeriv(object$reStruct, meanEff)
	result[[ "corStruct" ]] <- NULL
	result[[ "varStruct" ]] <- NULL
	result
}

pdDeriv.respStruct <-
 function(object, second = F, ...){
	if(second) return(object$varsecDer)
		else return(object$varDeriv)
}

pdDeriv.lmeObject <-
 function(object, second = F, ...){
	object <- respStruct(object)
	if(second) return(object$varsecDer)
		else return(object$varDeriv)
}
