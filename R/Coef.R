## file hoalme/R/Coef.R, v 0.0-1 2009/11/18
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


##Coef methods

Coef <-
 function(Struct, ...)
	UseMethod("Coef")

Coef.pdSymm <-
  function(Struct, ...){
	par <- coef(Struct, F)
	par
}

Coef.pdLogChol <-
  function(Struct, ...){
	par <- coef(Struct, F)
	par
}

Coef.pdCompSymm <-
  function(Struct, ...){
# 	nam <- Names(Struct)
	par <- coef(Struct, F)
	par[ 1 ] <- par[ 1 ]^2
	par[ 2 ] <- par[2]*par[1]
 	names(par) <- c("var", "cov")
	par
}

Coef.pdIdent <-
  function(Struct, ...){
# 	nam <- Names(Struct)
	par <- coef(Struct, F)^2
	names(par) <- "var"
	par
}

Coef.pdDiag <-
  function(Struct, ...){
	par <- coef(Struct, F)^2
	names(par) <- paste("var(", Names(Struct),")", sep ="")
	par
}

Coef.pdNatural <-
  function(Struct, ...){
	par <- coef(Struct, F)
	par
}

Coef.pdBlocked <-
  function(Struct, ...){
	par <- unlist(lapply(Struct, Coef))
	par
}


Coef.reStruct <-
  function(Struct, ...){
	nams <- names(Struct)
	pars <- list()
	for(i in seq_along(Struct)){
		newnams <- paste(nams[ i ], 1:length(Struct[[ i ]]), sep = "")
		par <- Coef(Struct[[ i ]])
		oldnams <- names(par)
		names(par) <- paste(newnams, oldnams, sep = ".")
		pars[[ i ]] <- par
	}
	unlist(pars)
}



Coef.modelStruct <-
  function(Struct, ...){
	par <- unlist(lapply(Struct, Coef))
	par
}

Coef.lmeStructInt <-
  function(Struct, ...){
	par <- unlist(lapply(Struct, Coef))
	par
}

Coef.lmeObject <-
  function(Struct, scaled = T, ...){
	lmeSt <- Struct$modelStruct
	if(scaled) par <- Coef(lmeSt)
		else{
			sigmasq <- Struct$sigma^2
			par <- Coef(lmeSt)*sigmasq
		}
	par
}


"Coef<-" <-
  function(object, ..., value)
	UseMethod("Coef<-")

"Coef<-.pdSymm" <-
   function(object, ..., value){
# print(value)
	npar <- length(value)
	ndiag.entries <- as.integer((-1 + sqrt(1 + 8*npar))*.5)
	mat <- diag(ndiag.entries)
# print(mat)
	mat[col(mat) >= row(mat)] <- value
	mat[col(mat) < row(mat)] <- mat[col(mat) > row(mat)]
# print(mat)
	coef(object) <- coef(pdConstruct(object, mat))
	object
}

"Coef<-.pdLogChol" <-
   function(object, ..., value){
# print(value)
	npar <- length(value)
	ndiag.entries <- as.integer((-1 + sqrt(1 + 8*npar))*.5)
	mat <- diag(ndiag.entries)
# print(mat)
	mat[col(mat) >= row(mat)] <- value
	mat[col(mat) < row(mat)] <- mat[col(mat) > row(mat)]
# print(mat)
	coef(object) <- coef(pdConstruct(object, mat))
	object
}


"Coef<-.pdCompSymm" <-
   function(object, ..., value){
	npar <- dim(as.matrix(object))[ 1 ]
	mat <- diag(npar)
	diag(mat) <- value[ 1 ]
	mat[col(mat) != row(mat)] <- value[ 2 ]
	coef(object) <- coef(pdConstruct(object, mat))
	object
}

"Coef<-.pdIdent" <-
    function(object, ..., value){
	mat <- diag(ncol(as.matrix(object)))
	diag(mat) <- value
 	coef(object) <- coef(pdConstruct(object, mat))
 	object
 }

"Coef<-.pdDiag" <-
   function(object, ..., value){
	mat <- diag(ncol(as.matrix(object)))
	diag(mat) <- value
	coef(object) <- coef(pdConstruct(object, mat))
	object
}

"Coef<-.pdBlocked" <-
   function(object, ..., value){
	nspar <- unlist(lapply( object, function(x) length(coef(x)) ))
	ends <- cumsum(nspar)
	starts <- 1 + c(0, ends[-length(ends)])
	for (i in seq_along(object)) {
    		Coef(object[[i]]) <- value[(starts[i]):(ends[i])]
	}
	object
}

"Coef<-.reStruct" <-
   function(object, ..., value){
	npar <- attr(object, "plen")
	ends <- cumsum(npar)
	starts <- 1 + c(0, ends[-length(ends)])
	for (i in seq_along(object)) {
    		Coef(object[[i]]) <- value[(starts[i]):(ends[i])]
	}
	object
}

"Coef<-.lmeStruct" <-
   function(object, ..., value){
	value <- as.numeric(value)
	parMap <- attr(object, "pmap")
	for(i in names(object)){
		if (any(parMap[,i])){
			Coef(object[[i]]) <- value[parMap[,i]]
		}
	}
	object
}

"Coef<-.modelStruct" <-
  function(object, ..., value){
	value <- as.numeric(value)
	parMap <- attr(object, "pmap")
	for(i in names(object)){
		if (any(parMap[,i])){
			Coef(object[[i]]) <- value[parMap[,i]]
		}
	}
	object
}

"Coef<-.lmeStructInt"<- 
   function(object, ..., value){
	value <- as.numeric(value)
	parMap <- attr(object, "pmap")
	for(i in names(object)){
		if (any(parMap[,i])){
			Coef(object[[i]]) <- value[parMap[,i]]
		}
	}
	object
}


"Coef<-.lmeObject" <-
   function(object, sigmasq = 1, ..., value){
	lmeSt <- object$modelStruct
	Coef(lmeSt) <- value/sigmasq
	object$modelStruct <- lmeSt
	object
}

