## file hoalme/R/Diag.R, v 0.0-1 2009/11/18
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

##Diag methods


Diag <-
  function(object, ...)
	UseMethod("Diag")
	
Diag.pdDiag <-
  function(object, ...){
 	rep(T, length(object))
}

Diag.pdIdent <-
  function(object, ...){
 	rep(T, length(object))
}

Diag.pdSymm <-
  function(object, ...){
 	len <- length(object)
 	ndiag.entries <- as.integer((-1 + sqrt(1 + 8*len))*.5)
	diag.entries <- cumsum(1:ndiag.entries)
	diag.dummy <- rep(0, len)
	diag.dummy[diag.entries] <- 1
	diag.dummy[-diag.entries] <- 0
	as.logical(diag.dummy)
}

Diag.pdLogChol <-
  function(object, ...){
 	len <- length(object)
 	ndiag.entries <- as.integer((-1 + sqrt(1 + 8*len))*.5)
	diag.entries <- cumsum(1:ndiag.entries)
	diag.dummy <- rep(0, len)
	diag.dummy[diag.entries] <- T
	diag.dummy[-diag.entries] <- F
	diag.dummy
}

Diag.pdCompSymm <-
  function(object, ...){
	c(T, F)
}

Diag.pdBlocked <-
  function(object, ...){
	unlist(lapply(object, Diag))
}

Diag.reStruct <-
  function(object, ...){
 	unlist(lapply(object, Diag))
}

Diag.lmeStructInt <-
  function(object, ...){
 	unlist(lapply(object, Diag))
}

Diag.lmeStruct <-
  function(object, ...){
 	unlist(lapply(object, Diag))
}

