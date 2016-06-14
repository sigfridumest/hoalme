## file hoalme/R/conslmeObjects.R, v 0.0-1 2009/11/18
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


##conslmeObjects methods

conslmeObjects <-
 function(object, scale = 2, precision = 10, control = list(returnObject = TRUE))
	UseMethod("conslmeObjects")

conslmeObjects.lmeObject <-
  function(object, scale = 2, precision = 10, control = list(returnObject = TRUE)){
	hypotheses <- NULL
	controlvals <- lmeControl()
	controlvals[names(control)] <- control

	##Constrained optimization on the fixed effects parameters
	npar.fixed <- ncol(object$evalEnv$X)
	if(npar.fixed > 1){
		if(!is.list(scale) && length(scale) == 1) 	results.fixed <- consFixedEffects(object, scale, precision, controlvals)
			else results.fixed <- consFixedEffects(object, scale$fixed, precision, controlvals)
	}
	else{
		results.fixed <- NULL
		npar.fixed <- 0
	}

	##
	##Getting the estimates as initial values for optimization and setting the control list in the "lme' environment
	##
	par <- coef(object)
	nampar <- names(par)
	npar <- length(nampar)
 	if(!is.even(precision)) precision <- precision + 1
	object$evalEnv$controlvals[names(control)] <- control
	
	##The constrained optimization function
	consOpt <- constrainedOptim(object$evalEnv)
	
###### HAS TO BE CHANGED TO HANDLE MORE THAN ONLY RANDOM EFFECTS STRUCTURES #################################
	##Extracting parameter information and constructing a set of plausible hypotheses
	##
	hyp.nams <- names(Coef(object$modelStruct)) #getting more readable names
	
	##Need to preserve this order: reStruct first
	reSt <- object$modelStruct$reStruct 
	if(!is.null(hypotheses))
		hyplist <- hypotheses
		else{
			if(!is.list(scale) && length(scale) == 1) hyplist <- hypotheses(reSt, object, scale = scale, precision = precision)
				else hyplist <- hypotheses(reSt, object, scale = scale$reStruct, precision = precision) #generating re covariance hypotheses
		}
	hypotheses <- matrix(unlist(hyplist), ncol = precision, byrow = T)
	## hypotheses come in certain order... reordering
	ord.hyp <- apply(hypotheses, 1, sort.int, index.return = T)
	ix <- matrix(unlist(lapply(ord.hyp, function(x) x[[ 2 ]])), ncol = precision, byrow = T)
	hypotheses <- matrix(unlist(lapply(ord.hyp, function(x) x[[ 1 ]])), ncol = precision, byrow = T)	
	if(npar == 1)
		index <- c(ix)
		else
			index <- c( t(ix + precision*c( 0,1:(npar-1) )) )

	##
	##Constrained optimization
	##
	results <- vector(length = precision*npar, mode = "list")
	l <- 1
 	
 	
 	for(i in seq_along(hyplist)){ #level
	
	if( inherits(reSt[[ i ]], "pdBlocked") ){
	 	m <- length(reSt[[ i ]])
		par.hyp <- unlist(hyplist[[ i ]], F)
	}
	else par.hyp <- hyplist[[ i ]]
	
	for(j in seq_along(par.hyp)) { #parameter
			hyp <- par.hyp[[ j ]]
			nam <- attr(hyp, "name")
			pos <- match(nam, nampar)
			par.n.hyp <- rep(0, npar)
			
print(nam)
# print(pos)
			
			##if optimization is required or only recalculating at the hypothesis
			if(npar == 1){
				for(k in seq_along(hyp)){ #value
					only.free.par <- hyp[ k ]
					results[[  l  ]] <- consOpt(only.free.par, hyp = NULL, par.n.hyp = NULL, pos = NULL, npar, nam)
					l <- l + 1
				}
			}
			else{
				only.free.par <- par[ -pos] #only free parameters
				for(k in seq_along(hyp)){#hyp is now vector
					results[[  l  ]] <- consOpt(only.free.par, hyp[ k ], par.n.hyp, pos, npar, nam)
					only.free.par <- results[[  l  ]]$unconsPar[ -pos]
					l <- l + 1
				}
			}
		}
	}
	results <- results[index]

	#saving hypothesis name in each result
	for(i in seq_along(hyp.nams)){
		new.index <- (precision*(i -1) + 1):(precision*i)
		for(j in new.index){
			attr(results[[ j ]], "hypothesis") <- hyp.nams[ i ]
		}
	}

	#### Finished constrained optimization
	
	
	
	##Merging all results
	if(npar.fixed > 0)
		results <- c(results.fixed, results)
	
	##Extracting hypotheses, log-Likelihoods and parameter values in 'natural model parameterization'

	logliks <- matrix(unlist(lapply(results, function(x) x$logLik)), ncol = precision, byrow = T )
	logliks <- cbind(object$logLik, logliks)
	
	##1. Random effects covariance parameter
	if(npar.fixed >0){
		hypotheses.fixed <- attr(results.fixed, "hypotheses")
	
		hypotheses <- rbind( cbind(object$coefficients$fixed, hypotheses.fixed),
						cbind(object$unconsPar, hypotheses) )
	}
	else
		hypotheses <- cbind(object$unconsPar, hypotheses)
	
	attr(hypotheses, "npar") <- attr(logliks, "npar") <- nrow(hypotheses)
	attr(hypotheses, "nhyp") <- attr(logliks, "nhyp") <- nhyp <- precision
	
	if(npar.fixed > 0)
		npar <- nrow(hypotheses) - npar.fixed
		else npar <- nrow(hypotheses)
		
	nat.hypotheses <- matrix(0, ncol = (nhyp + 1), nrow = npar)

	nampar <- names(Coef(results[[ 1 ]]))
	
	if(npar.fixed > 0)
		parameters <- array(unlist(lapply(results[-(1:(npar.fixed*precision))], Coef, F)), c(npar, nhyp, npar),
										dimnames = list(nampar, paste("hypothesis.", 1:nhyp , sep = ""), nampar))
		else
		parameters <- array(unlist(lapply(results, Coef, F)), c(npar, nhyp, npar),
										dimnames = list(nampar, paste("hypothesis.", 1:nhyp , sep = ""), nampar))
			
	for( i in 1:npar){
		nat.hypotheses[ i, -1] <- parameters[i , , i]
	}
	nat.hypotheses[ , 1] <- Coef(object, F)
	nat.hyp.nams <- dimnames(parameters)[ 1:2 ]

	nat.hyp.nams[[ 2 ]] <- c("mle", nat.hyp.nams[[ 2 ]])
	dimnames(nat.hypotheses) <- nat.hyp.nams
	dimnames(hypotheses)[[ 2 ]] <- nat.hyp.nams[[ 2 ]]
	
	if(npar.fixed > 1 )
		nat.hypotheses <- rbind(hypotheses[1:nrow(hypotheses.fixed), ] , nat.hypotheses)
	npar <- nrow(hypotheses)


	## wrapping up
	result <- list(logProfLiks = logliks, parameters = parameters, hypotheses = hypotheses, nat.hypotheses = nat.hypotheses,
						objects = results, npar = npar, nhyp = nhyp, npar.fixed = npar.fixed)
	
	#making the 'lmeObject' object available anytime
	attr(result, "lmeObject") <-  object
	attr(result, "npar.fixed") <- npar.fixed
	class(result) <- "conslmeObjects"
	object$conslmeObjects <- result
	result
}

plot.conslmeObjects <-
 function(x, uncons = TRUE, ...){
	if(uncons) hypotheses <- x$hypotheses
		else hypotheses <- x$nat.hypotheses
	npar <- x$npar
	
 	if(npar >= 3) par(mfrow = c(3, 1))
 		else par(mfrow = c(npar, 1))
	
	nams <- dimnames(hypotheses)[[ 1 ]]
	logliks <- x$logProfLiks
	max.loglik <- logliks[1,1]
	qa <- max.loglik - .5*1.281552^2
	qb <- max.loglik - .5*1.644854^2
	qc <- max.loglik - .5*2.326348^2

	if(npar > 3) pages <- floor(npar/3)
		else pages <- 1
	for(i in 1:npar){
		if( any(i == ((1:pages)*3 + 1) ) )
			readline("Press ENTER to continue ...")
 		s.hyp <- order(hypotheses[i, ])
 print(s.hyp)
  		plot(x = hypotheses[ i , s.hyp ], y = logliks[ i , s.hyp ], xlab = nams[ i ], ylab = "(log) Profile Likelihood", type = "l", ... )
#   		plot(x = hypotheses[ i ,  ], y = logliks[ i ,  ], xlab = nams[ i ], ylab = "(log) Profile Likelihood", type = "l", ... )
		abline(coef = c(qa, 0), col = "yellow"); abline(coef = c(qb, 0), col = "orange"); abline(coef = c(qc, 0), col = "red")
	}
}

print.conslmeObjects <-
  function(x, ...){
	cat("** Constrained Maximum Likelihood Estimation **\n")
	cat("                  for the\n")
	cat("       ** Linear Mixed Effects Model **\n\n")
	cat("Profile Likelihood computed for the following parameters:\n\n")
	cat(dimnames(hypotheses(x))[[1]], "\n")
	cat("(", x$nhyp, " hypotheses per parameter).\n", sep = "")
}

constrainedOptim <- function(lmeObjectEnv){
  local({
  ##
  ## constrainedOptim creates a function that has attached the evaluation environment of the lmeObject used to fit the LME model.
  ##
	  function(par, hypothesis, par.n.hyp, pos, npar, nam){
 		
 		if(controlvals$msVerbose)
 			cat("\nHypothesis:", nam, "\n", hypothesis, "\n")
				
 		if(npar== 1){
 			coef(lmeSt) <- par
 			attr(lmeSt, "lmeFit") <- nlme:::MEestimate(lmeSt, grps)
 			lmeSt <- update(lmeSt, dataMix)
 		}
		else {
			##Optimization
			## getting the linear mixed effects fit object,
			## possibly iterating for variance functions
			##
			numIter <- 0
			repeat {
			oldPars <- coef(lmeSt)
#  			oldPars <- par.n.hyp
			oldPars[ pos] <- hypothesis
			optRes <- if (controlvals$opt == "nlminb") {
				nlminb(par, function(par, hypothesis, par.n.hyp, pos, lmeSt){
						par.n.hyp[ pos ] <- hypothesis
						par.n.hyp[-pos ] <- par
						-logLik(lmeSt, par.n.hyp)
						}, gradient = NULL, hessian = NULL,
					hypothesis, par.n.hyp, pos, lmeSt,
					control = list(iter.max = controlvals$msMaxIter,
					eval.max = controlvals$msMaxEval,
					trace = controlvals$msVerbose)
				)
				} else {
					optim(par, function(par, hypothesis, par.n.hyp, pos){
						par.n.hyp[ pos ] <- hypothesis
						par.n.hyp[-pos ] <- par
						-logLik(lmeSt, par.n.hyp)
						}, , hypothesis, par.n.hyp, pos,
					control = list(trace = controlvals$msVerbose,
					maxit = controlvals$msMaxIter,
					reltol = if(numIter == 0) controlvals$msTol
					else 100*.Machine$double.eps),
					method = controlvals$optimMethod)
				}
				numIter0 <- NULL
				
				par.n.hyp[-pos ] <- optRes$par
				par.n.hyp[ pos ] <- hypothesis
				coef(lmeSt) <- par.n.hyp
				attr(lmeSt, "lmeFit") <- nlme:::MEestimate(lmeSt, grps)
				## checking if any updating is needed
				if (!needUpdate(lmeSt)) {
					if (optRes$convergence) {
					msg <- paste(controlvals$opt, " problem, convergence error code = ",
							optRes$convergence, "\n  message = ", optRes$message,
							sep='')
					if(!controlvals$returnObject)
						stop(msg)
					else
						warning(msg)
					}
					break
				}
				
				## updating the fit information
				numIter <- numIter + 1
				lmeSt <- update(lmeSt, dataMix)
				## calculating the convergence criterion
				aConv <- coef(lmeSt)
				conv <- abs((oldPars - aConv)/ifelse(aConv == 0, 1, aConv))
				aConv <- NULL
				for(i in names(lmeSt)) {
					if (any(parMap[,i])) {
					aConv <- c(aConv, max(conv[parMap[,i]]))
					names(aConv)[length(aConv)] <- i
					}
				}
				if (max(aConv) <= controlvals$tolerance) {
					break
				}
				if (numIter > controlvals$maxIter) {
					msg <- paste("Maximum number of iterations",
						"(lmeControl(maxIter)) reached without convergence.")
					if (controlvals$returnObject) {
					warning(msg)
					break
					} else
					stop(msg)
				}
			} ## end{repeat}
		}
			## wrapping up
			lmeFit <- attr(lmeSt, "lmeFit")
			names(lmeFit$beta) <- namBeta <- colnames(X)
			attr(fixDF, "varFixFact") <- varFix <- lmeFit$sigma * lmeFit$varFix
			varFix <- crossprod(varFix)
			dimnames(varFix) <- list(namBeta, namBeta)
			##
			## fitted.values and residuals (in original order)
			##
			Fitted <- fitted(lmeSt, level = 0:Q,
					conLin = if (decomp) oldConLin else attr(lmeSt, "conLin"))[
					revOrder, , drop = FALSE]
			Resid <- y[revOrder] - Fitted
			rownames(Resid) <- rownames(Fitted) <- origOrder
			attr(Resid, "std") <- lmeFit$sigma/(varWeights(lmeSt)[revOrder])
			## putting groups back in original order
			grps <- grps[revOrder, , drop = FALSE]
			
			## making random effects estimates consistently ordered
			#  for(i in names(lmeSt$reStruct)) {
			#    lmeFit$b[[i]] <- lmeFit$b[[i]][unique(as.character(grps[, i])),, drop = F]
			#    NULL
			#  }
			
			###
			###log-likelihood and unconstrained parameters and hessian, before **any** inverse transformation is made
			###
			unconsPar <- coef(lmeSt)
			hess <- list()
			hess$mean <- logLik(lmeSt, unconsPar)
			attributes(hess$mean) <- NULL
			hess$Hessian <- hessian(function(pars) logLik(lmeSt, pars), unconsPar)
			hess$gradient <- grad(function(pars) logLik(lmeSt, pars), unconsPar)
			dimnames(hess$Hessian) <- list(names(unconsPar), names(unconsPar))

			unsolved.reStruct <- lmeSt$reStruct

			## inverting back reStruct
			lmeSt$reStruct <- solve(lmeSt$reStruct)
			## saving part of dims
			dims <- attr(lmeSt, "conLin")$dims[c("N", "Q", "qvec", "ngrps", "ncol")]
			## getting the approximate var-cov of the parameters
# 			if (controlvals$apVar) {
# 			apVar <- lmeApVar(lmeSt, lmeFit$sigma,
# 					.relStep = controlvals[[".relStep"]],
# 					minAbsPar = controlvals[["minAbsParApVar"]],
# 					natural = controlvals[["natural"]])
# 			} else {
# 			apVar <- "Approximate variance-covariance matrix not available"
# 			}
			
			#Preserving condensed linear model and fit, as they will be used as attributes later in CMle
# 			# getting rid of condensed linear model and fit
# 			attr(lmeSt, "conLin") <- NULL
# 			attr(lmeSt, "lmeFit") <- NULL
			##
			## creating the  lme object
			##
			
			#putting into evaluation enviroment the model matrices from the 'lmeObject' object
			eval.env <- mget(ls(environment( )), environment( ))
			eval.env$X <- X
			eval.env$Z <- Z
			eval.env$y <- y
			
			estOut <- list(modelStruct = lmeSt,
				evalEnv = eval.env, #saving all local variables
				Hessian = hess,
				unconsPar = unconsPar,
				unsolved.reStruct = unsolved.reStruct,
				optRes = optRes,
				dims = dims,
				contrasts = contr,
				coefficients = list(
					fixed = lmeFit$beta,
					random = lmeFit$b),
				varFix = varFix,
				sigma = lmeFit$sigma,
				apVar = apVar,
				logLik = lmeFit$logLik,
				numIter = if (needUpdate(lmeSt)) numIter
				else numIter0,
				groups = grps,
				call = Call,
				terms = Terms,
				method = method,
				fitted = Fitted,
				residuals = Resid,
				fixDF = fixDF,
				na.action = attr(dataMix, "na.action"))
			if (keep.data && !miss.data) estOut$data <- data
			if (inherits(data, "groupedData")) {
			## saving labels and units for plots
			attr(estOut, "units") <- attr(data, "units")
			attr(estOut, "labels") <- attr(data, "labels")
			}
			class(estOut) <- c("conslmeObject","lmeObject")
			invisible(estOut)
		}
	}, lmeObjectEnv)
}

