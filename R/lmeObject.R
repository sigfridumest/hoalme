## file hoalme/R/lmeObject.R, v 0.0-1 2009/11/18
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


##lmeObject methods

lmeObject <-
  function(fixed, ...)
	UseMethod("lmeObject")

lmeObject.formula <-
  function(fixed,
	   data = sys.frame(sys.parent()),
	   random = pdSymm( eval( as.call( fixed[ -2 ] ) ) ),
	   correlation = NULL,
	   weights = NULL,
	   subset,
	   method = c("ML", "REML"), #sobra
	   na.action = na.fail,
	   control = list(natural = F),
       contrasts = NULL,
       keep.data = TRUE,
	   start = NULL, ...)
{
  Call <- match.call()
  miss.data <- missing(data) || !is.data.frame(data)

  ## control parameters
  controlvals <- lmeControl()
  if (!missing(control)) {
    if(!is.null(control$nlmStepMax) && control$nlmStepMax < 0) {
      warning("Negative control$nlmStepMax - using default value")
      control$nlmStepMax <- NULL
    }
    controlvals[names(control)] <- control
  }

  ##
  ## checking arguments
  ##
  if (!inherits(fixed, "formula") || length(fixed) != 3) {
    stop("\nFixed-effects model must be a formula of the form \"resp ~ pred\"")
  }
  if(is.null(Call$random)){
	Call$random <- call("pdSymm", eval( as.call( fixed[ -2 ] )))
	}
##coercing method to be "ML"
#   method <- match.arg(method)
  method <- "ML"
  REML <- method == "REML"
  reSt <- reStruct(random, REML = REML, data = NULL)
  groups <- getGroupsFormula(reSt)
  if (is.null(groups)) {
    if (inherits(data, "groupedData")) {
      groups <- getGroupsFormula(data)
      namGrp <- rev(names(getGroupsFormula(data, asList = TRUE)))
      Q <- length(namGrp)
      if (length(reSt) != Q) { # may need to repeat reSt
	if (length(reSt) != 1) {
	  stop("Incompatible lengths for \"random\" and grouping factors")
	}
        randL <- vector("list", Q)
        names(randL) <- rev(namGrp)
        for(i in 1:Q) randL[[i]] <- random
        randL <- as.list(randL)
	reSt <- reStruct(randL, REML = REML, data = NULL)
      } else {
	names(reSt) <- namGrp
      }
    } else {
      ## will assume single group
      groups <- ~ 1
      names(reSt) <- "1"
    }
  }

  ## check if corStruct is present and assign groups to its formula,
  ## if necessary
  if (!is.null(correlation)) {
    if(!is.null(corGrpsForm <- getGroupsFormula(correlation, asList = TRUE))) {
      corGrpsForm <- unlist(lapply(corGrpsForm,
                                   function(el) deparse(el[[2]])))
      corQ <- length(corGrpsForm)
      lmeGrpsForm <- unlist(lapply(splitFormula(groups),
                        function(el) deparse(el[[2]])))
      lmeQ <- length(lmeGrpsForm)
      if (corQ <= lmeQ) {
        if (any(corGrpsForm != lmeGrpsForm[1:corQ])) {
          stop(paste("Incompatible formulas for groups in \"random\"",
                     "and \"correlation\""))
        }
        if (corQ < lmeQ) {
          warning(paste("Cannot use smaller level of grouping for",
                        "\"correlation\" than for \"random\". Replacing",
                        "the former with the latter."))
          attr(correlation, "formula") <-
            eval(parse(text = paste("~",
                    nlme:::c_deparse(getCovariateFormula(formula(correlation))[[2]]),
                         "|", deparse(groups[[2]]))))
        }
      } else {
        if (any(lmeGrpsForm != corGrpsForm[1:lmeQ])) {
          stop(paste("Incompatible formulas for groups in \"random\"",
                     "and \"correlation\""))
        }
      }
    } else {
      ## using the same grouping as in random
      attr(correlation, "formula") <-
        eval(parse(text = paste("~",
		     c_deparse(getCovariateFormula(formula(correlation))[[2]]),
		     "|", deparse(groups[[2]]))))
      corQ <- lmeQ <- 1
    }
    } else {
    corQ <- lmeQ <- 1
  }
  ## create an lme structure containing the random effects model and plug-ins
  lmeSt <- lmeStruct(reStruct = reSt, corStruct = correlation,
		     varStruct = varFunc(weights))

  ## extract a data frame with enough information to evaluate
  ## fixed, groups, reStruct, corStruct, and varStruct
  mfArgs <- list(formula = asOneFormula(formula(lmeSt), fixed, groups),
		 data = data, na.action = na.action)
  if (!missing(subset)) {
    mfArgs[["subset"]] <- asOneSidedFormula(Call[["subset"]])[[2]]
  }
  mfArgs$drop.unused.levels <- TRUE
  dataMix <- do.call("model.frame", mfArgs)
  origOrder <- row.names(dataMix)	# preserve the original order
  for(i in names(contrasts))            # handle contrasts statement
      contrasts(dataMix[[i]]) = contrasts[[i]]
  ## sort the model.frame by groups and get the matrices and parameters
  ## used in the estimation procedures
  grps <- getGroups(dataMix, groups)
  ## ordering data by groups
  if (inherits(grps, "factor")) {	# single level
    ord <- order(grps)	#"order" treats a single named argument peculiarly
    grps <- data.frame(grps)
    row.names(grps) <- origOrder
    names(grps) <- as.character(deparse((groups[[2]])))
  } else {
    ord <- do.call("order", grps)
    ## making group levels unique
    for(i in 2:ncol(grps)) {
      grps[, i] <-
        as.factor(paste(as.character(grps[, i-1]), as.character(grps[,i]),
                        sep = "/"))
      NULL
    }
  }
  if (corQ > lmeQ) {
    ## may have to reorder by the correlation groups
    ord <- do.call("order", getGroups(dataMix,
                                 getGroupsFormula(correlation)))
  }
  grps <- grps[ord, , drop = FALSE]
  dataMix <- dataMix[ord, ,drop = FALSE]
  revOrder <- match(origOrder, row.names(dataMix)) # putting in orig. order

  ## obtaining basic model matrices
  N <- nrow(grps)
  Z <- model.matrix(reSt, dataMix)
  ncols <- attr(Z, "ncols")
  Names(lmeSt$reStruct) <- attr(Z, "nams")
  ## keeping the contrasts for later use in predict
  contr <- attr(Z, "contr")
  X <- model.frame(fixed, dataMix)
  Terms <- attr(X, "terms")
  auxContr <- lapply(X, function(el)
		     if (inherits(el, "factor") &&
                         length(levels(el)) > 1) contrasts(el))
  contr <- c(contr, auxContr[is.na(match(names(auxContr), names(contr)))])
  contr <- contr[!unlist(lapply(contr, is.null))]
  X <- model.matrix(fixed, data = X)
  y <- eval(fixed[[2]], dataMix)
  ncols <- c(ncols, dim(X)[2], 1)
  Q <- ncol(grps)
  ## creating the condensed linear model
  attr(lmeSt, "conLin") <-
    list(Xy = array(c(Z, X, y), c(N, sum(ncols)),
	     list(row.names(dataMix), c(colnames(Z), colnames(X),
					deparse(fixed[[2]])))),
	 dims = nlme:::MEdims(grps, ncols), logLik = 0)
 
  ## checking if enough observations per group to estimate ranef
  tmpDims <- attr(lmeSt, "conLin")$dims
  if (max(tmpDims$ZXlen[[1]]) < tmpDims$qvec[1]) {
    warning(paste("Fewer observations than random effects in all level",
                  Q,"groups"))
  }
  ## degrees of freedom for testing fixed effects
  fixDF <- nlme:::getFixDF(X, grps, attr(lmeSt, "conLin")$dims$ngrps,
                    terms = Terms)
  ## initialization
  lmeSt <- Initialize(lmeSt, dataMix, grps, control = controlvals)
  parMap <- attr(lmeSt, "pmap")
  ## Checking possibility of single decomposition
  if (length(lmeSt) == 1)  {	# reStruct only, can do one decomposition
    ## need to save conLin for calculating fitted values and residuals
    oldConLin <- attr(lmeSt, "conLin")
    decomp <- TRUE
    attr(lmeSt, "conLin") <- nlme:::MEdecomp(attr(lmeSt, "conLin"))
  } else decomp <- FALSE

# print(attr(lmeSt, "conLin"))
  ##
  ## getting the linear mixed effects fit object,
  ## possibly iterating for variance functions
  ##
# print(class(lmeSt))
  numIter <- 0
  repeat {
    oldPars <- coef(lmeSt)
    optRes <- if (controlvals$opt == "nlminb") {
        nlminb(c(coef(lmeSt)),
               function(lmePars) -logLik(lmeSt, lmePars),
               control = list(iter.max = controlvals$msMaxIter,
               eval.max = controlvals$msMaxEval,
               trace = controlvals$msVerbose))
    } else {
        optim(c(coef(lmeSt)),
              function(lmePars) -logLik(lmeSt, lmePars),
              control = list(trace = controlvals$msVerbose,
              maxit = controlvals$msMaxIter,
              reltol = if(numIter == 0) controlvals$msTol
              else 100*.Machine$double.eps),
              method = controlvals$optimMethod)
    }
    numIter0 <- NULL
    coef(lmeSt) <- optRes$par

#    cat("\n acabado de salir de la optimizacion \n")
#    print(logLik(lmeSt, coef(lmeSt)))

   #optimization over
  
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

###log-likelihood and unconstrained parameters and hessian, before **any** inverse transformation is made
	unconsPar <- coef(lmeSt)
	hess <- list()
	hess$mean <- logLik(lmeSt, unconsPar)
	attributes(hess$mean) <- NULL
	hess$Hessian <- hessian(function(pars) logLik(lmeSt, pars), unconsPar)
	hess$gradient <- grad(function(pars) logLik(lmeSt, pars), unconsPar)
	dimnames(hess$Hessian) <- list(names(unconsPar), names(unconsPar))

	unsolved.reStruct <- lmeSt$reStruct

  ## making random effects estimates consistently ordered
#
#
#
  ## inverting back reStruct
  lmeSt$reStruct <- solve(lmeSt$reStruct)

  ## saving part of dims
  dims <- attr(lmeSt, "conLin")$dims[c("N", "Q", "qvec", "ngrps", "ncol")]
  ## getting the approximate var-cov of the parameters
  if (controlvals$apVar) {
    apVar <- nlme:::lmeApVar(lmeSt, lmeFit$sigma,
		      .relStep = controlvals[[".relStep"]],
                      minAbsPar = controlvals[["minAbsParApVar"]],
		      natural = controlvals[["natural"]])
  } else {
    apVar <- "Approximate variance-covariance matrix not available"
  }
   ##Preserving condensed linear model and fit, as they will be used as attributes later in CMle
  # getting rid of condensed linear model and fit
  #attr(lmeSt, "conLin") <- NULL
  #attr(lmeSt, "lmeFit") <- NULL
  ##
  ## creating the  lme object
  ##

  estOut <- extend(Object(), "lmeObject", modelStruct = lmeSt, respStruct = NULL, conslmeObjs = NULL,
		 evalEnv = mget(ls(environment( )), environment( )), #saving all local variables
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
                 na.action = attr(dataMix, "na.action"),
		parentframe = parent.frame())
  if (keep.data && !miss.data) estOut$data <- data
	else
		if(is.environment(data)) estOut$data <- data
			else estOut$data <- deparse(substitute(data))
  if (inherits(data, "groupedData")) {
    ## saving labels and units for plots
    attr(estOut, "units") <- attr(data, "units")
    attr(estOut, "labels") <- attr(data, "labels")
  }

#   invisible(estOut)
     estOut

}

print.lmeObject <-
 function(x, ...){
	class(x) <- c("lme", class(x))
# 	print(class(x))
   	print(x, ...)
	dd <- x$dims
	cat("Number of Independent Vectors:", dd$ngrps[dd$Q], "\n")
	cat("\nFitting information:\n")
	optRes <- x$optRes
	cat("\n Convergence type: ", optRes$convergence, " (", optRes$message, ")", sep = "")
	cat("\n Iterations:", optRes$iterations)
	cat("\n Hessian: \n")
	print(Hessian(x))
}

summary.lmeObject <-
 function(object, adjustSigma = TRUE, verbose = FALSE, ...){
  class(object) <- c("lme", class(object))
  ##  variance-covariance estimates for fixed effects
  fixed <- fixef(object)
  stdFixed <- sqrt(diag(as.matrix(object$varFix)))
  object$corFixed <- array(t(object$varFix/stdFixed)/stdFixed,
                           dim(object$varFix), list(names(fixed),names(fixed)))
  if (object$method == "ML" && adjustSigma == TRUE) {
    stdFixed <-
      sqrt(object$dims$N/(object$dims$N - length(stdFixed))) * stdFixed
  }
  ## fixed effects coefficients, std. deviations and t-ratios
  ##
  tTable <- data.frame(fixed, stdFixed, object$fixDF[["X"]],
                       fixed/stdFixed, fixed)
  dimnames(tTable)<-
    list(names(fixed),c("Value", "Std.Error", "DF", "t-value", "p-value"))
  tTable[, "p-value"] <- 2 * pt(-abs(tTable[,"t-value"]), tTable[,"DF"])
  object$tTable <- as.matrix(tTable)
  ##
  ## residuals
  ##
  resd <- resid(object, type = "pearson")
  if (length(resd) > 5) {
    resd <- quantile(resd, na.rm = TRUE) # might have NAs from na.exclude
    names(resd) <- c("Min","Q1","Med","Q3","Max")
  }
  object$residuals <- resd
  ##
  ## generating the final object
  ##
  aux <- logLik(object)
  object$BIC <- BIC(aux)
  object$AIC <- AIC(aux)
  attr(object, "oClass") <- class(object)
  attr(object, "verbose") <- verbose
  class(object) <- c("summary.lme", class(object))
  object
  print(object)
  dd <- object$dims
  cat("\nNumber of Independent Vectors:", dd$ngrps[dd$Q], "\n")
}

coef.lmeObject <-
  function(object, ...){
	object$unconsPar
}

unconsVar <-
 function(lmeObj){
	##Extracting parameters and hessian
	par <- coef(lmeObj)
	apvar <- solve(-lmeObj$Hessian$Hessian)
	dimnames(apvar) <- list(names(par), names(par))
	attr(apvar, "par") <- par	
 	apvar
}

model.matrix.lmeObject <-
  function(object, ...){
	Z <- object$evalEnv$Z
	y <- object$evalEnv$y
	X <- object$evalEnv$X
	parFix <- matrix(object$coefficients$fixed, ncol = 1)
	r <- y - X%*%parFix
	
	dims <- object$dims
	ncols <- dims$ncol #columns: Z 1st level ... highest level, matrix X, response
	N <- dims$N #total number of observations
	Q <- dims$Q #number of groups

	G <- rev(object$groups) #factors for grouping; come froem highest to lowest: reversing
# 	GZXy <- data.frame(G, Z, X, y)
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

lmeObject.lme <-
  function(fixed, ...){
	call <- fixed$call
	call$method <- "ML"
	args <- as.list(call)[-1]
	do.call("lmeObject", args)
}

lmeObjectHessian <-
  function(object, conLin = attr(object$modelStruct, "conLin"))
{
  ## calculate approximate variance-covariance matrix of all parameters
  ## except the fixed effects. 
	lmeSt <- object$modelStruct
	pars <- c(object$unconsPar, residual = log(object$sigma^2))
	dims <- conLin$dims
	settings <- attr(lmeSt, "settings")
	N <- dims$N - settings[1] * dims$ncol[dims$Q + 1]
	settings[2:3] <- c(1, 0)			# asDelta = TRUE and no grad/Hess
 	conLin[["logLik"]] <- 0
	p <- length(pars)

	laf <- function(pars, lmeSt, dims, settings, conLin, N, p){
		sig <- sqrt(exp(pars[p]))
		pars <- pars[-p]
		lmeSt$reStruct <- solve(lmeSt$reStruct)
		coef(lmeSt) <- pars
	      	val <- .C(nlme:::mixed_loglik,
			as.double(conLin$Xy),
			as.integer(unlist(dims)),
			as.double(unlist(pdFactor(lmeSt$reStruct))),
			as.integer(settings),
			logLik = double(1),
			lRSS = double(1))[c("logLik", "lRSS")]
			aux <- (exp(val[["lRSS"]])/sig)^2
			conLin[["logLik"]] + val[["logLik"]] + (N * log(aux) - aux)/2
	}

	hess <- hessian(laf, pars, , list(eps=1e-4, d=1e-3,
					zero.tol=100*.Machine$double.eps,  r=4, show.details = F), lmeSt, dims, settings, conLin, N, p)
	lambda.ofAlpha <-
		function(pars, lmeSt, p, i){
			pars[ c(Diag(lmeSt), T) ] <- exp(pars[ c(Diag(lmeSt), T) ])
			sigmasq <- pars[p]
			pars <- pars[-p]
			Coef(lmeSt) <- pars/sigmasq
			lmeSt$reStruct <- solve(lmeSt$reStruct)
			pars <- c(coef(lmeSt), residual = log(sigmasq))
			pars[ i ]
		}
	derivs <- diag(p)
	nat.Pars <- Pars <- c(Coef(object, F), residual = object$sigma^2)
	Pars[ c(Diag(lmeSt), T) ] <- log(Pars[ c(Diag(lmeSt), T) ])
	for(i in 1:p){
		derivs[ i , ] <- grad(lambda.ofAlpha, Pars, , , , object$modelStruct, p, i)
	}
	nat.Pars <- 1/nat.Pars
	nat.Pars[ !c(Diag(lmeSt), T) ] <- 1
	nat.Pars <- diag(nat.Pars)
	nat.Pars%*%t(derivs)%*%hess%*%derivs%*%nat.Pars
}

