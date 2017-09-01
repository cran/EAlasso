#' @title Compute lasso estimator
#'
#' @description Computes lasso, group lasso, scaled lasso, or scaled group lasso solution.
#' The outputs are coefficient-estimate and subgradient. If \code{type = "slasso"}
#' or \code{type = "sgrlasso"}, the output will include the sigma-estimate.
#'
#' @param X predictor matrix.
#' @param Y response vector.
#' @param type type of penalty. Must be specified to be one of the following:
#'  \code{"lasso", "grlasso", "slasso"} or \code{"sgrlasso"}.
#' @param lbd penalty term of lasso. By letting this argument be \code{"cv.1se"} or
#' \code{"cv.min"}, users can have the cross-validated lambda that gives either minimum
#' squared error or that is within 1 std error bound.
#' @param weights weight vector with length equal to the number of groups. Default is
#' \code{rep(1, max(group))}.
#' @param group \code{p} x \code{1} vector of consecutive integers describing the group structure.
#' The number of groups should be the same as max(group). Default is \code{group = 1:p}
#' , where \code{p} is number of covariates.
#' @param verbose logical. Only available for \code{type = "slasso"} or \code{type = "sgrlasso"}.
#' @param ... auxiliary arguments for \code{lbd = "cv.min", lbd = "cv.1se"}.
#' See \code{\link{cv.lasso}} for details.
#' @details
#' Computes lasso, group lasso, scaled lasso, or scaled group lasso solution.
#' Users can specify the value of lbd or choose to run cross-validation to get
#' optimal lambda in term of mean squared error.
#'
#' @return \item{B0}{coefficient estimator.}
#' @return \item{S0}{subgradient.}
#' @return \item{lbd, weights, group}{same as input arguments.}
#' @examples
#' set.seed(123)
#' n <- 50
#' p <- 10
#' X <- matrix(rnorm(n*p), n)
#' Y <- X %*% c(1, 1, rep(0, p-2)) + rnorm(n)
#' #
#' # lasso
#' #
#' Lasso.MHLS(X = X, Y = Y, type = "lasso", lbd = .5)
#' #
#' # group lasso
#' #
#' Lasso.MHLS(X = X, Y = Y, type = "grlasso", lbd = .5, weights = rep(1,2),
#' group = rep(1:2, each=5))
#' @export
Lasso.MHLS <- function(X, Y, type, lbd,
  group=1:ncol(X), weights=rep(1,max(group)), verbose = FALSE, ...)
{
  n <- nrow(X)
  p <- ncol(X)
  X <- as.matrix(X)
  Y <- matrix(Y, , 1)

  if (!type %in% c("lasso", "grlasso", "slasso", "sgrlasso")) {
    stop("type has to be either lasso, grlasso, slasso or sgrlasso.")
  }

  if (!all(group==1:p) && (!type %in% c("grlasso", "sgrlasso"))) {
    stop("Choose type to be either grlasso or sgrlasso if group-structure exists.")
  }

  # if (all(group==1:p) && (!type %in% c("lasso", "slasso"))) {
  #   stop("Choose type to be either lasso or slasso if group-structure does not exist.")
  # }

  if (!lbd  %in% c("cv.1se", "cv.min")) {
    if (!is.numeric(lbd) || lbd <= 0) {stop("invalid lbd input.")}
  }

  if (verbose) {cat("# Cross-validation \n")}
  if (lbd %in% c("cv.1se", "cv.min")) {
    lbdTEMP <- cv.lasso(X = X, Y = Y, group = group, weights = weights,
                    type = type, verbose=verbose, ...)
    if (lbd == "cv.1se") {lbd <- lbdTEMP$lbd.1se} else {
      lbd <- lbdTEMP$lbd.min
    }
  }

  # if (missing(lbd)) {
  #   if (type %in% c("lasso", "grlasso")) {
  #     lbd <- .37
  #   } else {
  #     lbd <- .5
  #   }
  # }

  # if (lbd <= 0) {
  #   stop("lbd has to be positive.")
  # }
  if (length(group) != p) {
    stop("length(group) has to be the same with ncol(X)")
  }
  if (length(weights) != length(unique(group))) {
    stop("length(weights) has to be the same with the number of groups")
  }
  if (any(weights <= 0)) {
    stop("weights should be positive.")
  }
  if (any(!1:max(group) %in% group)) {
    stop("group index has to be a consecutive integer starting from 1.")
  }

  if (length(Y) != n) {
    stop("dimension of X and Y are not conformable.")
  }

  IndWeights <- rep(weights,table(group))
  # scale X with weights
  Xtilde   <- scale(X,FALSE,scale=IndWeights)

  # slassoLoss <- function(X,Y,beta,sig,lbd) {
  #   n <- nrow(X)
  #   crossprod(Y-X%*%beta) / 2 / n / sig + sig / 2 + lbd * sum(abs(beta))
  # }

  if (type %in% c("lasso", "grlasso")) {
    # compute group lasso estimator B0 and S0
    B0 <- coef(gglasso(Xtilde, Y, pf = rep(1,max(group)), group = group,
                       loss="ls", intercept=F, lambda=lbd))[-1] / IndWeights
    S0 <- (t(Xtilde) %*% Y - t(Xtilde) %*% Xtilde %*%
             (B0 * IndWeights)) / n / lbd
    #A <- which(B0!=0)
    return(list(B0=B0, S0=c(S0), lbd=lbd, weights=weights, group=group))
  } else {
    TEMP <- slassoFit.tilde(Xtilde = Xtilde, Y=Y, lbd=lbd, group=group, weights = weights, verbose = verbose)
    if (sum(TEMP$B0!=0) == (n-1)) {
      warning("Active set is too large. Try to increase the value of lbd.")
    }
    return(list(B0=TEMP$B0, S0=TEMP$S0, sigmaHat=TEMP$hsigma, lbd=lbd, weights=weights, group=group))
  }
}

#' @title Post-inference for lasso estimator
#'
#' @description Provides confidence intervals for the set of active coefficients
#' from lasso estimator using Metropolis-Hastings sampler.
#'
#' @param X predictor matrix.
#' @param Y response vector.
#' @param lbd penalty term of lasso. By letting this argument be \code{"cv.1se"} or
#' \code{"cv.min"}, users can have the cross-validated lambda that gives either minimum
#' squared error or that is within 1 std error bound.
#' @param weights weight vector with length equal to the number of coefficients.
#' Default is \code{rep(1, ncol(X))}.
#' @param tau numeric vector. Standard deviaion of proposal distribution
#'  for each beta. Adjust the value to get relevant level of acceptance rate.
#'  Default is \code{rep(1, ncol(X))}.
#' @param sig2.hat variance of error term.
#' @param alpha confidence level for confidence interval.
#' @param nChain the number of chains. For each chain, different plug-in beta will be generated
#' from its confidence region.
#' @param niterPerChain the number of iterations per chain.
#' @param parallel logical. If \code{parallel = TRUE}, uses parallelization.
#' Default is \code{parallel = FALSE}.
#' @param ncores integer. The number of cores to use for parallelization.
#' @param returnSamples logical. If \code{returnSamples = TRUE}, print Metropolis-Hastings samples.
#' @param ... auxiliary \code{\link{MHLS}} arguments.
#' @details
#' This function provides post-selection inference for lasso estimator.
#' Using Metropolis-Hastings sampler with multiple chains, generates \code{(1-alpha)}
#' confidence interval for each active coefficients.
#' Set \code{returnSamples = TRUE} to check the samples.
#' Check the acceptance rate and adjust \code{tau} accordingly.
#' We recommend to set \code{nChain >= 10} and \code{niterPerChain >= 500}.
#'
#' @return \item{MHsamples}{a list of class MHLS.}
#' @return \item{confidenceInterval}{(1-alpha) confidence interval
#' for each active coefficient.}
#' @examples
#' set.seed(123)
#' n <- 5
#' p <- 10
#' X <- matrix(rnorm(n*p),n)
#' Y <- X %*% rep(1,p) + rnorm(n)
#' sig2 <- 1
#' lbd <- .37
#' weights <- rep(1,p)
#' Postinference.MHLS(X = X, Y = Y, lbd = lbd, sig2.hat = 1, alpha = .05,
#' nChain = 3, niterPerChain = 20, parallel = TRUE)
#' Postinference.MHLS(X = X, Y = Y, lbd = lbd, sig2.hat = 1, alpha = .05,
#' nChain = 3, niterPerChain = 20, parallel = TRUE, returnSamples = TRUE)
#' @export
Postinference.MHLS <- function(X, Y, lbd, weights = rep(1, ncol(X)),
  tau = rep(1, ncol(X)), sig2.hat, alpha = .05, nChain = 10,
  niterPerChain = 500, parallel = FALSE, ncores = 2L, returnSamples=FALSE, ...)
{
  # nChain : the number of MH chains
  # niterPerChain : the number of iteration for each chain
  # B0, S0 : The lasso estimator
  # tau : same as in MHLS function

  LassoEst <- Lasso.MHLS(X=X, Y=Y, type = "lasso", lbd=lbd, weights=weights)
  B0 <- LassoEst$B0
  S0 <- LassoEst$S0
  lbd <- LassoEst$lbd
  A <- which(B0!=0)

  if (length(A)==0) {
    stop("Given lbd, active set is empty.")
  }

  Y <- matrix(Y, , 1)
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)

  nChain <- as.integer(nChain)
  niterPerChain <- as.integer(niterPerChain)

  #--------------------
  # Error Handling
  #--------------------
  parallelTemp <- ErrorParallel(parallel,ncores)
  parallel <- parallelTemp$parallel
  ncores <- parallelTemp$ncores

  if (nrow(X) != nrow(Y)) {
    stop("The dimension of X and Y are not conformable.")
  }
  if (sig2.hat <=0 || lbd <= 0) {
    stop("sig2.hat and/or lbd have to be positive.")
  }
  if (length(weights) != ncol(X)) {
    stop("length(weights) has to be the same with col(X).")
  }
  if (any(weights <= 0)) {
    stop("weights should be positive.")
  }
  if (alpha <=0 || alpha >=1) {
    stop("alpha needs to be between 0 and 1.")
  }
  if (any(c(nChain,niterPerChain) <= 0)) {
    stop("nChain & niterPerChain have to be a positive integer.")
  }
  # if (all.equal(coef(gglasso(X, Y, pf = weights, group = 1:p, loss="ls",
  #                             intercept=F, lambda=lbd))[-1],B0) != TRUE ||
  #     all.equal(c(((t(X)/weights)%*%Y - (t(X) /weights) %*% X %*% B0) / n / lbd)
  #                , S0) != TRUE) {
  #   stop("Invalid B0 or S0, use Lasso.MHLS to get a valid lasso solution.")
  # }
  # Draw samples of pluginbeta from the 95% confidence
  #  region boundary of restricted lse.
  # If nChain ==1, we just use restricted lse.
  if (missing(sig2.hat)) {
    if (length(A) >= nrow(X)) {
      stop("If size of active set matches with nrow(X), sig2.hat needs to be provided.")
    }
    sig2.hat <- summary((lm(Y~X[,A]+0)))$sigma^2
  }

  Pluginbeta.seq <- Pluginbeta.MHLS(X,Y,A,nChain,sqrt(sig2.hat))

  FF <- function(x) {
    MHLS(X = X, PE = Pluginbeta.seq[x,], sig2 = sig2.hat, lbd = lbd,
         weights = weights, niter=niterPerChain,
         burnin = 0, B0 = B0, S0 = S0, tau = tau, PEtype = "coeff", verbose=FALSE, ...)
  }

  if (parallel) {
    TEMP <- parallel::mclapply(1:nChain,FF, mc.cores = ncores)
  } else {
    TEMP <- lapply(1:nChain,FF)
  }
  names(TEMP) <- paste("Chain",1:nChain,sep="")

  MCSAMPLE <- TEMP[[1]]
  if (nChain > 1) {
    for (i in 2:nChain) {
      MCSAMPLE$beta <- rbind(MCSAMPLE$beta, TEMP[[i]]$beta)
      MCSAMPLE$subgrad <- rbind(MCSAMPLE$subgrad, TEMP[[i]]$subgrad)
      MCSAMPLE$acceptHistory <- MCSAMPLE$acceptHistory + TEMP[[i]]$acceptHistory
      MCSAMPLE$niteration <- MCSAMPLE$niteration + TEMP[[i]]$niteration
      MCSAMPLE$burnin <- MCSAMPLE$burnin + TEMP[[i]]$burnin
    }
  }

  # Using MH samples, refit the coeff.
  RefitBeta <- Refit.MHLS(X,weights,lbd,MCSAMPLE)
  if (returnSamples) {
    return(list(MHsamples=TEMP,pluginbeta=Pluginbeta.seq,
                          confidenceInterval=CI.MHLS(betaRefit = RefitBeta,
                          pluginbeta = Pluginbeta.seq, alpha=alpha)))
  } else {
    return(CI.MHLS(betaRefit = RefitBeta, pluginbeta = Pluginbeta.seq,
                   alpha=alpha))
  }
}

# Refit the beta estimator to remove the bias
Refit.MHLS <- function(X,weights,lbd,MHsample) {
  # W : diag(weights)
  # MHsample : MH samples from MHLS function
  n <- nrow(X)
  # Active set
  A <- which(MHsample$beta[1,]!=0)
  # Recover y using KKT condition
  Y.MH <- solve(X%*%t(X))%*%X %*% (crossprod(X) %*% t(MHsample$beta) / n +
                                     lbd * weights * t(MHsample$subgrad)) * n
  # Refit y to the restricted set of X
  beta.refit <- solve(t(X[,A])%*%X[,A])%*%t(X[,A]) %*% Y.MH
  return(beta.refit)
}

# Generate 1-alpha Confidence Interval based on the deviation
CI.MHLS <- function(betaRefit, pluginbeta, alpha=.05) {
  # pluginbeta : a nPlugin x |A| matrix of pluginbeta,
  #  note that the size is not nPlugin x p.
  # beta.refit : refitted beta via Refit.MHLS, a niter x |A| matrix.
  # alpha : significant level.
  A <- which(pluginbeta[1,]!=0)
  Quantile <- apply(betaRefit - pluginbeta[1,A], 1, function(x)
    {quantile(x,prob=c(alpha/2, 1 - alpha/2))})
  Result <- rbind(LowerBound = -Quantile[2,] + pluginbeta[1,A] ,
                  UpperBound = -Quantile[1,] + pluginbeta[1,A])
  colnames(Result) <- paste("beta", A, sep="")
  return(Result)
}

# Generate pluginbeta's from 95% confidence region
Pluginbeta.MHLS <- function(X,Y,A,nPlugin,sigma.hat) {
  # nPlugin : number of pluginbeta's want to generate
  # sigma.hat : estimator of sigma , \epsilon ~ N(0, sigma^2)
  #             If missing, use default way to generate it
  if (length(A)==0) stop("The lasso solution has an empty active set.")

  if (missing(sigma.hat)) {
    sigma.hat <- summary((lm(Y~X[,A]+0)))$sigma
  }
  beta.refit <- coef(lm(Y~X[,A]+0))

  if (nPlugin == 1) {
    return(matrix(beta.refit,1))
  } else {
    xy <- matrix(rnorm(length(A)*(nPlugin-1)), nPlugin-1)
    lambda <- 1 / sqrt(rowSums(xy^2))
    xy <- xy * lambda * sqrt(qchisq(0.95, df=length(A)))
    coeff.seq <- matrix(0,nPlugin,ncol(X))
    coeff.seq[,A]  <- rbind(beta.refit,t(t(xy%*%chol(solve(crossprod(X[,A])))) *
                                           sigma.hat + beta.refit))
    return(coeff.seq=coeff.seq)
  }
}


