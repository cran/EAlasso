#' @title Compute K-fold cross-validated mean squared error for lasso
#'
#' @description Computes K-fold cross-validated mean squared error
#' to propose a lambda value for lasso, group lasso, scaled lasso or scaled
#' group lasso.
#'
#' @param X predictor matrix.
#' @param Y response vector.
#' @param group \code{p} x \code{1} vector of consecutive integers describing the group structure.
#' The number of groups should be the same as max(group). Default is \code{group = 1:p}
#' , where \code{p} is number of covariates. See examples for a guideline.
#' @param weights weight vector with length equal to the number of groups. Default is
#' \code{rep(1, max(group))}.
#' @param type type of penalty. Must be specified to be one of the following:
#'  \code{"lasso", "grlasso", "slasso"} or \code{"sgrlasso"}.
#' @param K integer. Number of folds
#' @param minlbd numeric. Minumum value of the lambda sequence.
#' @param maxlbd numeric. Maximum value of the lambda sequence.
#' @param num.lbdseq integer. Length of the lambda sequence.
#' @param parallel logical. If \code{parallel = TRUE}, uses parallelization.
#' Default is \code{parallel = FALSE}.
#' @param ncores integer. The number of cores to use for parallelization.
#' @param plot.it logical. If true, plots the squared error curve.
#' @param verbose logical.
#'
#' @return \item{lbd.min}{a value of lambda which gives a minimum squared error.}
#' @return \item{lbd.1se}{a largest lambda within 1 standard error from \code{lbd.min}.}
#' @return \item{lbd.seq}{lambda sequence.}
#' @return \item{cv}{mean squared error at each lambda value.}
#' @return \item{cvsd}{the standard deviation of cv.}
#'
#' @examples
#' set.seed(123)
#' n <- 30
#' p <- 50
#' group <- rep(1:(p/10),each=10)
#' weights <- rep(1, max(group))
#' X <- matrix(rnorm(n*p),n)
#' truebeta <- c(rep(1,5),rep(0,p-5))
#' Y <- X%*%truebeta + rnorm(n)
#'
#' # To accelerate the computational time, we set K=2 and num.lbdseq=2.
#' # However, in practice, Allowing K=10 and num.lbdseq > 100 is recommended.
#' cv.lasso(X = X, Y = Y, group = group, weights = weights, K = 2,
#' type = "grlasso", num.lbdseq = 2, plot.it = FALSE)
#' cv.lasso(X = X, Y = Y, group = group, weights = weights, K = 2,
#' type = "sgrlasso", num.lbdseq = 2, plot.it = FALSE)
#' @export
cv.lasso <- function(
  X,
  Y,
  group = 1:ncol(X),
  weights = rep(1,max(group)),
  type,
  K = 10L,
  minlbd,
  maxlbd,
  num.lbdseq = 100L,
  parallel = FALSE,
  ncores = 2L,
  plot.it = FALSE,
  verbose = FALSE)
{
  n <- nrow(X)
  p <- ncol(X)

  Y <- as.vector(Y)
  K <- as.integer(K)
  num.lbdseq <- as.integer(num.lbdseq)

  if (!parallel) {ncores <- 1}

  if(missing(minlbd)) {minlbd <- 0}
  if(missing(maxlbd)) {
    maxlbd <- if (type == "lasso")
      {
        max(abs(1/weights * t(X) %*% Y))/n
      } else if (type == "grlasso")
      {
        lbdTEMP <- c()
        for (i in 1:max(group)) {
          lbdTEMP[i] <- sqrt(crossprod(t(X[,group==i])%*%Y))/weights[i]
        }
        max(lbdTEMP / n)
      } else {2}
  }
  #--------------------
  # Error Handling
  #--------------------
  if (n < 10) {stop("Sample size is too small to run cross-validation.")}
  if (n < 30) {warning("Insufficient sample size. The result can be unstable.")}

  if (length(Y) != n) {
    stop("dimension of X and Y are not conformable.")
  }
  if (!type %in% c("lasso", "grlasso", "slasso", "sgrlasso")) {
    stop("type has to be either lasso, grlasso, slasso or sgrlasso.")
  }
  if (length(group) != p) {
    stop("group must have a same length with the number of X columns")
  }
  if (!all(group==1:p) && (!type %in% c("grlasso", "sgrlasso"))) {
    stop("Choose type to be either grlasso or sgrlasso if group-structure exists.")
  }
  parallelTemp <- ErrorParallel(parallel,ncores)
  parallel <- parallelTemp$parallel
  ncores <- parallelTemp$ncores
  if (length(weights) != length(unique(group))) {
    stop("weights has to have a same length as the number of groups")
  }
  if (any(weights <= 0)) {
    stop("weights should be positive.")
  }
  if (any(!1:max(group) %in% group)) {
    stop("group index has to be a consecutive integer starting from 1.")
  }
  if (any(c(minlbd,maxlbd) < 0)) {stop("minlbd/maxlbd should be non-negative")}
  if (minlbd >= maxlbd) {stop("minlbd is too large compared to maxlbd.")}
  if (num.lbdseq <= 0) {
    stop("num.lbdseq should be non-negative")
  }
  if (K <= 0) {
    stop("K should be a positive integer.")
  }

  all.folds <- split(sample(1:n), rep(1:K, length = n))

  #index=seq(0,max(t(X)%*%Y),length=100)
  index <- seq(minlbd,maxlbd,length=num.lbdseq)[-1]
  residmat <- matrix(0, length(index), K)

  FF <- function(x,omit) {
    fit <- Lasso.MHLS(X=X[-omit,,drop=FALSE],Y=Y[-omit],type=type,
                      lbd=index[x],group=group,weights=weights)$B0
    fit <- X[omit,,drop=FALSE]%*%fit
    return(mean((Y[omit]-fit)^2))
  }

  for (i in seq(K)) {
    omit <- all.folds[[i]]
    residmat[,i] <- do.call(rbind,parallel::mclapply(1:length(index),
      FF,omit = omit, mc.cores = ncores))
    if(verbose) {cat("\n CV Fold", i, "\n\n")}
  }
  #apply(residmat, 2, which.min)
  cv <- apply(residmat, 1, mean)
  cvsd <- apply(residmat, 1, sd)
  index.min.cv <- which.min(cv)

  err.1se <- cvsd[index.min.cv] + cv[index.min.cv]

  index.min.cv:num.lbdseq

  lbd.1se <- index[which.min(abs(cv - err.1se)[index.min.cv:num.lbdseq]) + index.min.cv - 1]
  lbd.min <- index[index.min.cv]

  if (plot.it) {
    matplot(x=index, y=cbind(cv,cv-cvsd,cv+cvsd), type="l",
            lty=c(1,2,2), col=c(1,2,2),
            xlab="lambda",ylab="mean squared error",
            main="cross validation")
    abline(v=lbd.min,lty=2)
    abline(v=lbd.1se,lty=3)
  }
  return(list(lbd.seq=index, cv=cv, cvsd=cvsd,
              lbd.min=lbd.min, lbd.1se=lbd.1se,
              type=type))
}
