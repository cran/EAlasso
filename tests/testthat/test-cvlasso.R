
context("cross validation")

set.seed(1234)
n <- 30
p <- 50
Niter <-  10
group <- rep(1:(p/10), each = 10)
weights <- rep(1, p/10)
X <- matrix(rnorm(n*p), n)
Y <- X%*%rep(1,p) + rnorm(n)
parallel <- (.Platform$OS.type != "windows")

test_that("invalid parameters", {
  expect_error(cv.lasso(X = X, Y = Y, group = group, weights = weights, K=-1,
    type="sgrlasso", num.lbdseq=10, parallel = parallel),
    "positive integer")
  expect_error(cv.lasso(X = X, Y = Y, group = group, weights = weights, K=5,
                        minlbd = 2, maxlbd = 1,
                        type="sgrlasso", num.lbdseq=10, parallel = parallel),
               "too large")
  expect_error(cv.lasso(X = X, Y = Y, group = group, weights = weights, K=5,
                        minlbd = -2, maxlbd = 1, num.lbdseq = 10,
                        type="sgrlasso",  parallel = parallel),
               "non-negative")
  expect_error(cv.lasso(X = X, Y = Y, group = group, weights = weights, K=5,
                        num.lbdseq = -10,
                        type="sgrlasso", parallel = parallel),
               "non-negative")
})

test_that("NA", {
  expect_error(cv.lasso(X = X, Y = Y, group = group, weights = weights, K=3,
                        type="sgrlasso", num.lbdseq=10, parallel = parallel)
               , NA)
  expect_error(cv.lasso(X = X, Y = Y, group = group, weights = weights, K=3,
                        type="grlasso", num.lbdseq=10, parallel = parallel)
               , NA)
  expect_error(cv.lasso(X = X, Y = Y, K=3,
                        type="slasso", num.lbdseq=10, parallel = parallel)
               , NA)
  expect_error(cv.lasso(X = X, Y = Y, K=3,
                        type="lasso", num.lbdseq=10, parallel = parallel)
               , NA)
})
