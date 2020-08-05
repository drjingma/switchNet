
#' Regime switching function to estimate the change points
#' @param  X A p by N data matrix
#' @param  tau0 initialization of tau
#' @param  eps nbhd size of approximation
#' @param  lambda The tuning parameters in a vector of length 2
#' @param  method The optimization method used. Could be "Newton" or "search",
#' @param  maxIter The maximum number of iterations
#' @param  search Whether to perform grid search after Newton's update
#' @param  grid.tau1
#' @param  grid.tau2 the search grid for tau1 and tau2. Must be provided if not using NR optimization
#' @export
findChangePoints <-
  function(X,
           tau0,
           eps = c(5,5),
           lambda = NULL,
           method = c("newton", "search"),
           maxIter = 30,
           both = TRUE,
           grid.tau1 = NULL,
           grid.tau2 = NULL
  ) {
    # pause(0.1)
    N <- ncol(X)
    p <- nrow(X)

    method <- match.arg(method)

    if (method == "newton") {
      r <- find.cp.NR.only(X, tau0, eps = eps, lambda=lambda, type="pseudo", maxIter = maxIter)
      tau.hat <- r$tau
      fit <- r$fit
      if (both){
        grid.tau1 <- seq(r$tau[1] - eps[1], r$tau[1] + eps[1])
        grid.tau2 <- seq(r$tau[2] - eps[1], r$tau[2] + eps[1])
        fit <- find.cp.search(X, grid.tau1, grid.tau2)
        temp <- which(fit == min(fit), arr.ind = TRUE)
        tau.hat <- c(grid.tau1[temp[1]], grid.tau2[temp[2]])
      }

    } else {
      #Search based method
      fit <- find.cp.search(X, grid.tau1, grid.tau2, type = 'pseudo')
      temp <- which(fit == min(fit), arr.ind = TRUE)
      tau.hat <- c(grid.tau1[temp[1]], grid.tau2[temp[2]])
    }
    return(list(tau=tau.hat,fit=fit))
  }



#'   Once tau.hat is available, refit the model to obtain the precision matrices
#' @param  X A p by N data matrix
#' @param  tau Estimate of tau
#' @param  lambda Optimal tuning parameters for each regime
#' @export
switchNet <- function(X, tau, lambda = NULL, cstar = c(1, 1) ) {
  # pause(0.1)
  p <- nrow(X)
  N <- ncol(X)
  X <- t(scale(t(X)))

  listX <- vector("list", length(tau))
  listX[[1]] <- X[,1:tau[1]]
  listX[[2]] <- X[,(1 + tau[2]):N]

  if (is.null(lambda)) {
    lambda <- cstar * c(sqrt(tau[1] * log(N * p * (p+1)/2)) / N, sqrt((N - tau[2]) * log(N * p * (p+1)/2)) / N)
  }

  Omega.hat <- lapply(1:length(tau), function(j) concord(t(listX[[j]]), lambda[j], maxit = 30))

  return(Omega.hat)
}
#'   Compute the bic.score over a grid of lambda_1 and lambda_2
#' @param X A p by N data matrix
#' @param tau A length 2 vector of change points
#' @param lambda1 A candidate sequence of tuning parameters for regime 1
#' @param lambda2 A candidate sequence of tuning parameters for regime 2
#' @param tol Tolerance for declaring an edge
#' @details This function
#' @export
switchNetBIC <- function(X,  tau,  lambda1,  lambda2,  tol = 1e-06) {
  p <- nrow(X)
  N <- ncol(X)
  Ip <- diag(p)
  nregime <- length(tau)
  listX <- vector("list", nregime)
  listX[[1]] <- X[,1:tau[1]]
  listX[[2]] <- X[,(1 + tau[2]):N]

  bic.score <- matrix(0, length(lambda1), length(lambda2))

  for (loop_lam1 in 1:length(lambda1)) {
    cat("The ", loop_lam1, "-th step in tuning... \n")
    for (loop_lam2 in 1:length(lambda2)) {
      for (loop_regime in 1:nregime) {
        X1 <- scale(t(listX[[loop_regime]]))
        empcov <- crossprod(X1)/(nrow(X1)-1)
        rho <- sqrt(nrow(X1) * log(N * p * (p+1)/2)) / N
        rho <- ifelse(loop_regime==1, rho*lambda1[loop_lam1],rho*lambda2[loop_lam2])

        sigInv <- concord(X1, rho, maxit = 30)
        # sigInv <- solve(cov2cor(fit1))
        no.edge <- sum(abs(sigInv) > tol) - p
        bic.score[loop_lam1, loop_lam2] <- bic.score[loop_lam1, loop_lam2] +
          matTr(empcov, sigInv) - log(det(sigInv)) + log(nrow(X1)) * no.edge /(2 * nrow(X1))
      }

    }
  }
  opt.index <- which(bic.score==min(bic.score),arr.ind = TRUE)
  return(list(lambda=c(lambda1[opt.index[1]],lambda2[opt.index[2]]),bic.score=bic.score))
}


