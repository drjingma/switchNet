require(glasso)
require(gconcord)

#' To generate a scale free network with p variables.
#' @param p number of variables
#' @param m numeric constant that controls the density of the scale-free network.
#' @details See barabasi.game() for details. Default is NULL, the same as 1.
#' @return  The adjacency matrix A, the graph g and its density.
sf.net <- function(p, m = NULL) {
  g <-
    barabasi.game(
      n = p,
      m = m,
      directed = F,
      algorithm = "psumtree"
    )
  adjm <- as.matrix(get.adjacency(g))
  d <- graph.density(g)
  return(list(A = adjm, g = g, density = d))
}

#' To generate a pd matrix and return the matrix and its inverse
#' @param A A symmetric matrix
pd <- function(A, zeta = 0.01) {
  if (sum(A != t(A)) > 0) {
    stop("This method only works for symmetric A!")
  }

  p <- dim(A)[1]
  diag(A) <- rep(0, p)
  diag(A) <- abs(min(eigen(A)$values)) + zeta
  Ainv <- chol2inv(chol(A))
  Ainv <- cov2cor(Ainv)
  A <- chol2inv(chol(Ainv))
  return(list(A = A, Ainv = Ainv))
}

#' To symmetrize a matrix
#' @param A The initial matrix
symmetrize <- function(A, eps = 1e-06) {
  A <- (A + t(A)) / 2
  A[abs(A) < eps] <- 0
  diag(A) <- 0
  return(A)
}

# log_det <- function(X,eta=0.1){
#   n <- nrow(X)
#   p <- ncol(X)
#   k_v <- 1:p
#   empcov <- 1/n*(t(X)%*%X)
#   if (p<=n){
#     tmp <- sum(digamma((n-k_v+1)/2) - log(n/2))
#   } else {
#     tmp <- -Inf
#   }
#   out <- determinant(empcov+eta*diag(1,p), logarithm = TRUE)$modulus - tmp
# }


# KL_divergence <- function(X, N, delta=0){
#   X <- scale(X)
#   p <- ncol(X)
#   Ip <- diag(rep(1,p))
#   n1 <- N[1]
#   n2 <- N[2]
#   X1 <- X[1:n1,]
#   X2 <- X[-(1:n1),]
#   empcov1 <- 1/n1*(t(X1)%*%X1)
#   empcov2 <- 1/n2*(t(X2)%*%X2)
#   d <- matTr(solve(empcov2+delta*Ip)%*%(empcov1+delta*Ip)) - p + determinant(empcov2+delta*Ip, logarithm = TRUE)$modu - determinant(empcov1+delta*Ip, logarithm = TRUE)$modu
#   return(d/2)
# }



#' Function that describe the transition around the two change points
f <- function(t, tau) {
  if (t < tau[1])
    val <- 0
  else if (t < tau[2])
    val <- (t - tau[1]) / (tau[2] - tau[1])
  else
    val <- 1
  return(val)
}

#' Function that approximates f around the first change point using LQA
g <- function(t, tau, eps) {
  if (t < tau[1] - eps)
    val <- f(t, tau)
  else if (t > tau[1] + eps)
    val <- f(t, tau)
  else
    val <- (t - tau[1] + eps) ^ 2 / (4 * (tau[2] - tau[1]) * eps)
  return(val)
}


#' Function that approximates f around the second change point using LQA
h <- function(t, tau, eps) {
  if (t < tau[2] - eps)
    val <- f(t, tau)
  else if (t > tau[2] + eps)
    val <- f(t, tau)
  else
    val <- 1 - (tau[2] - t + eps) ^ 2 / (4 * (tau[2] - tau[1]) * eps)
  return(val)
}

f.tilde <- function(t, tau, eps) {
  if (t <= tau[1] - eps[1])
    val <- 0
  else if ( (t > tau[1] - eps[1]) && (t <= tau[1] + eps[1]))
    val <- (t - tau[1] + eps[1]) ^ 2 / (4 * (tau[2] - tau[1]) * eps[1])
  else if ( (t > tau[1] + eps[1]) && (t <= tau[2] - eps[2]))
    val <- f(t, tau)
  else if ( (t > tau[2] - eps[2]) && (t <= tau[2] + eps[2]) )
    val <- 1 - (tau[2] - t + eps[2]) ^ 2 / (4 * (tau[2] - tau[1]) * eps[2])
  else
    val <- 1
  return(val)
}

#' @details Calculate the trace of a %*% b; dimensions of a and b should match.
matTr <- function(a,b) {
  if (!identical(dim(a), dim(b))){
    stop("Dimensions do not match!")
  } else {
    crossprod(as.vector(a), as.vector(t(b)))
  }
}

#' #' Calculate the pseudo likelihood using the method in Khare et al. 2015 for a single precision matrix.
#' #' @details This version does not require the inverse covariance matrix to be positive definite.
#' #' @param X The n by p data matrix
#' #' @param Omega A single p by p precision matrix
#' #' @return The pseudo likelihood value
#' pseudo.lkl.single <- function(X, Omega) {
#'   p <- ncol(X)
#'   n <- nrow(X)
#'   h <- function(i) {
#'     temp <- Omega[i, i] * X[, i] + base::crossprod(t(X[, -i]), Omega[-i, i])
#'     return(sum(temp ^ 2) / 2 - n * log(Omega[i, i]))
#'   }
#'   loss <- sapply(1:p, h)
#'   return(sum(loss))
#' }

#' @details This function takes as input a list of precision matrices (\code{Omega}) and returns the pseudo likelihood.
#' The length of \code{Omega} is the same as the number of distinct regimes.
#' @param X A p by N data matrix
#' @param tau The change points
#' @param Omega list(Omega_1,Omega_2)
pseudo.lkl <- function(X, tau, Omega) {
  p <- nrow(X)
  N <- ncol(X)
  K <- length(Omega) # K should be 2.
  Ip <- diag(1,p)

  if (!identical(K,length(tau))){
    stop("Number of change points does not match Number of distinct precision matrices!")
  }
  Omega.sq <- tcrossprod(Omega[[1]])
  Delta.sq <- tcrossprod(Omega[[2]] - Omega[[1]])
  Delta.Omega <- tcrossprod(Omega[[2]] - Omega[[1]],Omega[[1]])

  ind.term <- function(t) {
    S.hat <- tcrossprod(X[,t])
    val <- 0.5*matTr(S.hat,Omega.sq) + 0.5*f(t,tau)^2 * matTr(S.hat,Delta.sq) + f(t,tau)*matTr(S.hat,Delta.Omega)
    val - sum(log(diag(Omega[[1]]) + f(t,tau)*(diag(Omega[[2]]-Omega[[1]]))))
  }

  lkl <-  sum(sapply(1:N, ind.term))/N
  return(lkl)
}

#' The approximated loss based on local quardratic approximation with eps parameters
#' @param X The p by N data matrix
#' @param tau The change points of length 2
#' @param Omega list(Omega_1,Omega_2)
#' @param eps The neighborhood size for approximation, length 2, default to c(5,5)
surrogate.lkl <- function(X, tau, Omega, eps=c(5,5)) {
  p <- nrow(X)
  N <- ncol(X)
  K <- length(Omega) # K should be 2.
  Ip <- diag(1,p)

  if (!identical(K,length(tau))){
    stop("Number of change points does not match Number of distinct precision matrices!")
  }

  Omega.sq <- tcrossprod(Omega[[1]])
  Delta.sq <- tcrossprod(Omega[[2]] - Omega[[1]])
  Delta.Omega <- tcrossprod(Omega[[2]] - Omega[[1]],Omega[[1]])

  ind.term <- function(t) {
    S.hat <- tcrossprod(X[,t])
    val <- 0.5*matTr(S.hat,Omega.sq) + 0.5*f.tilde(t,tau,eps)^2 * matTr(S.hat,Delta.sq) + f.tilde(t,tau,eps)*matTr(S.hat,Delta.Omega)
    val - sum(log(diag(Omega[[1]]) + f.tilde(t,tau,eps)*(diag(Omega[[2]]-Omega[[1]]))))
  }

  lkl <-  sum(sapply(1:N, ind.term))/N
  return(lkl)
}



#-----Optimization Steps-----

#' Function for newton raphlson optimization
#' @param x0 Initial estimation
#' @param lb The lower bound of the search domain
#' @param ub The upper bound of the search domain
#' @param f The objective function
#' @param g The gradient function
#' @param h The hessian function
#' @param alpha Parameter to determine the step size. Default 0.25.
#' @param beta Parameter to determine the step size. Default 0.5.
#' @param maxIter Maximum number of iterations. Default 30.
#' @param tol The tolerance parameter for convergence. Default 0.01.
newton <-
  function(x0,
           lb,
           ub,
           f,
           g,
           h,
           alpha = 0.25,
           beta = 0.5,
           maxIter = 30,
           tol = 0.01,
           verbose = FALSE) {
    count <- 0
    x <- x0

    repeat {
      count <- count + 1

      # print(count)

      ## Newton's step
      delta <- -g(x) / h(x)
      # print(delta)
      ## line search to pick the stepsize
      size <- 1
      while ((x + size * delta <= 0) ||
             (log(f(x + size * delta)) > log(f(x) + alpha * size * g(x) * delta))) {
        size <- beta * size
      }
      # print(size)
      ## Update
      x.new <- x + size * delta

      if (count >= maxIter || abs(x - x.new) < tol || x.new > ub || x.new < lb) {
        if (count == maxIter)
          warning("Maximum number of iterations reached!")
        break
      }
      x <- x.new
      if (verbose) {
        print(count)
      }
    }

    return(list(
      solution = x,
      iter = count,
      stepSize = size
    ))
  }


#' Function to update tau1 using the Newton Raphlson algorithm
#' @param tau0 Initialization, which is a vector of (tau1_0, tau2_0)
#' @param X The p by N data matrix
#' @param Omega A list of precision matrices, length 2
#' @return The initial tau0 and the one-step Newton update.
update.tau1.NR <- function(X, tau0, eps = c(5,5), Omega, maxIter = 30,tol = 5, beta = 0.5,verbose = F) {
  p <- nrow(X)
  N <- ncol(X)

  if (N < 3)
    stop("The time series is too short!")

  if (length(eps)==1){
    eps <- rep(eps,2)
  }

  if (verbose) {
    cat('Currently updating: change point 1...\n')
  }

  Delta <- Omega[[2]]-Omega[[1]]
  Omega.sq <- tcrossprod(Omega[[1]])
  Delta.sq <- tcrossprod(Delta)
  Delta.Omega <- tcrossprod(Delta,Omega[[1]])

  #Function to calculate the first order partial derivative with f.tilde
  pdf <- function(t, tau, eps) {
    d <- (tau[2] - tau[1]) ^ 2
    if ((t > tau[1] - eps[1]) && (t < tau[1] + eps[1] + 1)){
      val <- (t - tau[1] + eps[1]) * (t + tau[1] + eps[1] - 2 * tau[2]) / (4 * eps[1] * d)
    } else if ( (t > tau[1] + eps[1]) && (t < tau[2] - eps[2] + 1))
      val <- (t - tau[2]) / d
    else if ( (t > tau[2] - eps[2]) && (t < tau[2] + eps[2] + 1))
      val <- - (tau[2] - t + eps[2])^2 / (4 * eps[2] * d)
    else
      val <- 0
    return(val)
  }

  #Function to calculate the Second order partial derivative with f.tilde
  pdf2 <- function(t, tau, eps) {
    d <- (tau[2] - tau[1]) ^ 3
    if ( (t > tau[1] - eps[1]) && (t < tau[1] + eps[1] + 1) ){
      val <- (tau[2] - t - eps[1])^2 / (2 * eps[1] * d)
    } else if ( (t > tau[1] + eps[1]) && (t < tau[2] - eps[2] + 1) )
      val <- 2 * (t - tau[2]) / d
    else if ( (t > tau[2] - eps[2]) && (t < tau[2] + eps[2] + 1))
      val <- - (tau[2] - t + eps[2])^2 / (2 * eps[2] * d)
    else
      val <- 0
    return(val)
  }

  #Function to calculate the First order derivative of the loss function
  pdL <- function(tau1) {
    deriv <- function(t) {
      S.hat <- tcrossprod(X[,t])
      val <- matTr(S.hat, Delta.Omega) + f.tilde(t,c(tau1,tau0[2]),eps) * matTr(S.hat,Delta.sq) -
        sum( diag(Delta)/(diag(Omega[[1]])+f.tilde(t,c(tau1,tau0[2]),eps)*diag(Delta)) )

      val <- val * pdf(t,c(tau1,tau0[2]),eps)
      return(val)
    }
    res <- sum(sapply(1:N, deriv))/N

    return(res)
  }


  #Function to calculate the Second order derivative of the loss function
  pdL2 <- function(tau1) {
    deriv <- function(t) {
      S.hat <- tcrossprod(X[,t])
      val <- matTr(S.hat, Delta.Omega) + f.tilde(t,c(tau1,tau0[2]),eps) * matTr(S.hat,Delta.sq) -
        sum( diag(Delta)/(diag(Omega[[1]])+f.tilde(t,c(tau1,tau0[2]),eps)*diag(Delta)) )
      val <- val * pdf2(t,c(tau1,tau0[2]),eps)
      val <- val + matTr(S.hat,Delta.sq) * pdf(t,c(tau1,tau0[2]),eps) * pdf(t,c(tau1,tau0[2]),eps)
      val <- val + pdf(t,c(tau1,tau0[2]),eps) * pdf(t,c(tau1,tau0[2]),eps) * sum(diag(Delta)^2/(diag(Omega[[1]])+f.tilde(t,c(tau1,tau0[2]),eps)*diag(Delta))^2)
      return(val)
    }
    res <- sum(sapply(1:N, deriv))/N

    return(res)
  }

  ## One step Newton update
  z_old <- tau0[1]
  # z_new <- newton(z_old, 5, tau0[2], L, pdL, pdL2)$solution
  z_new <- z_old - beta*pdL(z_old) / pdL2(z_old)
  cat("2nd order derivative wrt tau1 is ...",pdL2(z_old),"...\n")
  # print(z_new)
  if (z_new < 1) {
    stop('The estimated change point tau_1 is negative!')
  }

  z_old <- z_new

  return(list(tau = z_old, tau0 = tau0))

}

#' Function to update tau2 using the Newton Raphlson algorithm
#' @param tau0 Initialization, which is a vector of (tau1_0, tau2_0)
#' @param X The p by N data matrix
#' @param Omega A list of precision matrices, length 2, each of p by p
#' @return The initial tau0 and the one-step Newton update
update.tau2.NR <- function(X,tau0,eps = c(5,5),Omega,maxIter = 30,tol=5,beta=0.5,verbose=F) {
  p <- nrow(X)
  N <- ncol(X)

  if (N < 3)
    stop("The time series is too short!")

  if (length(eps)==1){
    eps <- rep(eps,2)
  }

  if (verbose) {
    cat('Currently updating: change point 1...\n')
  }

  Delta <- Omega[[2]]-Omega[[1]]
  Omega.sq <- tcrossprod(Omega[[1]])
  Delta.sq <- tcrossprod(Delta)
  Delta.Omega <- tcrossprod(Delta,Omega[[1]])

  #Function to calculate the first order partial derivative with f.tilde
  pdf <- function(t, tau, eps) {
    d <- (tau[2] - tau[1]) ^ 2
    if ( (t > tau[1] - eps[1]) && (t < tau[1] + eps[1] + 1)){
      val <- - (t - tau[1] + eps[1]) ^2 / (4 * eps[1] * d)
    } else if ( (t > tau[1] + eps[1]) && (t < tau[2] - eps[2] + 1))
      val <- - (t - tau[1]) / d
    else if ( (t > tau[2] - eps[2]) && (t < tau[2] + eps[2] + 1))
      val <-  (tau[2] - t + eps[2])  * (-tau[2] - t + eps[2] + 2*tau[1]) / (4 * eps[2] * d)
    else
      val <- 0
    return(val)
  }

  #Function to calculate the Second order partial derivative with f
  pdf2 <- function(t, tau, eps) {
    d <- (tau[2] - tau[1]) ^ 3
    if ( (t > tau[1] - eps[1]) && (t < tau[1] + eps[1] + 1)){
      val <- (t - tau[1] + eps[1])^2 / (2*eps[1] *d)
    } else if ( (t > tau[1] + eps[1]) && (t < tau[2] - eps[2] + 1))
      val <- 2 * (t - tau[1]) / d
    else if ( (t > tau[2] - eps[2]) && (t < tau[2] + eps[2] + 1))
      val <- - (t - tau[1] - eps[2])^2 / (2*eps[2] *d)
    else
      val <- 0
    return(val)
  }


  #Function to calculate the First order derivative of the loss function
  pdL <- function(tau2) {
    deriv <- function(t) {
      S.hat <- tcrossprod(X[,t])
      val <- matTr(S.hat, Delta.Omega) + f.tilde(t,c(tau0[1],tau2),eps) * matTr(S.hat,Delta.sq) -
        sum( diag(Delta)/(diag(Omega[[1]])+f.tilde(t,c(tau0[1],tau2),eps)*diag(Delta)) )

      val <- val * pdf(t,c(tau0[1],tau2),eps)
      return(val)
    }
    res <- sum(sapply(1:N, deriv))/N

    return(res)
  }

  # Function to calculate the Second order derivative of the loss function
  pdL2 <- function(tau2) {
    deriv <- function(t) {
      S.hat <- tcrossprod(X[,t])
      val <- matTr(S.hat, Delta.Omega) + f.tilde(t,c(tau0[1],tau2),eps) * matTr(S.hat,Delta.sq) -
        sum( diag(Delta)/(diag(Omega[[1]]) + f.tilde(t,c(tau0[1],tau2),eps) * diag(Delta)) )
      val <- val * pdf2(t,c(tau0[1],tau2),eps)
      val <- val + matTr(S.hat,Delta.sq) * pdf(t,c(tau0[1],tau2),eps) * pdf(t,c(tau0[1],tau2),eps)
      val <- val + pdf(t,c(tau0[1],tau2),eps) * pdf(t,c(tau0[1],tau2),eps) * sum(diag(Delta)^2/(diag(Omega[[1]])+f.tilde(t,c(tau0[1],tau2),eps)*diag(Delta))^2)
      return(val)
    }
    res <- sum(sapply(1:N, deriv))/N

    return(res)
  }

  ## One-step newton update
  z_old <- tau0[2]
  # z_new <- newton(z_old, tau0[1], N-5, L, pdL, pdL2)$solution
  z_new <- z_old - beta*pdL(z_old) / pdL2(z_old)
  cat("2nd order derivative wrt tau2 is ...",pdL2(z_old),"...\n")
  # print(z_new)
  if (z_new > N) {
    stop('The estimated change point 2 is out of bounds!')
  }

  z_old <- z_new
  return(list(tau = z_old, tau0 = tau0))
}



#' For fixed change points, compute the partial correlations in regime 1 and 2 using Concord.
#' Compute the partial correlations in the intermediate regime using only nearby observations.
#' @param X The p by N data matrix
#' @param tau The change points c(tau1, tau2)
#' @param eps A vector of neighborhood sizes (eps1, eps2) used in estimating the change points
#' @param lambda A vector of tuning parameters (lambda1, lambda2) for estimating the two inverse covariances
#'               with Concord
#' @return The estimated precision matrices
update.SigInv.concord <- function(
  X,
  tau,
  eps,
  lambda=NULL
) {
  # pause(0.1)
  N <- ncol(X)
  p <- nrow(X)
  Ip <- diag(1, p)

  if (is.null(lambda)){
    lambda <- c(sqrt(tau[1] * log(N * p * (p + 1) / 2)) / N, sqrt((N - tau[2]) * log(N * p * (p + 1) / 2)) / N)
  }

  listX <- vector("list",length(tau))
  listX[[1]] <- t(X[,1:(tau[1]-eps[1])])
  listX[[2]] <- t(X[,(1 + tau[2] + eps[2]):N])

  Omega <- lapply(1:length(tau), function(j) concord(scale(listX[[j]]), lambda[j], maxit=100))
  return(Omega)
}


#' Find the initial estimate of change points using NR optimization
#' @param  X  A p by N data matrix
#' @param  tau0 Initialization of the change points
#' @param  eps Step size for approximating the loss function
#' @param  lambda The tuning parameters used c(lam1, lam2) for estimating partial correlations.
#' @param  type Method for calculating the loss. Default "approx".
#' @param  maxIter Maximum number of iteration. Default 20.
#' @param  tol Tolerance level for changes in loss.
find.changepoints.NR <-
  function(X,
           tau0,
           eps = c(5,5),
           lambda = NULL,
           type = c("pseudo","approx"),
           maxIter = 10,
           tol = 0.01,
           verbose=FALSE
  ) {
    # pause(0.1)
    type <- match.arg(type)

    N <- ncol(X)
    p <- nrow(X)
    X <- t(scale(t(X)))

    iter <- 0
    tau_new <- c(0,0)
    tau_old <- tau0
    Omega_old <- update.SigInv.concord(X, tau_old, eps, lambda)
    loss_old <- 10^8
    change_loss <- 100
    while (change_loss > tol && iter < maxIter) {
      if (verbose){
        cat("iter...", iter, "...\n")
        cat("loss", loss_old, "...\n")
        cat("tau", tau_old, "...\n")
      }
      ##Update tau
      tau_new[1] <- update.tau1.NR(X, floor(tau_old), eps, Omega_old, beta=1)$tau
      tau_new[2] <- update.tau2.NR(X, floor(c(tau_new[1], tau_old[2])), eps, Omega_old, beta=1)$tau
      tau_new <- floor(tau_new)

      ##Update SigInv
      Omega_new <- update.SigInv.concord(X, tau_new, eps)

      ##Change in loss
      ##We expect the loss to decrease and this quantity to be positive
      loss_new <- switch(type,
                         approx = surrogate.lkl(X,tau_new,eps,Omega_new),
                         pseudo = pseudo.lkl(X,tau_new,Omega_new))
      change_loss <- loss_old - loss_new
      cat("change in loss is...", change_loss, "...\n")

      if (change_loss>tol){
        tau_old <- tau_new
        loss_old <- loss_new
        Omega_old <- Omega_new
      } else{ # change in loss is negative
        tau_old <- tau_old
        loss_old <- loss_old
        Omega_old <- Omega_old
      }

      # ##Update lambda
      # lambda <- c(sqrt(tau_old[1] * log(N * p * (p + 1) / 2)) / N, sqrt((N - tau_old[2]) * log(N * p * (p + 1) / 2)) / N)
      #
      iter <- iter + 1
    }

    out <- list(
      tau = round(tau_old),
      fit = loss_old,
      Omega = Omega_old,
      Iter = iter
    )
    return(out)
  }



#' Conventional search over the 2-d grid for optimal tau
#' @param X A p by N data matrix
#' @param grid.tau1 The grid for the first change point
#' @param grid.tau2 The grid for the second change point
find.changepoints.search <- function(
  X,
  grid.tau1,
  grid.tau2,
  eps=c(5,5),
  type = c("pseudo","approx")
) {
  # pause(0.1)
  n1 <- length(grid.tau1)
  n2 <- length(grid.tau2)
  p <- nrow(X)
  X <- t(scale(t(X)))
  type <- match.arg(type)

  loss <- matrix(0, n1, n2)
  for (loop.grid1 in 1:n1) {
    for (loop.grid2 in 1:n2){
      tau <- c(grid.tau1[loop.grid1], grid.tau2[loop.grid2])
      Omega <- update.SigInv.concord(X, tau, eps=c(0,0))
      loss[loop.grid1,loop.grid2] <- switch(type,
                                            approx = surrogate.lkl(X,tau,eps,Omega),
                                            pseudo = pseudo.lkl(X,tau,Omega))
    }
  }

  return(loss)
}

#----- Regime switching function to estimate the change points ----
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
switchNet.cp <-
  function(X,
           tau0,
           eps = c(5,5),
           lambda = NULL,
           method = c("newton", "search"),
           maxIter = 30,
           search = FALSE,
           grid.tau1 = NULL,
           grid.tau2 = NULL
  ) {
    # pause(0.1)
    N <- ncol(X)
    p <- nrow(X)

    method <- match.arg(method)

    if (method == "newton") {
      r <- find.changepoints.NR(X, tau0, eps = eps, lambda=lambda, type="pseudo", maxIter = maxIter)
      tau.hat <- r$tau
      fit <- r$fit
      if (search){
        grid.tau1 <- seq(r$tau[1] - eps[1], r$tau[1] + eps[1])
        grid.tau2 <- seq(r$tau[2] - eps[1], r$tau[2] + eps[1])
        fit <- find.changepoints.search(X, grid.tau1, grid.tau2)
        temp <- which(fit == min(fit), arr.ind = TRUE)
        tau.hat <- c(grid.tau1[temp[1]], grid.tau2[temp[2]])
      }
    } else {
      #Search based method
      fit <- find.changepoints.search(X, grid.tau1, grid.tau2)
      temp <- which(fit == min(fit), arr.ind = TRUE)
      tau.hat <- c(grid.tau1[temp[1]], grid.tau2[temp[2]])
    }
    return(list(tau=tau.hat,fit=fit))
  }



#-----Model Selection-----
##Once tau.hat is available, refit the model to obtain the precision matrices
#' @param  X A p by N data matrix
#' @param  tau Estimate of tau
#' @param  lambda Optimal tuning parameters for each regime
#' @export
switchNet <- function(X,
                      tau,
                      lambda = NULL,
                      cstar = c(1, 1) ) {
  # pause(0.1)
  p <- nrow(X)
  N <- ncol(X)
  X <- scale(t(X))

  listX <- vector("list", length(tau))
  listX[[1]] <- X[,1:tau[1]]
  listX[[2]] <- X[,(1 + tau[2]):N]

  if (is.null(lambda)) {
    lambda <- cstar * c(sqrt(tau[1] * log(N * p * (p + 1) / 2)) / N, sqrt((N - tau[2]) * log(N * p * (p + 1) / 2)) / N)
  }

  Omega.hat <- lapply(1:length(tau), function(j) concord(t(listX[[j]]), lambda[j], maxit = 30))

  return(Omega.hat)
}
#' Compute the bic.score over a grid of lambda_1 and lambda_2
#' @param X
#' @param tau
#' @param lambda1 A candidate sequence of tuning parameters for regime 1
#' @param lambda2 A candidate sequence of tuning parameters for regime 2
#' @param tol
#' @export
switchNet.bic <- function(
  X,
  tau,
  lambda1,
  lambda2,
  tol = 1e-06
) {
  p <- nrow(X)
  N <- ncol(X)
  Ip <- diag(p)
  listX <- vector("list", length(tau))
  listX[[1]] <- X[,1:tau[1]]
  listX[[2]] <- X[,(1 + tau[2]):N]

  bic.score <- matrix(0, length(lambda1), length(lambda2))

  for (loop_lam1 in 1:length(lambda1)) {
    cat("The ", loop_lam1, "-th step in tuning... \n")
    for (loop_lam2 in 1:length(lambda2)) {
      for (loop_regime in 1:nregime) {
        X1 <- scale(listX[[loop_regime]])

        rho <- lambda1[loop_lam1]
        if (loop_regime == 2) {
          rho <- lambda2[loop_lam2]
        }

        sigInv <- concord(X1, rho, maxit = 30)
        # sigInv <- solve(cov2cor(fit1))
        no.edge <- sum(abs(sigInv) > tol) - p
        bic.score[loop_lam1, loop_lam2] <- bic.score[loop_lam1, loop_lam2] +
          matTr(empcov, sigInv) - log(det(sigInv)) +
          log(nrow(X1)) * no.edge /(2 * nrow(X1))
      }

    }
  }
  return(bic.score)
}


#' Metrics to compare the estimated and the truth matrices.
#' @param Ahat a list of the estimated matrices.
#' @param Amat a list of the corresponding true matrices.
#' @return
#'  \item{FPrate}{false positive rate}
#'  \item{FNrate}{false negative rate}
#'  \item{SHD}{structural hamming distance}
#'  \item{Floss}{Frobenius norm loss}
#' @details  The above measures are as defined in the paper, which are the average deviance
#' from all pairs of matrices. The following is specifically for calculating the ROC curve.
#'
StructDiff <- function(Ahat, Amat, eps = 1e-06) {
  K <- length(Amat)

  if (is.null(dim(Amat)) == F) {
    # indicating there's only one matrix compared.
    K <- 1
    Ahat <- list(Ahat)
    Amat <- list(Amat)
  }

  p <- dim(Amat[[1]])[1]
  TP <- rep(0, K)
  FP <- rep(0, K)
  TN <- rep(0, K)
  FN <- rep(0, K)
  SHD <- rep(0, K)
  FPrate <- rep(0, K)
  TPrate <- rep(0, K)
  FNrate <- rep(0, K)
  Pr <- rep(0, K)
  Re <- rep(0, K)
  F1 <- rep(0, K)
  Floss <- rep(0, K)
  for (k in 1:K) {
    Floss[k] <- Matrix::norm(Ahat[[k]] - Amat[[k]], "F") / Matrix::norm(Amat[[k]], "F")
    diag(Ahat[[k]]) <- 0
    diag(Amat[[k]]) <- 0

    TP[k] <- sum((abs(Ahat[[k]]) > eps) * (abs(Amat[[k]]) > eps))
    TN[k] <- sum((abs(Ahat[[k]]) <= eps) * (abs(Amat[[k]]) <= eps))
    FP[k] <- sum((abs(Ahat[[k]]) > eps) * (abs(Amat[[k]]) <= eps))
    FN[k] <- sum((abs(Ahat[[k]]) <= eps) * (abs(Amat[[k]]) > eps))
    SHD[k] <- FP[k] + FN[k]

    P <- TP[k] + FN[k]
    N <- TN[k] + FP[k]
    TPrate[k] <- TP[k] / (P + eps)
    FPrate[k] <- FP[k] / (N + eps)
    FNrate[k] <- FN[k] / (P + eps)

    Re[k] <- TP[k] / (P + eps) ## Recall
    Pr[k] <- TP[k] / (TP[k] + FP[k] + eps) ## Precision

    F1[k] <- (2 * Pr[k] * Re[k]) / (Pr[k] + Re[k] + eps)

  }

  dev <- data.frame(
    TP = TP,
    FP = FP,
    TN = TN,
    FN = FN,
    TPrate = TPrate,
    FPrate = FPrate,
    FNrate = FNrate,
    SHD = SHD,
    Precision = Pr,
    Recall = Re,
    F1 = F1,
    FL = Floss
  )
  return(dev)
}
