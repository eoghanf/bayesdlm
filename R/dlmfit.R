#------------------------------------------------------------------------------
#
#                                 bayesdlm
#
#                                    by
#
#                              Eoghan Flanagan
#
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------

#' Checks whether a provided value is a scalar.
#'
#' \code{is.scalar} returns true if the provided argument is a scalar and false
#' otherwise.
#'
#' @param x is an object of unknown type
#' @examples is.scalar('A')
#'is.scalar(matrix(5))
#'is.scalar(sum)


is.scalar <- function(x){
  if(is.vector(x) & length(x)==1 & is.numeric(x)){
    return(TRUE)
  }else{
    return(FALSE)
  }
}
#------------------------------------------------------------------------------

#' Helper function which converts arguments and result of DLM fitting
#' into a function suitable for passing to optim
#'
#' \code{is.scalar} returns the negative log-likelihood of the data between time
#' t.start and t.stop under the DLM model with W discount factor W.disc.factor
#' and observation volatility v contained in the vector x
#'
#' @param x is a two element vector containing W.disc.factor and v
#' @param y is
#' @param t.start
#' @param t.stop
#' @examples
#'
#'


llike.hyperparams <- function(hyperparams, data, t.start=1, t.stop=length(y)){
  dlm.output=dlmfilter(y=data$y, X=data$X, W.disc.factor=hyperparams[1],
                   v=hyperparams[2], prior.dist = "N", prior.mean=0, prior.var=10e+7,
                   ng.alpha=1, ng.beta=1, calc.hist.vars=FALSE)

  if (is.finite(-sum(dlm.output$llike[t.start:t.stop]))){

    return(-sum(dlm.output$llike[t.start:t.stop]))}

  else{
    stop("The MLE calculation diverged. Sorry.")
  }

  }


#------------------------------------------------------------------------------
#
#' \code{prodnorm} returns the mean and variance of the normal PDF which is the
#product of two Normal PDFs
#' @param mu1 is a vector of means
#' @param sigma is a covariance matrix
#' @param m2 is a vector of means
#' @param sigma2 is a covariance matrix '
#'
#' @return A list [mu, sigma] of the mean and covariance matrix of the product PDF

prodnorm <- function(mu1, sigma1, mu2, sigma2){

  if (!(is.vector(mu1) & is.vector(mu2))){
    stop('mu1 and mu2 must be vectors of means')
  }

  if (!(is.numeric(mu1) & is.numeric(mu2))){
    stop('mu1 and mu2 must be vectors of means')
  }

  if (!(length(mu1) == length(mu2))){
    stop('mu1 and mu2 must have the same length')
  }

  n <- length(mu1)

  if (!(nrow(sigma1) == n & ncol(sigma1)==n)){
    stop(paste('sigma 1 must be a square matrix of size', n, 'x', n))
  }

  if (!(nrow(sigma2) == n & ncol(sigma2)==n)){
    stop(paste('sigma 1 must be a square matrix of size', n, 'x', n))
  }

  if (!(matrixcalc::is.positive.definite(sigma1) & (matrixcalc::is.positive.definite(sigma2)))){
    stop('sigma1 and sigma2 must both be positive definite')
  }

  # Transform mu1, mu2 to column vectors
  mu1 <- as.matrix(mu1, ncol=1)
  mu2 <- as.matrix(mu2, ncol=1)

  # The precision matrix (the inverse of the covariance matrix) is additive
  sigma3 <- solve(solve(sigma1) + solve(sigma2))

  # The mean of the product pdf is the precision-weighted mean of the means
  mu3 <- (sigma3 %*% (solve(sigma1) %*% mu1)) + (sigma3 %*% (solve(sigma2) %*% mu2))

  ans <- list(mu=mu3, sigma=sigma3)

  return(ans)
}
#
#------------------------------------------------------------------------------
#' Performs estimation of a bayesian dynamic linear model (DLM) of the (univariate) series y to the
#' covariates covars by Kalman filtering
#'
#' \code{dlmtfilter} fits a dynamic linear model of the time series y to the
#' covariate time series in covars.
#'
#' @param y is a vector of equally-spaced observations of the dependent series y
#' @param covars is a matrix of equally-spaced observations of covariates. Each covariate
#' is contained in a column. The number of rows of the covars matrix must be equal
#' to the length of the vector y
#' @param W.disc.factor is the discount factor applied in the estimation of W
#' @param V is blah
#'
#' @return coeffs is an (n+1) x p matrix containing the means of the posterior distribution of the coefficients
#' @return vars is a length (n+1) list of pxp covariance matrices of the posterior distributions of the coefficients.
#
#The function fits a 'local-level' model to the covariates and the dependent
#series y
#
#A local level model consists of two coupled equations
#
#y(t) = X(t) *theta(t) + e(t)      (1) theta(t+1) = theta(t) + W(t)      (2)
#
#Equation (1) is the classic linear regression model relating the observed
#series y(t) to the covariates X(t) via a linear relationship with noise.
#e(t)~N(0, v)
#
#Equation (2) is the distinctive feature of the DLM model - the coefficients
#theta(t) follow a random walk rather than being constant in the classic
#regression paradigm. Note that theta is a latent variable and is assumed never
#to be observed.
#
#

dlmfilter <- function (y, X, W.disc.factor=0.9, v=1.0, prior.dist = "N",
                    prior.mean=0, prior.var=10e+7, ng.alpha=1, ng.beta=1,
                    calc.hist.vars=TRUE, ...){

  # first perform checks on the function arguments y is the vector or matrix of
  # independent values y must be a vector or a 1-row (column) matrix

  if (!(is.vector(y))){
    if (!(is.matrix(y))){
      stop("y must be a vector of observations or 1xn matrix")
    }
    if (ncol(y) != 1){
      stop("y must be a vector of observations or 1xn matrix")
    }
  }

  if (!(is.numeric(y))){
    stop ("y must be a vector or matrix of real values")
  }

  # The covars argument must be a matrix of real values

  if (!(is.matrix(X))){
    stop ("X must be a matrix of covariates")
  }

  if (!(is.numeric(X))){
    stop ("X must be a matrix of real values")
  }

  # define n and p as the number of rows, columns of the
  # covars matrix as specified above.

  n <- nrow(X)
  p <- ncol(X)

  # y and X must be compatible

  if (length(y) != n){
    stop ("y and X must have the same number of rows")
  }

  # W.disc.factor must be a scalar between 0 and 1

  if (!(is.scalar(W.disc.factor) & W.disc.factor>0  & W.disc.factor<1)){
    stop ("W.disc.factor must be a scalar strictly between
          0 and 1")
  }

  # V must be a positive scalar

  if (!(is.scalar(v) & v > 0)){
    stop ("V must be a positive scalar")
  }

  # The prior specified must be "N" for Normal or "NG" for Normal Inverse-Gamma
  # specification

  if (!(is.element(prior.dist, c("N", "NG")))){
    stop ("The 'prior' argument must be 'N' for Normal
        or 'NG' for Normal Inverse-Gamma")
  }

  # the default value of m is 0. If m is specified it must
  # be a p-element vector or 1xp matrix

  if (!(is.numeric(prior.mean))){
     stop ("The prior.mean argument must be a real number
        or a p x p matrix")
  }

  if (is.scalar(prior.mean) | is.vector(prior.mean)){
    prior.mean<-matrix(prior.mean, nrow=1, ncol=p)
  }

  if (!(is.matrix(prior.mean) & is.numeric(prior.mean) &
        nrow(prior.mean) == 1 & ncol(prior.mean) == p)){
    stop("If prior.mean is specified it must be a vector of
         prior means on coefficients")
  }



  # If the prior variance is specified it must be a positive scalar
  # or a PSD matrix

  if (is.scalar(prior.var)) {
    if (prior.var <= 0){
      stop('The prior.sd argument must be positive scalar or a PSD matrix')
    }

    prior.var <- diag(prior.var, nrow=p)
  }

  if (is.matrix(prior.var)) {
    if (!(matrixcalc::is.positive.definite(prior.var))){
      stop ("The prior.sd argument must be a real positive number
          or a p x p positive semi-definite matrix")
    }
    if ((!(nrow(prior.var) == p) & (ncol(prior.var) == p))){
      stop ("The prior.sd argument must be a real positive num ber
            or a p x p positive semi-definite matrix")
    }
  }


  # If the columns do not have names, create numeric values
  # as the names
  if (is.null(colnames(X))){
    colnames(X) <- as.character(seq(1, ncol(X)))
  }



  # a is the list containing the prior expected values of theta ie
  # theta(t)|y(1), y(2)... y(t-1)
  a <- list()

  # m is the list containing the posterior expected value of theta, given a(t)
  # and y(t)
  m <- list()

  # f is the list containing the one-step-ahead predictive value of y ie
  # y(t+1)|y(1), y(2)... y(t). Note: This will be used to calculate the (semi?)
  # log likelihood
  f <- list()
  e <- list()

  # q is the list containing the one-step-ahead variance of the predicted
  # y value
  q <- list()

  # m is the posterior on theta, given a(t) and y(t)
  # f is the one-step-ahead predictive value of y
  # The m$0 distribution (on theta)
  # is the "posterior after 0 datapoints".
  m[[1]]<-matrix(prior.mean, nrow=p, ncol=1)

  # Create column names for the output coefficient matrix
  # that are generated from the covariate column names
  rownames(m[[1]]) <- paste("Coeff ", colnames(X))

  # C[t] is the variance matrix for the prior on theta at time t-1 We store C[t]
  # in Singular Value Decomposition form
  #
  # The SVD decomposition of a general matrix M is M = UDV' If M is n x p then U
  # is an orthogonal n x n matrix D is a diagonal n x p matrix (ie D[i,j]=0 for
  # i!=j) V is an orthogonal p x p matrix.
  #
  # M is square and PSD (as a covariance matrix is) then we can take V=U' and
  # D[i,i]>0. Therefore we have M = U S S' U'

  # We will store the U,S components of this decomposition of C

  U.C <- list()
  S.C <- list()
  U.R <- list()
  S.R <- list()

  # Compute the U,S composition of C
  svdC0 <- svd(prior.var, nu=0)
  U.C[[1]] <- svdC0$v
  S.C[[1]] <- diag(sqrt(svdC0$d), nrow=p)

  # Create a n x 1 matrix of the log-likelihoods

  llike<-vector(mode="numeric", length = n)

  # Iterate through time from t=1 to end
  for (tt in seq(length(y))){

    tX.Vinv <- matrix(X[tt, ], nrow=p) * (1/v) # calculate t(X)*V^-1

    a[tt] <- m[tt] # see equation (2) above

    Z <- sqrt((1-W.disc.factor)/W.disc.factor)
    sqrtW <- S.C[[tt]] %*% t(U.C[[tt]]) * Z


    # Compute f(t) - the predictive mean of y[[tt]]
    f[[tt]] <- X[tt,] %*% a[[tt]]     # See 2.8b p. 53



    # Compute R(t+1) in SVD form For an explanation of the below, see p.238 of
    # DLM or alternatively, equations 13-19 of Zhang, Li (1990) They refer to R
    # as P, W as Q. Note that their PHI=1 since this is a local level model

    M.1 <- svd(rbind((S.C[[tt]] %*% U.C[[tt]]), sqrtW), nu=0)

    U.R[[tt]] <- M.1$v
    S.R[[tt]] <- diag(M.1$d, nrow = length(M.1$d))

    S.Rinv <- diag(1/M.1$d, nrow = length(M.1$d))


    # Compute C(t+1) in SVD form FOr an explanation of the below, see p.238.239
    # of DLM book or equations 20-25 of Zhang, Li (1990)

    M.2 <- svd(rbind(v^(-1/2) * X[tt,] %*% U.R[[tt]] ,
               S.Rinv), nu=0)
    U.C[[tt+1]] <- U.R[[tt]] %*% M.2$v
    S.C[[tt+1]] <- diag(1/M.2$d, nrow=length(M.2$d))

    e[[tt]] <- y[[tt]] - f[[tt]]

    # Compute the predictive variance of f[[tt]]. The 'drop' function converts a
    # 1x1 matrix into a scalar
    q[[tt]] <- drop(t(X[tt,]) %*% crossprod(t(U.R[[tt]]) %*% S.R[[tt]]) %*% X[tt,]) + v

    # we calculate the loglikelihood of each of the terms
    llike[tt] <- - 0.5 * (log(2 * pi * q[[tt]]) + (e[[tt]]^2/q[[tt]]))

    # Compute the posterior on theta using the 'Kalman gain'
    # version of the filtering equation
    m[[tt+1]] <- a[[tt]] + crossprod(S.C[[tt+1]] %*% t(U.C[[tt+1]])) %*% tX.Vinv %*% e[[tt]]

  }

  vars<-list()

  # If cal.hist.vars is TRUE compute the variance matrices of the coefficients
  # from the U,S decomposition

  if (calc.hist.vars){
    for (tt in (seq(n+1))){
      vars[[tt]] <- crossprod(S.C[[tt]] %*% t(U.C[[tt]]))
      }
  }

  # Convert m from a list of vectors to a matrix and output as coeffs

  coeffs <- matrix(0, nrow=n+1, ncol=p)
  for (i in seq(n+1)){
    coeffs[i,] <- m[[i]]
  }

  # Returns coefficient mean posterior values as a list
  # and coefficient variance matrices as a list of matrices
  # Also returns the log likelihood

  ans <- list(coeffs=coeffs, vars=vars, llike=llike)

  return(ans)

}

#------------------------------------------------------------------------------
#
#' This function performs estimation of a bayesian dynamic linear model (DLM) of the (univariate) series y to the
#' covariates covars by Kalman smoothing
#'
#' \code{dlmsmooth} fits a dynamic linear model of the time series y to the
#' covariate time series in covars.
#'
#' @param y is a vector of equally-spaced observations of the dependent series y
#' @param covars is a matrix of equally-spaced observations of covariates. Each covariate
#' is contained in a column. The number of rows of the covars matrix must be equal
#' to the length of the vector y
#' @param W.disc.factor is the discount factor applied in the estimation of W
#' @param V is the observation noise variance (a positive scalar)
#' @param prior.dist can be set to
#'
#' @return coeffs a (n+1) x p matrix of posterior means on coefficients
#' @return vars a list of n+1 covariance matrices of the posterior distributions of the coefficients
#------------------------------------------------------------------------------

dlmsmooth <- function (y, X, W.disc.factor=0.9, v=1.0, prior.dist = "N",
                       prior.mean=0, prior.var=10e+7, ng.alpha=1, ng.beta=1,
                       calc.hist.vars=TRUE, ...){

    if (prior.dist!='N'){
      stop('Smoothing functionality available for normal distribution only')
    }

    n <- length(y)
    p <- ncol(X)

    filter.output <- dlmfilter(y=y, X=X, ...)

    # The posterior estimates for the forward filtering step become the prior
    # estimates for the backward step
    forward.coeff.means <- filter.output$coeffs
    forward.coeff.vars <- filter.output$vars
    forward.y.means <- filter.output$pred.means
    forward.y.vars <- filter.output$pred.vars

    posterior.coeff.means <- filter.output$coeffs[,n+1]
    posterior.coeff.vars <-filter.output$vars[[n+1]]

    fsllike <- filter.output$sllike

    y.rev <- y[n:1]
    X.rev <- X[c(n:1),]

    backward.filter.output <- dlmfilter(y=y.rev, X=X.rev, prior.mean=posterior.coeff.means,
                                        prior.var=posterior.coeff.vars, ...)

    backward.coeff.means <- backward.filter.output$coeffs[,c(eval(n+1):1)]
    backward.coeff.vars <- rev(backward.filter.output$vars)

    backward.y.means <-backward.filter.output$pred.means[c(n:1)]
    backward.y.vars <-rev(backward.filter.output$pred.vars)

    # The total log likelihood is the sum of the foward and reverse semi-loglikelihoods
    llike <- backward.filter.output$sllike + fsllike

    smoothed.dist <- list()
    smoothed.dist.means <- matrix(0, nrow=p, ncol = n+1)
    smoothed.dist.vars <- list()


    for (i in seq(n+1)){
      smoothed.dist[[i]] <- prodnorm(
        mu1=forward.coeff.means[,i], sigma1=forward.coeff.vars[[i]],
        mu2=backward.coeff.means[,i], sigma2=backward.coeff.vars[[i]])
      smoothed.dist.means[,i] <- smoothed.dist[[i]]$mu
      smoothed.dist.vars[[i]] <- smoothed.dist[[i]]$sigma
    }


    ans<-list(coeffs=smoothed.dist.means, vars=smoothed.dist.vars, llike=llike)

    return(ans)

}

#------------------------------------------------------------------------------
#
#' This function generates samples of coefficients from a DLM fitted to data y
#' and covariates covars
#'
#' \code{dlmffbs} fits a dynamic linear model of the time series y to the
#' covariate time series in covars.
#'
#' @param y is a vector of equally-spaced observations of the dependent series y
#' @param covars is an n x p matrix of equally-spaced observations of covariates.
#' Each covariate is contained in a column. The number of rows of the covars
#' matrix must be equal to the length of the vector y
#' @param W.disc.factor is the discount factor applied in the estimation of W
#' @param V is the observation noise variance (a positive scalar)
#' @param prior.dist can be set to
#'
#' @returns a list of length num.samples with each element being an n x p matrix, which is
#' a draw from the DLM
#'
#------------------------------------------------------------------------------

dlmffbs <- function (y, X, W.disc.factor = 0.9, v =1.0, prior.dist = "N",
                       prior.mean = 0, prior.var = 10e+7, num.samples = 1000){

  if (prior.dist != 'N'){
    stop('Smoothing functionality available for normal distribution only')
  }

  n <- length(y)
  p <- ncol(X)

  filter.output <- dlmfilter(y=y, X=X, y, X, W.disc.factor = W.disc.factor,
                             v = v, prior.dist = prior.dist,
                             prior.mean = prior.mean, prior.var = prior.var)

  samples.out <- list()

  for (sample.count in seq(num.samples)){

    sample <- matrix(0, nrow = n, ncol = p)

    # Sample from time n filtering distribution
    # We use the multivariate normal calculation function rmvnorm
    # from the package mvtnorm
    sample[n, ] <- mvtnorm::rmvnorm(1, mean=filter.output$coeffs[n,],
                                    sigma=filter.output$vars[[n]])

    for (tt in seq(n-1, 1)){

      # Backwards sample to time t=1 again using rmvnorm
      sample.mean <- filter.output$coeffs[tt,] +
                      W.disc.factor*(sample[tt+1,] - filter.output$coeffs[tt,])
      sample.covariance <- filter.output$vars[[tt]] * (1-W.disc.factor)

      sample[tt, ] <- mvtnorm::rmvnorm(1, mean=sample.mean, sigma=sample.covariance)
    }

    samples.out[[sample.count]] <- sample
  }

  return(samples.out)
}
