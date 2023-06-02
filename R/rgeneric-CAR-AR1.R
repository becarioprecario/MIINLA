#' @name inla.rgeneric.micarar1.model
#' @rdname micarar1
#'
#' @title Multiple imputation with INLA using a spatio-temporal model
#'
#' @description Multiple imputation using a spatial CAR regression model. The response may have missing values to be imputed, but the covariates must be fully observed.
#'
#' @param cmd Arguments used by latent effects defined using the 'rgeneric' latent effect.
#' @param theta Vector of hyperparameters.
#'
#' @return This is used internally by the 'INLA::inla()'.
#'
#' @details This function used is used to define a latent effect that is a linear term on some covariates with missing observations. However, multiple imputation is performed on the missing values of the covariates internally using another linear model.  For this reason, this package requires the following arguments
#' when defining the latent effect:
#' \itemize{
#'   \item \emph{x} Vector of covariates (with missing observations).
#'
#'   \item \emph{W} SCALED (i.e., divided by its maximum eigenvalue) adjacency SPARSE matrix for spatial effect.
#'
#'   \item \emph{n} Number of areas.
#'
#'   \item \emph{p} Number of temporal units.
#'
#'   \item \emph{idx.na} Index with the positions of the missing observations.
#'
#'}
#'
#' This model is defined using the 'f()' function and an index of \code{NA}'s
#' is set so that imputation is done but the covariate not included in the
#' actual model. Then, this latent effect is 'copied', which makes the
#' covariates (observed and imputed values) into a linear term in the liner
#' predictor. See the example.
#'
#' @examples
#'
#' library(spdep)
#' library(rgdal)
#' library(sp)
#' library(INLA)
#' 
#' # Load data
#' nc.sids <- readOGR(system.file("shapes/sids.shp", package="spData")[1])
#' proj4string(nc.sids) <- CRS("+proj=longlat +ellps=clrk66")

#' # Compute covariate and expected counts
#' nc.sids$NWPROP74 <- (nc.sids$NWBIR74 / nc.sids$BIR74)
#' nc.sids$EXP74 <- nc.sids$BIR74 * sum(nc.sids$SID74) / sum(nc.sids$BIR74)
#' 
#' # Create covariate with missing observations
#' nc.sids$NWPROP74M <- nc.sids$NWPROP74
#' idx.na <- sample(1:100, 20) #1:10 #seq(1, 100, by = 10)
#' nc.sids$NWPROP74M [idx.na] <- NA
#'
#' # Standard model
#' m0 <- inla(SID74 ~ NWPROP74, family = "poisson", data = as.data.frame(nc.sids),
#'  E = EXP74)
#' summary(m0)
#' 
#' # Imputation model
#' adj <- poly2nb(nc.sids)
#' W <- as(nb2mat(adj, style = "B"), "Matrix")
#' W.scaled <- W / max(eigen(W)$values)
#'
#' nc.sids$idx <- 1:nrow(nc.sids)
#' r.imp <- inla(NWPROP74M ~ 0 + f(idx, model = "generic1", Cmatrix = W.scaled),
#'   data = as.data.frame(nc.sids),
#'   control.predictor = list(compute = TRUE))
#'
#' plot(r.imp$summary.random$idx[idx.na, "mean"], nc.sids$NWPROP74[idx.na])
#' abline(0, 1)
#' 
#' 
#' model = inla.rgeneric.define(inla.rgeneric.micar.model, debug = TRUE,
#'  n = nrow(nc.sids),
#'  x = nc.sids$NWPROP74M,
#'  idx.na = which(is.na(nc.sids$NWPROP74M)),
#'  W = W.scaled)
#' 
#' nc.sids$idxNA <- NA
#' formula = SID74 ~ 1 + f(idxNA, model = model) +
#'   f(idx, copy = "idxNA", fixed = FALSE,
#'    hyper = list(beta = list(prior = "normal", param = c(0, 0.001))))
#'
#' r = inla(formula, data = as.data.frame(nc.sids),
#'  family = "poisson", E = EXP74,
#'  verbose = TRUE)
#' summary(r)
#'
#' r.imp$summary.fitted.values[idx.na, "mean"]
#' r$summary.random$idx[idx.na, "mean"]
#' nc.sids$NWPROP74[idx.na]
#'
#' @export

'inla.rgeneric.micarar1.model' = function(cmd = c("graph", "Q", "mu", "initial",
  "log.norm.const", "log.prior", "quit"), theta = NULL)
{
  interpret.theta <- function()
  {
    print(theta)

    ## internal helper-function to map the parameters from the internal-scale to the
    ## user-scale
    return (list(prec = exp(theta[1L]),#theta1 = log(tau)
      rho1 = 1 / (1 + exp(-theta[2L])), #theta2 = logit(rho1) SPATIAL
      alpha = theta[3L], #theta3 = alpha
      rho2 = 1 / (1 + exp(-theta[4L])))) #theta4 = logit(rho2) TEMPORAL
  }

  #Structure of precision of CAR
  prec_CAR <- function(rho) {
    return(Diagonal(n, x = 1) - rho * W)
  }

  #Structure of precision of AR1
  # Modified from https://haakonbakkagit.github.io/btopic120.html
  prec_AR1 <- function(rho) {
    PREC <- Matrix(0, p, p)
    diag(PREC) <- 1 + rho^2
    for (i in 1:(p - 1)) {
      PREC[i, i + 1] = -rho
      PREC[i + 1, i] = -rho
    }

    PREC[1, 1] <- 1
    PREC[p, p] <- 1

    return(PREC)
  }

  # Structure of spatio-temporal precision matrix
  # Spatial index moves faster.
  prec_CAR_AR1 <- function(rho1, rho2) {
    return(Matrix::kronecker(prec_AR1(rho2), prec_CAR(rho1)))
  }
 
  graph <- function()
  {
    G <- Matrix::Diagonal(n * p, 1)
    # rho1 = rho2 = 0.5 because actual values are not needed
    Qst <- prec_CAR_AR1(0.5, 0.5)
    G[idx.na, idx.na] <- Qst[idx.na, idx.na]
    return (G) 
  }

  Q <- function()
  {
    ## returns the precision matrix for given parameters
    param = interpret.theta()
    # Matrix of full s-t effect
    Qst <- prec_CAR_AR1(param$rho1, param$rho2)

    # Precion matrix for imputation
    Q.diag <- rep(0, n * p)
    Q.diag[-idx.na] <- 10^10
    Q <- Matrix::Diagonal(n * p , Q.diag)
    Q[idx.na, idx.na] <- param$prec * Qst[idx.na, idx.na]

    return (Q)
  }

  mu <- function() {
    param <- interpret.theta()
    Qst <- prec_CAR_AR1(param$rho1, param$rho2)
    
    # Matrix of full s-t effect
    print(determinant(Qst))

    #Mean is the observed values
    mu.x <- x
    # Fill NA's according to the imputation model
    mu.x [idx.na] <- param$alpha -
      (Matrix::solve(Qst[idx.na, idx.na], Qst[idx.na, -idx.na]) %*% 
        matrix(x[-idx.na] - param$alpha, ncol = 1))[, 1]
 
    return(mu.x)
  }

log.norm.const <- function() {
  ## return the log(normalising constant) for the model
  return (numeric(0))
}
    log.prior <- function()
    {
      ## return the log-prior for the hyperparameters. the ’+theta[1L]’ is the log(Jacobian)
      ## for having a gamma prior on the precision and convert it into the prior for the
      ## log(precision).
      param <- interpret.theta()
      # Matrix of full s-t effect
      Qst <- prec_CAR_AR1(param$rho1, param$rho2)

      val <- (dgamma(param$prec, shape = 0.01, rate = 0.01, log = TRUE) +
        theta[1L])
      val <- val + dnorm(theta[2L], 0, sqrt(10), log = TRUE)
      val <- val + dnorm(param$alpha, 0, sqrt(1000), log = TRUE)
      val <- val + dnorm(theta[4L], 0, sqrt(10), log = TRUE)

  #Number of non-NA values
  k <- n * p - length(idx.na)
  Q <- param$prec * Qst[-idx.na, -idx.na] 

  # Values of x minus the estimate of the mean 'alpha'
  x.mat <- matrix(x[-idx.na] - param$alpha, ncol = 1)

  log_detQ <- Matrix::determinant(Q)
  if(log_detQ$sign < 0) stop("Negative determinant of Q")

  #val <- val - k * 0.5 * log(2 * pi) + 0.5 * log(Matrix::det(Q)) - 
  val <- val - k * 0.5 * log(2 * pi) + 0.5 * log_detQ$modulus - 
    0.5 * (t(x.mat) %*% Q %*% x.mat)
      return (as.numeric(val))
    }

    initial <- function()
    {
        ## return initial values
        ntheta = 4
        return (rep(0, ntheta))
    }

    quit <- function()
    {
        return (invisible())
    }

    # FIX for rgeneric to work on R >= 4
    # Provided by E. T. Krainski
    if (as.integer(R.version$major) > 3) {
      if (!length(theta))
        theta <- initial()
    } else {
      if (is.null(theta)) {
        theta <- initial()
      }
    }

    val <- do.call(match.arg(cmd), args = list())
    return (val)
}

