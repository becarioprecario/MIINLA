#' @name inla.rgeneric.milm.model
#' @rdname milm
#'
#' @title Multiple imputation using a linear regression.
#'
#' @description Multiple imputation using a linear regression model. The response may have missing values to be imputed, but the covariates must be fully observed.
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
#'   \item \emph{XX} Matrix of covariates to be used in the multiple imputation linear model. Must be fully observed.
#'
#'   \item \emph{n} Number of observations.
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
#' library(INLA)
#' library(MIINLA)
#' library(mice)
#' data(nhanes)
#'
#' nhanes$age1 <- as.numeric(nhanes$age == 1)
#' nhanes$age2 <- as.numeric(nhanes$age == 2)
#' nhanes$age3 <- as.numeric(nhanes$age == 3)
#' 
#' # Standard model (if bmi == NA then the term is not included in l.pred.)
#' m0 <- inla(chl ~ bmi, data = nhanes, control.predictor = list(compute = TRUE))
#' summary(m0)

#' # Imputation model on 'bmi'
#' r.imp <- inla(bmi ~ age2 + age3, data = nhanes, control.predictor = list(compute = TRUE))
#' summary(r.imp)
#'
#' # Define imputation model
#'  model <- inla.rgeneric.define(inla.rgeneric.milm.model, debug = TRUE,
#'    x = nhanes$bmi,
#'    XX = cbind(1, nhanes$age2, nhanes$age3),
#'    n = nrow(nhanes),
#'    idx.na = which(is.na(nhanes$bmi)))
#'
#' nhanes$idx <- rep(NA, nrow(nhanes))
#' nhanes$idx2 <- 1:nrow(nhanes)
#' formula = chl ~ 1 + f(idx, model = model) + f(idx2, copy = "idx", fixed = FALSE,
#'    hyper = list(beta = list(prior = "normal", param = c(0, 0.001))))
#' 
#' #Joint model: main model + imputation of missing covariates
#' r = inla(formula, data = nhanes[, c("chl", "idx", "idx2")],
#'   family = "gaussian",
#'   verbose = TRUE,
#'   control.family = list(hyper = list(prec = list(param = c(0.01, 0.01)))),
#'   control.fixed = list(prec.intercept = 0.001))
#' summary(r)
#'
#' r.imp$summary.fitted.values[, "mean"]
#' r$summary.random$idx[, "mean"]
#'
#' @export

'inla.rgeneric.milm.model' = function(cmd = c("graph", "Q", "mu", "initial",
  "log.norm.const", "log.prior", "quit"), theta = NULL)
{
  interpret.theta = function()
  {
  ## internal helper-function to map the parameters from the internal-scale to the
  ## user-scale
  p <- ncol(XX)
  return (list(beta = theta[as.integer(1:p)],
    prec = exp(theta[as.integer(p + 1)])))
}

graph = function()
{
  G = Matrix::Diagonal(n, 1)
  return (G) 
}

Q = function()
{
    ## returns the precision matrix for given parameters
    param = interpret.theta()
    Q.diag <- rep(1, n)
    Q.diag[-idx.na] <- 10^10
    Q.diag[idx.na] <- param$prec
    Q <- Matrix::Diagonal(n, Q.diag)
  return (Q)
}

mu = function() {
  param = interpret.theta()
  mu.x <- x
  
  mu.x [idx.na] <- XX [idx.na, ] %*% matrix(param$beta, ncol = 1)
  return(mu.x)
}

log.norm.const = function() {
  ## return the log(normalising constant) for the model
  param = interpret.theta()
  return (numeric(0))
}
    log.prior = function()
    {
      ## return the log-prior for the hyperparameters. the ’+theta[1L]’ is the log(Jacobian)
      ## for having a gamma prior on the precision and convert it into the prior for the
      ## log(precision).
      param <- interpret.theta()
      #val <- dnorm(param$alpha, 25.161, 1 / 100000, log = TRUE) 
      #val <- val + dnorm(param$beta, 3.203, 1 / 100000, log = TRUE) 
      val <- sum(dnorm(param$beta, 0, 1000, log = TRUE))
      val <- val + (dgamma(param$prec, shape = 0.01, rate = 0.01, log = TRUE) +
        theta[1L])
      val <- val + 0.5 * n *log(param$prec) -0.5 * param$prec * 
        sum((x[-idx.na] - (XX [-idx.na, ] %*% matrix(param$beta, ncol = 1))[, 1])^2)
        return (val)
    }

    initial = function()
    {
        ## return initial values
        ntheta = ncol(XX) + 1
        return (rep(0, ntheta))
}

    quit = function()
    {
        return (invisible())
    }
    val = do.call(match.arg(cmd), args = list())
    return (val)
}

