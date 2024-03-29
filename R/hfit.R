#' @include hspec.R hmoment.R hllf.R
NULL

#' Perform maximum likelihood estimation
#'
#' Generic function hfit.
#' A method for estimating the parameters of the exponential Hawkes model.
#' The reason for being constructed as the S4 method is as follows.
#' First, to represent the structure of the model as an hspec object.
#' There are numerous variations on the multivariate marked Hawkes model.
#' Second, to convey the starting point of numerical optimization.
#' The parameter values assigned to the hspec slots become initial values.
#' This function uses \code{\link[maxLik]{maxLik}} for the optimizer.
#'
#' @param object \code{\link{hspec-class}}. This object includes the parameter values
#' @param inter_arrival Inter-arrival times of events which includes inter-arrival for events that occur in all dimensions. Start with zero.
#' @param type A vector of dimensions. Distinguished by numbers, 1, 2, 3, and so on. Start with zero.
#' @param mark A vector of mark (jump) sizes. Start with zero.
#' @param N A matrix of counting processes.
#' @param Nc A matrix of counting processes weighted by mark.
#' @param lambda_component0 Initial values of lambda component. It must have the same dimensional matrix (n by n) with \code{object}.
#' @param N0 Initial values of N.
#' @param mylogLik User defined log-likelihood function. `mylogLik` function should have `object` argument consistent with \code{object}.
#' @param reduced When `TRUE`, reduced estimation performed.
#' @param constraint Constraint matrices. For more information, see \code{\link[maxLik]{maxLik}}.
#' @param method A Method for optimization. For more information, see \code{\link[maxLik]{maxLik}}.
#' @param grad A Gradient matrix for the likelihood function. For more information, see \code{\link[maxLik]{maxLik}}.
#' @param hess A Hessian matrix for the likelihood function. For more information, see \code{\link[maxLik]{maxLik}}.
#' @param verbose If `TRUE`, print the progress of the estimation.
#' @param ... Other parameters for optimization. For more information, see \code{\link[maxLik]{maxLik}}.
#'
#' @return \code{\link{maxLik}} object
#'
#' @docType methods
#' @rdname hfit
#' @export
#'
#' @seealso \code{\link{hspec-class}}, \code{\link{hsim,hspec-method}}
#'
#' @examples
#'
#' # example 1
#' mu <- c(0.1, 0.1)
#' alpha <- matrix(c(0.2, 0.1, 0.1, 0.2), nrow=2, byrow=TRUE)
#' beta <- matrix(c(0.9, 0.9, 0.9, 0.9), nrow=2, byrow=TRUE)
#' h <- new("hspec", mu=mu, alpha=alpha, beta=beta)
#' res <- hsim(h, size=100)
#' summary(hfit(h, inter_arrival=res$inter_arrival, type=res$type))
#'
#'
#' # example 2
#' \donttest{
#' mu <- matrix(c(0.08, 0.08, 0.05, 0.05), nrow = 4)
#' alpha <- function(param = c(alpha11 = 0, alpha12 = 0.4, alpha33 = 0.5, alpha34 = 0.3)){
#'   matrix(c(param["alpha11"], param["alpha12"], 0, 0,
#'            param["alpha12"], param["alpha11"], 0, 0,
#'            0, 0, param["alpha33"], param["alpha34"],
#'            0, 0, param["alpha34"], param["alpha33"]), nrow = 4, byrow = TRUE)
#' }
#' beta <- matrix(c(rep(0.6, 8), rep(1.2, 8)), nrow = 4, byrow = TRUE)
#'
#' impact <- function(param = c(alpha1n=0, alpha1w=0.2, alpha2n=0.001, alpha2w=0.1),
#'                    n=n, N=N, ...){
#'
#'   Psi <- matrix(c(0, 0, param['alpha1w'], param['alpha1n'],
#'                   0, 0, param['alpha1n'], param['alpha1w'],
#'                   param['alpha2w'], param['alpha2n'], 0, 0,
#'                   param['alpha2n'], param['alpha2w'], 0, 0), nrow=4, byrow=TRUE)
#'
#'   ind <- N[,"N1"][n] - N[,"N2"][n] > N[,"N3"][n] - N[,"N4"][n] + 0.5
#'
#'   km <- matrix(c(!ind, !ind, !ind, !ind,
#'                  ind, ind, ind, ind,
#'                  ind, ind, ind, ind,
#'                  !ind, !ind, !ind, !ind), nrow = 4, byrow = TRUE)
#'
#'   km * Psi
#' }
#' h <- new("hspec",
#'          mu = mu, alpha = alpha, beta = beta, impact = impact)
#' hr <- hsim(h, size=100)
#' plot(hr$arrival, hr$N[,'N1'] - hr$N[,'N2'], type='s')
#' lines(hr$N[,'N3'] - hr$N[,'N4'], type='s', col='red')
#' fit <- hfit(h, hr$inter_arrival, hr$type)
#' summary(fit)
#' }
#'
#' # example 3
#' \donttest{
#' mu <- c(0.15, 0.15)
#' alpha <- matrix(c(0.75, 0.6, 0.6, 0.75), nrow=2, byrow=TRUE)
#' beta <- matrix(c(2.6, 2.6, 2.6, 2.6), nrow=2, byrow=TRUE)
#' rmark <- function(param = c(p=0.65), ...){
#'   rgeom(1, p=param[1]) + 1
#' }
#' impact <- function(param = c(eta1=0.2), alpha, n, mark, ...){
#'   ma <- matrix(rep(mark[n]-1, 4), nrow = 2)
#'   alpha * ma * matrix( rep(param["eta1"], 4), nrow=2)
#' }
#' h1 <- new("hspec", mu=mu, alpha=alpha, beta=beta,
#'           rmark = rmark,
#'           impact=impact)
#' res <- hsim(h1, size=100, lambda_component0 = matrix(rep(0.1,4), nrow=2))
#'
#' fit <- hfit(h1,
#'             inter_arrival = res$inter_arrival,
#'             type = res$type,
#'             mark = res$mark,
#'             lambda_component0 = matrix(rep(0.1,4), nrow=2))
#' summary(fit)
#' }
#'# For more information, please see vignettes.
setGeneric("hfit", function(object, inter_arrival = NULL,
                            type = NULL, mark = NULL,
                            N = NULL, Nc = NULL,
                            lambda_component0 = NULL, N0 = NULL,
                            mylogLik = NULL,
                            reduced = TRUE,
                            grad = NULL, hess = NULL, constraint = NULL, method = "BFGS",
                            verbose = FALSE, ...) standardGeneric("hfit"))

#' @rdname hfit
setMethod(
  f="hfit",
  signature(object="hspec"),
  function(object, inter_arrival = NULL,
           type = NULL, mark = NULL,
           N = NULL, Nc = NULL,
           lambda_component0 = NULL, N0 = NULL,
           mylogLik = NULL,
           reduced = TRUE,
           grad = NULL, hess = NULL, constraint = NULL, method = "BFGS",
           verbose = FALSE, ...){

    additional_argument <- list(...)
    if ("lambda0" %in% names(additional_argument)) {

      warning("lambda0 is deprecated; instead use lambda_component0")

      lambda_component0 <- additional_argument[["lambda0"]]

    }

    if(is.null(lambda_component0)){
      message("The initial values for intensity processes are not provided. Internally determined initial values are used for estimation.\n")
    }

    # parameter setting
    plist <- setting(object)
    mu <- plist$mu
    alpha <- plist$alpha
    beta <- plist$beta
    eta <- plist$eta
    impact <- plist$impact
    rmark <- plist$rmark
    dmark <- plist$dmark
    dimens <- plist$dimens

    # parameter setting
    if(is.function(mu)){
      # When mu is provided by a function,
      # the first formal arguments is a parameter vector.
      # For, example, formals(object@mu)[[1]] is equal to c("mu1"=0.1, "mu2"=0.2)
      pr_mus <- eval(formals(object@mu)[[1]])
    }else{
      # When mu is a matrix, select only unique values as parameters (if reduced is TRUE)
      pr_mus <- as.param(object@mu, "mu", reduced)
    }

    pr_alphas <- as.param(object@alpha, "alpha", reduced)
    pr_betas <- as.param(object@beta, "beta", reduced)

    if(!is.null(eta)){
      pr_etas <- as.param(object@eta, "eta", reduced)
    } else {
      pr_etas <- NULL
    }

    if(!is.null(impact)){
      #get default parameter values
      pr_impact <- eval(formals(object@impact)[[1]])
    } else {
      pr_impact <- NULL
    }

    if(!is.null(dmark)){
      #get default parameter values
      pr_mark <- eval(formals(object@dmark)[[1]])
    } else {
      pr_mark <- NULL
    }

    names_mu <- names(pr_mus)
    names_alpha <- names(pr_alphas)
    names_beta <- names(pr_betas)
    names_eta <- names(pr_etas)
    names_impact <- names(pr_impact)
    names_mark <- names(pr_mark)

    # same name same parameter rule
    input_parameters <- c(pr_mus, pr_alphas, pr_betas, pr_etas, pr_impact, pr_mark)
    starting_point <- input_parameters[unique(names(input_parameters))]

    # loglikelihood function for maxLik

    llh_function <- function(param){

      # These values may be generated by maxLik repeatedly.
      # name based query

      pr_mus <- param[names_mu]
      pr_alphas <- param[names_alpha]
      pr_betas <- param[names_beta]
      pr_etas <- param[names_eta]

      pr_impact <- param[names_impact]
      pr_mark <- param[names_mark]

      # Reconstruct hspec
      # Convert to matrix
      if (is.function(object@mu)){
        mu0 <- hijack(object@mu, param = pr_mus)
      } else{
        if(reduced){
          mu0 <- matrix(pr_mus[look_up_mtrx(mu, "mu")], nrow=dimens)
        } else {
          mu0 <- matrix(pr_mus, nrow=dimens)
        }
      }
      if (is.function(object@alpha)){
        alpha0 <- hijack(object@alpha, param = pr_alphas)
      } else{
        if(reduced){
          alpha0 <- matrix(pr_alphas[look_up_mtrx(alpha, "alpha")], nrow=dimens)
        } else {
          alpha0 <- matrix(pr_alphas, nrow=dimens)
        }
      }
      if (is.function(object@beta)){

        beta0 <- hijack(object@beta, param = pr_betas)

      } else{
        if(reduced){
          beta0 <- matrix(pr_betas[look_up_mtrx(beta, "beta")], nrow=dimens)
        } else {
          beta0 <- matrix(pr_betas, nrow=dimens)
        }
      }

      if (!is.null(object@eta)){
        if(is.function(object@eta)) {
          eta0 <- hijack(object@eta, param = pr_etas)
        } else{
          if(reduced){
            eta0 <- matrix(pr_etas[look_up_mtrx(eta, "eta")], nrow=dimens)
          } else {
            eta0 <- matrix(pr_etas, nrow=dimens)
          }
        }
      } else {
        eta0 <- NULL
      }

      if (!is.null(object@dmark)){
        dmark0 <- hijack(object@dmark, param = pr_mark)
      } else {
        dmark0 <- NULL
      }

      #object@impact is user defined impact function
      if (!is.null(object@impact)){
        impact0 <- hijack(object@impact, param = pr_impact)
      } else {
        impact0 <- NULL
      }

      hspec0 <- methods::new("hspec", mu = mu0, alpha = alpha0, beta = beta0, eta = eta0,
                             impact = impact0, dmark = dmark0,
                             rmark = object@rmark, type_col_map = object@type_col_map)


      #this_flag_represents_binding_env_is_hfit <- TRUE


      if (is.null(mylogLik)){

        logl <- logLik(hspec0, inter_arrival = inter_arrival, type = type,
                       mark = mark, N = N, Nc = Nc, N0 = N0, lambda_component0 = lambda_component0,
                       showWarning = FALSE)

      } else {
        # arguments names for mylogLik
        args_needed <- names(formals(mylogLik))
        if ( !("object" %in% args_needed)) {
          stop('mylogLik needs \'object\' arguments with hspec')
        }
        args_logLik <- vector("list", length(args_needed))
        for (i in seq_along(args_needed)){
          # find necessary object for mylogLik arguments
          if(args_needed[i] == "object") args_logLik[[i]] <- hspec0
          else args_logLik[[i]] <- get(args_needed[i])
        }
        names(args_logLik) <- args_needed
        logl <- do.call(mylogLik, args = args_logLik)
      }

      if(verbose){
        cat("Parameters : ", param, "\n")
        cat("Log likelihood : ", logl, "\n")
      }

      logl

    }


    #llh_function(starting_point)
    maxLik::maxLik(logLik=llh_function,
                   start=starting_point, grad, hess, method = method, constraint = constraint, ...)

  }
)

