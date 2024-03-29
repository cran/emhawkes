#' @include hspec.R hmoment.R utilities.R
NULL

#' Simulate multivariate Hawkes process with exponential kernel.
#'
#' @description
#' The method simulate multivariate Hawkes processes.
#' The object \code{\link{hspec-class}} contains the parameter values such as \code{mu}, \code{alpha}, \code{beta}.
#' The mark (jump) structure may or may not be included.
#' It returns an object of class \code{\link{hreal}} which contains \code{inter_arrival}, \code{arrival},
#' \code{type}, \code{mark}, \code{N}, \code{Nc}, \code{lambda}, \code{lambda_component}, \code{rambda}, \code{rambda_component}.
#'
#'
#' @param object \code{\link{hspec-class}}. S4 object that specifies the parameter values.
#' @param lambda_component0 Starting values of lambda component. numeric or matrix.
#' @param N0 Starting values of N with default value 0.
#' @param size Number of observations.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return \code{\link{hreal}} S3-object, summary of the Hawkes process realization.
#'
#' @rdname hsim
#' @docType methods
#' @examples
#'
#' # example 1
#'
#' mu <- 1; alpha <- 1; beta <- 2
#' h <- new("hspec", mu=mu, alpha=alpha, beta=beta)
#' hsim(h, size=100)
#'
#'
#' # example 2
#' mu <- matrix(c(0.1, 0.1), nrow=2)
#' alpha <- matrix(c(0.2, 0.1, 0.1, 0.2), nrow=2, byrow=TRUE)
#' beta <- matrix(c(0.9, 0.9, 0.9, 0.9), nrow=2, byrow=TRUE)
#' h <- new("hspec", mu=mu, alpha=alpha, beta=beta)
#' res <- hsim(h, size=100)
#' print(res)
#'
#' @export
setGeneric("hsim", function(object, size = 100, lambda_component0 = NULL,  N0 = NULL, ...)
  standardGeneric("hsim"))

#' @rdname hsim
setMethod(
  f="hsim",
  signature(object = "hspec"),
  definition = function(object, size = 100, lambda_component0 = NULL,  N0 = NULL, ...){

    additional_argument <- list(...)
    if ("lambda0" %in% names(additional_argument)) {

      warning("lambda0 is deprecated; instead use lambda_component0.")

      lambda_component0 <- additional_argument[["lambda0"]]

    }

    # parameter setting
    plist <- setting(object)
    mu <- plist$mu
    alpha <- plist$alpha
    beta <- plist$beta
    eta <- plist$eta
    impact <- plist$impact
    rmark <- plist$rmark
    dimens <- plist$dimens

    type_col_map <- object@type_col_map

    N <- matrix(numeric(length = dimens * size), ncol = dimens)
    Nc  <- matrix(numeric(length = dimens * size), ncol = dimens)
    colnames(Nc)  <- paste0("Nc", 1:dimens)
    colnames(N) <- paste0("N", 1:dimens)

    type <- numeric(length = size)
    inter_arrival <- numeric(length = size)
    mark <- c(0, rep(1, size-1))

    if (!is.null(N0)){
      N[1,] <- N0
      Nc[1,] <- N0
    }

    if (is.function(mu)){
      mu0 <- mu(n = 1, mark = mark, type = type, inter_arrival = inter_arrival,
                 N = N, Nc = Nc,
                 alpha = alpha, beta = beta)
    } else {
      mu0 <- mu
    }

    # default lambda_component0
    if(!is.null(type_col_map)){
      if(length(type_col_map) > 0 & is.null(lambda_component0)){
        stop("In this model, please provide the initial values of lambda_component0.")
      }
    }


    if(!is.null(lambda_component0)){

      # If the dimensions of model and lambda_component0 do not match, lambda_component0 will be adjusted
      if (dimens * ncol(beta) > length(lambda_component0)){
        warning("The size of lambda_component0 does not match to the dimension of the model and is adjusted. \n
                lambda_component0 is now :")
        lambda_component0 <- rep(lambda_component0, dimens * ncol(beta))
      }
      if (dimens * ncol(beta) < length(lambda_component0)){
        warning("The size of lambda_component0 does not match to the dimension of the model and is adjusted.\n
                lambda_component0 is now :")
        lambda_component0 <- lambda_component0[1:dimens * ncol(beta)]
      }

      lambda_component0 <- as.matrix(lambda_component0, nrow = dimens)

    } else {
      # default lambda_component0
      message("The initial values for intensity processes are not provided. Internally determined initial values are used for simulation.\n")
      lambda_component0 <- get_lambda0(object, mark = mark, type = type, inter_arrival = inter_arrival,
                             N = N, Nc = Nc,
                             mu, alpha = alpha, beta = beta)
    }

    # Preallocation for lambdas and Ns and set initial values for lambdas
    lambda_component <- matrix(sapply(lambda_component0, c, numeric(length = size - 1)), ncol = dimens * ncol(beta))
    rowSums_lambda_component0 <- rowSums(matrix(lambda_component0, nrow=dimens))

    lambda <- matrix(sapply(mu0 + rowSums_lambda_component0, c, numeric(length = size - 1)), ncol = dimens)

    rambda <- lambda
    rambda_component <- lambda_component

    # Set column names
    colnames(lambda) <- paste0("lambda", 1:dimens)
    indxM <- matrix(rep(1:ncol(beta), dimens), byrow = TRUE, nrow = dimens)
    colnames(lambda_component) <- as.vector(t(outer(as.vector(outer("lambda", 1:dimens, FUN = paste0)),
                                                    1:ncol(beta), FUN=paste0)))

    colnames(rambda) <- paste0("rambda", 1:dimens)
    colnames(rambda_component) <- as.vector(t(outer(as.vector(outer("rambda", 1:dimens, FUN = paste0)),
                                                    1:ncol(beta), FUN=paste0)))


    if (is.function(mu)){
      rmu_n <- mu(n = 2, mark = mark, type = type, inter_arrival = inter_arrival,
                 N = N, Nc = Nc, lambda = lambda, lambda_component = lambda_component,
                 alpha = alpha, beta = beta)
    } else{
      rmu_n <- mu
    }

    #current_rambda_component <- lambda_component0
    for (n in 2:size) {

      # only piecewise constant mu is available, mu is left continuous.
      # rmu_n is previously defined right continuous version.
      mu_n <- rmu_n


      ### Determine the next arrival, type and mark #######################################################
      res <- rarrival(n = n, mark = mark, type = type, inter_arrival = inter_arrival,
               N = N, Nc = Nc, lambda = lambda, lambda_component = lambda_component,
               rambda = rambda, rambda_component = rambda_component,
               mu = mu_n, alpha = alpha, beta = beta, dimens = dimens, type_col_map = object@type_col_map)
      inter_arrival[n] <- res["inter_arrival"]
      type[n] <- res["type"]


      N[n, ] <- N[n-1, ]
      N[n, type[n]] <- N[n-1, type[n]] + 1

      # generate a mark for Hawkes
      if( !is.null(rmark) ){
        # mark may depends on other variables
        mark[n] <- rmark(n = n, Nc = Nc, N = N,
                         lambda = lambda, lambda_component = lambda_component,
                         type = type)
      }

      Nc[n, ] <- Nc[n-1, ]
      Nc[n, type[n]] <- Nc[n-1, type[n]] + mark[n]
      ######################################################################################################


      # lambda decayed due to time, impact due to mark is not added yet
      decayed <- exp(-beta * inter_arrival[n])
      #decayed_lambda <- current_rambda_component * decayed
      decayed_lambda <- lambda_component_n <- matrix(rambda_component[n-1,], nrow = dimens, byrow = TRUE) * decayed

      # update lambda
      lambda_component[n, ] <- t(decayed_lambda)
      lambda[n, ] <- mu_n + rowSums(decayed_lambda)


      # impact by alpha
      impact_alpha <- matrix(rep(0, dimens * ncol(beta)), nrow = dimens)

      if( is.null(object@type_col_map)){
        types <- type[n]
      } else{
        types <- object@type_col_map[[type[n]]]
      }

      impact_alpha[ , types] <- alpha[ , types]

      new_lambda <- decayed_lambda + impact_alpha

      # additional impact by eta
      if(!is.null(eta)){

        impact_eta <- matrix(rep(0, dimens * ncol(beta)), nrow = dimens)
        impact_eta[ , types] <- eta[ , types] * (mark[n] - 1)

        new_lambda <- new_lambda + impact_eta
      }


      # new_lambda = [[lambda11, lambda12, ...], [lambda21, lambda22, ...], ...]
      if(!is.null(impact)){
        # impact by mark
        impact_mark <- matrix(rep(0, dimens * ncol(beta)), nrow = dimens)

        impact_res <- impact(n = n, mark = mark, type = type, inter_arrival = inter_arrival,
                             N = N, Nc = Nc, lambda = lambda, lambda_component = lambda_component,
                             lambda_component_n = lambda_component_n,
                             mu = mu, alpha = alpha, beta = beta)

        impact_mark[ , types] <- impact_res[ , types]

        new_lambda <- new_lambda + impact_mark

      }

      # update rambda
      # rambda_component = {"rambda11", "rambda12", ..., "rambda21", "rambda22", ...}
      rambda_component[n, ] <- t(new_lambda)

      if (is.function(mu)){
        rmu_n <- mu(n = n + 1, mark = mark, type = type, inter_arrival = inter_arrival,
                    N = N, Nc = Nc, lambda = lambda, lambda_component = lambda_component,
                    alpha = alpha, beta = beta)
      } else{
        rmu_n <- mu
      }
      rambda[n, ] <- rmu_n + rowSums(new_lambda)

      #current_rambda_component <- new_lambda
    }

    realization <- list(object, inter_arrival, cumsum(inter_arrival), type, mark,
                        N, Nc, lambda, lambda_component, rambda, rambda_component)
    names(realization) <- c("hspec", "inter_arrival", "arrival", "type", "mark",
                            "N", "Nc", "lambda", "lambda_component", "rambda", "rambda_component")
    class(realization) <- c("hreal")

    realization
  }
)


