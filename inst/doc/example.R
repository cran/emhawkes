## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
knitr::opts_chunk$set(eval = requireNamespace("maxLik", quietly = TRUE))


## -----------------------------------------------------------------------------
library(emhawkes)

## -----------------------------------------------------------------------------
mu1 <- 0.3; alpha1 <- 1.2; beta1 <- 1.5
hspec1 <- new("hspec", mu = mu1, alpha = alpha1, beta = beta1)
show(hspec1)

## ----warning=FALSE------------------------------------------------------------
set.seed(1107)
res1 <- hsim(hspec1, size = 1000)
summary(res1)

## -----------------------------------------------------------------------------
# The first and third columns are equal.
cbind(res1$lambda[1:5], res1$lambda_component[1:5], mu1 + res1$lambda_component[1:5])

## -----------------------------------------------------------------------------
# The second and third columns are equal except in the initial row.
cbind(res1$lambda[1:5], res1$rambda[1:5], res1$lambda[1:5] + alpha1)

## -----------------------------------------------------------------------------
# By definition, the following two are equal:
res1$lambda[2:6]
mu1 + (res1$rambda[1:5] - mu1) * exp(-beta1 * res1$inter_arrival[2:6])

## ----warning=FALSE------------------------------------------------------------
logLik(hspec1, inter_arrival = res1$inter_arrival)

## ----warning=FALSE------------------------------------------------------------
# Initial parameter values for numerical optimization.
mu0 <- 0.5; alpha0 <- 1.0; beta0 <- 1.8
hspec0 <- new("hspec", mu = mu0, alpha = alpha0, beta = beta0)
# Supply the initial values through hspec0.
mle <- hfit(hspec0, inter_arrival = res1$inter_arrival)
summary(mle)

## ----fit-reconstruction, message=FALSE----------------------------------------
estimated <- fitted_hspec(mle)
initial <- fit_initial_state(mle)
inferred <- do.call(logLik, c(list(object = estimated,
  inter_arrival = res1$inter_arrival, infer = TRUE), initial))
c(fitted = as.numeric(logLik(mle)), restored = inferred$loglikelihood)
head(inferred$lambda)
head(inferred$res_process1)

## ----residual-check, fig.height=4, fig.width=4--------------------------------
p <- ppoints(100)
q <- quantile(inferred$res_process1, probs = p)
plot(qexp(p), q, xlab = "Theoretical quantiles", ylab = "Sample quantiles")
abline(0, 1, col = "blue", lty = 2)

## -----------------------------------------------------------------------------
mu2 <- matrix(c(0.2), nrow = 2)
alpha2 <- matrix(c(0.5, 0.9, 0.9, 0.5), nrow = 2, byrow = TRUE)
beta2 <- matrix(c(2.25, 2.25, 2.25, 2.25), nrow = 2, byrow = TRUE)
hspec2 <- new("hspec", mu=mu2, alpha=alpha2, beta=beta2)
print(hspec2)

## -----------------------------------------------------------------------------
set.seed(1107)
res2 <- hsim(hspec2,  size=1000)
summary(res2)

## -----------------------------------------------------------------------------
# A bivariate model has two event types: 1 and 2.
res2$type[1:10]

## -----------------------------------------------------------------------------
res2$N[1:3, ]

## -----------------------------------------------------------------------------
res2$lambda[1:3, ]

## -----------------------------------------------------------------------------
res2$lambda_component[1:3, ]

## -----------------------------------------------------------------------------
mu2[1] + rowSums(res2$lambda_component[1:5, c("lambda11", "lambda12")])
res2$lambda[1:5, "lambda1"]

## -----------------------------------------------------------------------------
inter_arrival2 <- res2$inter_arrival
type2 <- res2$type

## -----------------------------------------------------------------------------
logLik(hspec2, inter_arrival = inter_arrival2, type = type2)

## ----warning=FALSE------------------------------------------------------------
hspec0 <- hspec2
mle <- hfit(hspec0, inter_arrival = inter_arrival2, type = type2)
summary(mle)

## -----------------------------------------------------------------------------
coef(mle)

## -----------------------------------------------------------------------------
miscTools::stdEr(mle)

## ----shared-matrices----------------------------------------------------------
estimated2 <- fitted_hspec(mle)
estimated2@alpha
estimated2@beta

## ----eval=FALSE---------------------------------------------------------------
# inter_arrival <- c(0, diff(c(0, event_times)))
# type <- c(0, event_types)
# mark <- c(0, event_marks)

## ----stats-optimizer, eval=TRUE-----------------------------------------------
library(emhawkes)
poisson <- new("hspec", mu = 0.5, alpha = 0, beta = 1)
fit_stats <- hfit(poisson, c(0, 0.2, 0.3, 0.5), optimizer = "stats",
                  lambda_component0 = 0, fixed = c(FALSE, TRUE, TRUE),
                  control = list(maxit = 500, reltol = 1e-12))
summary(fit_stats)
fit_diagnostics(fit_stats)$computation$optimizer

