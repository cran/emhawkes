## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


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
# first and third columns are the same
cbind(res1$lambda[1:5], res1$lambda_component[1:5], mu1 + res1$lambda_component[1:5])

## -----------------------------------------------------------------------------
# second and third columns are the same
cbind(res1$lambda[1:5], res1$rambda[1:5], res1$lambda[1:5] + alpha1)

## -----------------------------------------------------------------------------
# By definition, the following two are equal:
res1$lambda[2:6]
mu1 + (res1$rambda[1:5] - mu1) * exp(-beta1 * res1$inter_arrival[2:6])

## ----warning=FALSE------------------------------------------------------------
logLik(hspec1, inter_arrival = res1$inter_arrival)

## ----warning=FALSE------------------------------------------------------------
# initial value for numerical optimization
mu0 <- 0.5; alpha0 <- 1.0; beta0 <- 1.8
hspec0 <- new("hspec", mu = mu0, alpha = alpha0, beta = beta0)
# the intial values are provided through hspec
mle <- hfit(hspec0, inter_arrival = res1$inter_arrival)
summary(mle)

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
# Under bi-variate model, there are two types, 1 or 2.
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

