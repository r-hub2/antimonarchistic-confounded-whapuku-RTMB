## ----setup, include=FALSE, message=FALSE--------------------------------------
##library(knitr)
library(RTMB)
set.seed(1)
formals(MakeADFun)$silent <- TRUE

## -----------------------------------------------------------------------------
data(ChickWeight)

## -----------------------------------------------------------------------------
parameters <- list(
    mua=0,          ## Mean slope
    sda=1,          ## Std of slopes
    mub=0,          ## Mean intercept
    sdb=1,          ## Std of intercepts
    sdeps=1,        ## Residual Std
    a=rep(0, 50),   ## Random slope by chick
    b=rep(0, 50)    ## Random intercept by chick
)

## -----------------------------------------------------------------------------
f <- function(parms) {
    getAll(ChickWeight, parms, warn=FALSE)
    ## Optional (enables extra RTMB features)
    weight <- OBS(weight)
    ## Initialize joint negative log likelihood
    nll <- 0
    ## Random slopes
    nll <- nll - sum(dnorm(a, mean=mua, sd=sda, log=TRUE))
    ## Random intercepts
    nll <- nll - sum(dnorm(b, mean=mub, sd=sdb, log=TRUE))
    ## Data
    predWeight <- a[Chick] * Time + b[Chick]
    nll <- nll - sum(dnorm(weight, predWeight, sd=sdeps, log=TRUE))
    ## Get predicted weight uncertainties
    ADREPORT(predWeight)
    ## Return
    nll
}

## -----------------------------------------------------------------------------
obj <- MakeADFun(f, parameters, random=c("a", "b"))

## -----------------------------------------------------------------------------
opt <- nlminb(obj$par, obj$fn, obj$gr)

## -----------------------------------------------------------------------------
sdr <- sdreport(obj)
sdr

## -----------------------------------------------------------------------------
set.seed(1)
chk <- checkConsistency(obj)
chk

## -----------------------------------------------------------------------------
osa <- oneStepPredict(obj, method="fullGaussian", discrete=FALSE)
qqnorm(osa$res); abline(0,1)

## -----------------------------------------------------------------------------
f2 <- function(parms) {
    getAll(ChickWeight, parms, warn=FALSE)
    ## Optional (enables extra RTMB features)
    weight <- OBS(weight)
    ## Random slopes
    a %~% dnorm(mean=mua, sd=sda)
    ## Random intercepts
    b %~% dnorm(mean=mub, sd=sdb)
    ## Data
    predWeight <- a[Chick] * Time + b[Chick]
    weight %~% dnorm(predWeight, sd=sdeps)
    ## Get predicted weight uncertainties
    ADREPORT(predWeight)
}

## -----------------------------------------------------------------------------
obj <- MakeADFun(f2, parameters, random=c("a", "b"))

