## ----setup, include=FALSE, message=FALSE--------------------------------------
##library(knitr)
knitr::opts_chunk$set(fig.align='center')
library(RTMB)
set.seed(1)
formals(MakeADFun)$silent <- TRUE
## Utilites for graph visulization
library(igraph)
## To tweak the vertex names
iformat <- function(n, to=c("superscript", "subscript")) {
    to <- match.arg(to)
    digits <- strsplit(as.character(n),"")[[1]]
    if (to=="superscript")
        x <- paste0("0x", c("2070", "00B9", "00B2", "00B3", paste0("207",4:9)))
    else
        x <- paste0("208", 0:9)
    names(x) <- 0:9
    intToUtf8(x[digits])
}
addindex <- function(x, to) {
    paste0(x, sapply(seq_along(x), iformat, to))
}
showGraph <- function(F) {
    g <- F$graph()
    colnames(g)[colnames(g)=="Inv"] <- "X"
    colnames(g)[colnames(g)=="Dep"] <- "Y"
    colnames(g) <- addindex(colnames(g), "sup")
    G <- graph_from_adjacency_matrix(g)
    oldpar <- par(mar=c(0,0,0,0),oma=c(0,0,0,0))
    on.exit(par(oldpar))
    ##plot(G, vertex.size=17, layout=layout_as_tree)
    plot(G, vertex.size=23, layout=layout_as_tree)
}

## -----------------------------------------------------------------------------
f <- function(x) exp( x[1] + 1.23 * x[2] )

## -----------------------------------------------------------------------------
F <- MakeTape(f, numeric(2))
F

## -----------------------------------------------------------------------------
F$methods()

## -----------------------------------------------------------------------------
F$print()

## -----------------------------------------------------------------------------
F(3:4)

## -----------------------------------------------------------------------------
F$print()

## -----------------------------------------------------------------------------
F$jacobian(3:4)

## -----------------------------------------------------------------------------
F$print()

## -----------------------------------------------------------------------------
F$graph()

## ----fig.cap="Operator graph of test function"--------------------------------
showGraph(F)

## -----------------------------------------------------------------------------
G <- MakeTape(function(x) c( F(x) , F(x) ), numeric(2))

## ----fig.cap="Operator graph of test function evaluated twice"----------------
showGraph(G)

## -----------------------------------------------------------------------------
F <- F$atomic()

## ----fig.cap="Operator graph of *atomic* test function evaluated twice"-------
G <- MakeTape(function(x) c( F(x) , F(x) ), numeric(2))
showGraph(G)

## -----------------------------------------------------------------------------
F <- MakeTape(function(x) {
    y1 <- sin(x[1] + x[2])
    y2 <- sin(x[1] + x[2])
    y3 <- cos(x[1] + x[2])
    y1+y2
}, numeric(2))

## ----fig.cap="Tape of function"-----------------------------------------------
showGraph(F)

## -----------------------------------------------------------------------------
F$simplify("eliminate")

## ----fig.cap="Tape of function after eliminate"-------------------------------
showGraph(F)

## -----------------------------------------------------------------------------
F$simplify("optimize")

## ----fig.cap="Tape of function after optimize"--------------------------------
showGraph(F)

## -----------------------------------------------------------------------------
D <- MakeTape(diff, numeric(5))
D

## -----------------------------------------------------------------------------
J <- D$jacfun()
J

## -----------------------------------------------------------------------------
J(1:5)

## -----------------------------------------------------------------------------
Js <- D$jacfun(sparse=TRUE)
Js

## -----------------------------------------------------------------------------
Js(1:5)

## -----------------------------------------------------------------------------
f <- function(X) X %*% X

## ----fig.cap="Plain matrix multiply: Expands all operations"------------------
TapeConfig(atomic="disable")
F <- MakeTape(f, matrix(0, 2, 2))
showGraph(F)

## ----fig.cap="Atomic matrix multiply: Collapses to a single operation. The constants are the matrix dimensions which are represented as additional inputs."----
TapeConfig(atomic="enable")
F <- MakeTape(f, matrix(0, 2, 2))
showGraph(F)

## -----------------------------------------------------------------------------
F <- MakeTape(function(x) {
    a <- sin(x)
    G <- MakeTape(function(y) {
        a * y
    }, numeric(1))
    DG <- G$jacfun()
    DG(x * x)
}, numeric(1))

## -----------------------------------------------------------------------------
## Negative log of the integrand
f <- MakeTape(function(x)x[1]^2+x[2]^2, numeric(2))
show(f)

## -----------------------------------------------------------------------------
## Integrate x[2] out and return the negative log of the result
F <- f$laplace(2)
show(F)

## -----------------------------------------------------------------------------
F(3)
-log(integrate(function(x2)exp(-3^2-x2^2), -Inf, Inf)$value)

## -----------------------------------------------------------------------------
## minimize wrt. x[2] and return optimum as function of x[1]
S <- f$newton(2)
show(S)

## -----------------------------------------------------------------------------
S(3)

## -----------------------------------------------------------------------------
rdblpois <- function(n, muX=5, muN=10) {
    replicate(n, sum(rpois(rpois(1, muN), muX)))
}
set.seed(1)
x <- rdblpois(100)

## -----------------------------------------------------------------------------
Kpois <- function(t, mu) mu * (exp(t) - 1)

## -----------------------------------------------------------------------------
Kpois2 <- function(s, muX, muN) Kpois(Kpois(s, muX), muN)

## -----------------------------------------------------------------------------
nldens <- function(obs, muX, muN) {
    ## 1
    F <- MakeTape(function(s) {
        K <- Kpois2(s, muX, muN)  ## CGF
        K <- K - s * obs          ## SPA adjustment
        sum(K)
    }, rep(1, length(obs)))
    ## 2
    L <- F$laplace(1:length(obs), SPA=TRUE)
    ## 3
    L(numeric(0))
}

## -----------------------------------------------------------------------------
obj <- MakeADFun(function(p) nldens(x, p$muX, p$muN), list(muX=1, muN=1), silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdreport(obj)

## -----------------------------------------------------------------------------
f <- function(x) {
    ## Get real/imag part
    xreal <- x[1]
    ximag <- x[2]
    ## Construct AD complex number
    z <- as.complex(xreal) + 1i * as.complex(ximag)
    ## Return real numbers
    Mod(exp(z))
}

## -----------------------------------------------------------------------------
F <- MakeTape(f, numeric(2))
F

## -----------------------------------------------------------------------------
## Using R complex arithmetic
f(1:2)

## -----------------------------------------------------------------------------
## Using the tape representation
F(1:2)

## -----------------------------------------------------------------------------
F$print()

