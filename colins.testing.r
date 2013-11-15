
library(knitr)
knit2pdf("a4aAssessmentMethodology.Rnw")




library(FLa4a)

data(ple4)
data(ple4.indices)

stock <- ple4
#attr(catch.n(stock), "var") <- replace(catch.n(stock), TRUE, NA)

indices <- ple4.indices
#index.var(indices[[2]])[,1] <- 1
#index.var(indices[[3]])[] <- 2

fmodel <- ~ s(age, k = 3) + factor(year)
#fmodel <- ~ s(replace(age, age > 5, 5), k = 3, by = breakpts(year, 2000))
#fmodel <- ~ trawl(plateau = 5, selectivity = "variable") + atype


#qmodel = ~ fleet + s(age, k = 3, by = factor(year >= 2000, labels = c("y<2000", "y>=2000")))
#qmodel = ~ atype
qmodel = ~ fleet + s(age, k = 4, by = fleet)


rmodel = ~ bevholt(amodel = ~ s(year, k = 3), CV = 0.3)


vmodel = ~ fleet

wkdir = "test"
verbose = FALSE

covar <- list(btype = replace(stock.n(stock), TRUE, 1:10), 
              ytype  = replace(catch(stock), TRUE, 1:52), 
              atype = replace(stock.n(stock)[1:7,1], TRUE, c(1:7)))


fmodel <- ~ trawl(plateau = 5, selectivity = "variable") + btype
qmodel <- ~ fleet + s(age, k = 3, by = breakpts(year, 2000))
rmodel <- ~ bevholt(amodel = ~ s(year, k = 3), CV = 0.3) + year + ytype
vmodel <- ~ fleet

a4a(fmodel, qmodel, rmodel, vmodel, stock = stock, indices = indices, covar = covar, wkdir = wkdir)



Rprof(tmp <- tempfile())
a4aFit(fmodel, qmodel, rmodel, vmodel, stock = stock, indices = indices)
Rprof()
flatProfile(readProfileData(tmp))
unlink(tmp)



fmodel <- ~ 1
qmodel <- list( ~ 1)
rmodel <- ~ factor(year) 
stock <- hakeGSA7
indices <- hakeGSA7.idx
srrCV = -1
fmodel.extra = NULL
qmodel.extra = NULL
rmodel.extra = NULL
vmodel = list(~1, ~1)
vmodel.extra = NULL
wkdir = NULL
verbose = FALSE
MCMC = FALSE
NMCMC = 1000





library(FLa4a)
options(width = 140)

data(ple4)
data(ple4.indices)

stk <- ple4
idx <- ple4.indices

out <- a4a(~ 1, list(~ 1, ~ 1, ~ 1), stock = stk, indices = idx)
index(out)


catch.res <- log(catch.n(stk)) - log(catch.n(out))
index.res <- lapply(1:length(idx), function(i) log(index(idx[[i]])) - log(index(out)[[i]]))



# residuals about what the model predicts the catch should be
xyplot(data ~ year | paste("age", age), data = catch.lres(object), type = c("g","p","smooth"), 
     lty = 2, col = 1, ylab = "standardised residual", as.table = TRUE, main = "log (Catch) Residuals")























reclogSE <- function(object) 
{
  opts <- options(contrasts = c(unordered = "contr.sum", ordered = "contr.poly")) 

  x <- stock.n(object)
  Xr <- getX(object @ models $ rmodel, 
             data.frame(year = as.numeric(colnames(x)), x = 1))

  options(opts) # reset options
  
  coefs <- coef(object)      

  # get variance matrix of recruit params
  Vr <- vcov(object)[grep("rpar", names(coefs)), grep("rpar", names(coefs))]
  
  # linearly transform variance using Xr to get variance mat of recruitment
  # then return standard error (on log scale recruitment)
  sqrt(diag(Xr %*% (Vr %*% t(Xr)))) 
}

FlogSE <- function(object) 
{
  # set up observation data frame
  x <- stock.n(object)
  f.df <- expand.grid(age = as.numeric(rownames(x)),
                         year = as.numeric(colnames(x)), 
                          x = 1)[c(2,1,3)]

  opts <- options(contrasts = c(unordered = "contr.sum", ordered = "contr.poly")) 

  Xf <- getX(object @ models $ fmodel, f.df)

  options(opts) # reset options
  
  coefs <- coef(object)      
  f <- coefs[grep("fpar", names(coefs))]

  Vf <- vcov(object)[grep("fpar", names(coefs)), grep("fpar", names(coefs))]
  
  out <- sqrt(diag(Xf %*% (Vf %*% t(Xf))))
  matrix(out, dim(x)[1], dim(x)[2])
}


ssblogSE(object)
fbarlogSE(object)
reclogSE(object)



covmat <- cov2cor(vcov(object))
diag(covmat) <- 0
image(Matrix(covmat))





calcLogLik(coef(object), object, stock = hakeGSA7, indices = hakeGSA7.idx)
logLik(object)

calcLogLik(t(sapply(1:20, function(x) coef(object))), object, stock = hakeGSA7, indices = hakeGSA7.idx)



# very simple model averaging

mcmc <- a4aFit(~ s(age, k = 3) + s(year, k = 4), list(~ s(age, k = 3)), stock = hakeGSA7, indices = hakeGSA7.idx, MCMC = TRUE, NMCMC = 5000, verbose = TRUE)
colnames(mcmc) <- names(coef(object))

# vectorised
llik <- calcLogLik(mcmc, object = object, stock = hakeGSA7, indices = hakeGSA7.idx)

dmc <- as.data.frame(mcmc)
dmc $ llik <- llik

plot(dmc[grep("logSd|llik", names(dmc))])

psim <- mvrnorm(100, coef(object), vcov(object))




#data(ple4)
#data(ple4.indices)

#object <- a4aFit(~ s(age, k = 3) + s(year, k = 4), list(~ s(age, k = 3), ~ s(age, k = 3), ~ age), stock = ple4, indices = ple4.indices)
#plot(object, ple4, what = "F", Ftext = TRUE)



if (0) {
filen <- file("bma1/a4a.psv", "rb")
nopar <- readBin(filen, what = integer(), n = 1)
mcmc <- readBin(filen, what = numeric(), n = nopar * 10000)
mcmc <- matrix(mcmc, byrow = TRUE, ncol = nopar)
close(filen)

fit2 <- a4aFit(~ s(age, k=3) + s(year, k = 6), list(~ s(age, k=3)), stock = hke, indices = hke.idx, wkdir = "bma2")

filen <- file("bma2/a4a.psv", "rb")
nopar <- readBin(filen, what = integer(), n = 1)
mcmc2 <- readBin(filen, what = numeric(), n = nopar * 10000)
mcmc2 <- matrix(mcmc2, byrow = TRUE, ncol = nopar)
close(filen)


#==============================================================================
#   FLData package testing
#==============================================================================



#source("http://www.math.ntnu.no/inla/givemeINLA.R")

library(devtools)
library(testthat)
library(roxygen2)

roxygenize("../packages/FLData")
pkg <- as.package("../packages/FLData")
build(pkg)
#check(pkg)
install(pkg)

library(FLData)


data(ple4)
object <- harvest(ple4[,1:5])

x1 <- genFLQuant(object, method = "ac", cv = 0.1)
x2 <- genFLQuant(object, method = "rw", cv = 0.1)

log(x1)
log(x2)




mvrnormQ <- function (n, mu, Q, tol = 1e-06, return.density = FALSE) {
  # 
  #  Sample from a GMRF with mean mu and precision Q
  #  based on mvrnorm from library MASS
  #
  p <- length(mu)
  if (!all(dim(Q) == c(p, p))) stop("incompatible arguments")
  eQ <- eigen(Q, symmetric = TRUE, EISPACK = TRUE)
  ev <- eQ $ values
  if (!all(ev >= -tol * abs(ev[1L])))
    stop("'Q' is not positive definite")
  X <- matrix(rnorm(p * n), nrow = n)
  evinv <- ifelse(ev < tol, 0, 1/ev)
  
  X <- drop(mu) + (eQ $ vectors %*% diag(sqrt(evinv), p)) %*% t(X)

  if (return.density) {
    ldet <- sum(log(ev[ev>0]))
    x <- X - mu
    X <- list(x = X, log.dens = 0.5 * ldet - 0.5 * t(x) %*% (Q%*%x))
  }
  
  X
}

Dfunc <- function(n) {
  D <- diag(n)
  diag(D[-1,]) <- -1
  D[-1,]
}


Qfunc <- function(n, type = "rw1", weights = 1) {
  type <- match.arg(type, c("rw1", "rw2"))
  
  switch(type,
      rw1 = {
        if (!(length(weights) %in% c(1, n-1))) 
          stop("supplied weights should be of length n - 1 forRW1 model")
        D <- Dfunc(n)
        Q <- t(D) %*% (diag(weights, nrow = n-1) %*% D)  
      },
      rw2 = {
        if (!(length(weights) %in% c(1, n-2)))
          stop("supplied weights should be of length n - 2 for RW2 model")
        D <- Dfunc(n-1) %*% Dfunc(n)
        Q <- t(D) %*% (diag(weights, nrow = n-2) %*% D) 
      })
  Q
}


selQfunc <- function(n, type = "smooth", edf) {
  type <- match.arg(type, c("smooth", "flattop"))
  
  switch(type,
      smooth = {
        Q <- Qfunc(n, type = "rw2")  
      },
      flattop = {
        weights <- rep(0, n - 1)
        filt <- floor(n/2):(n-1)
        weights[filt] <-  exp( 5 * (seq_along(filt) - 1) /(length(filt - 1)) )
        Q <- Qfunc(n, type = "rw2", weights = weights[-1] + 1) + Qfunc(n, type = "rw1", weights = weights) 
      })
  Q * n * n
}

plotit <- function(x, size = c(0.5, 1, 0.5), lty = c(2, 1, 2), 
                      probs = c(.025, .5, .975)) {
  df1 <- data.frame(y     = c(t(apply(x, 1, quantile, probs))), 
                   x     = rep(1:nrow(x), length(probs)), 
                   group = factor(rep(probs, each = nrow(x))))
  worm <- sample(ncol(x), 5)
  df2 <- data.frame(y    = c(x[,worm]), 
                   x     = rep(1:nrow(x), length(worm)), 
                   group = factor(rep(worm, each = nrow(x))))             
  p <- ggplot(df1)
  p <- p + geom_line(aes(x = x, y = y, 
              group = group, size = group, lty = group)) + 
              scale_size_manual(values = size, name = "Quantile") + 
              scale_linetype_manual(values = lty, name = "Quantile")
  
  p <- p + xlab("Year") + ylab("")
  p <- p + geom_line(aes(x, y, group = group, colour = group), data = df2)
  p
}

Q <- selQfunc(100, type = "flattop")
x <- mvrnormQ(1000, rep(0, nrow(Q)), Q, return.density = TRUE)
plotit(x $ x)



  

# selectivity function

getSelectivityFunction <-
    function(npar = 2, ages = 1:9, age.plat = max(ages), age.ref = 2:4)
{
  X <- model.matrix(gam(y ~ s(ages, k = max(3, npar + 1)), data = list(y = ages, ages = ages)))
  X <- X[,-1]
  if (npar == 1) X <- X[,-2]
  plat <- which(ages >= age.plat)
  X[plat,] <- matrix(X[plat[1],], ncol = npar, nrow = length(plat), byrow = TRUE)
  function(par) exp(drop(X %*% par))
}






selFun(c(1,1))



qset <- matrix(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.5, 0.5, 0.5,
        0.1, 0.15, 0.35, 0.5, 0.2, 0.1, 0.1, 0.1), ncol = 2)

set <- 2
g1 <- gam(log(y) ~ s(x, k = 5, fx = TRUE), data = list(y = qset[,set], x = 1:8))

i1 <- inla(log(y) ~ f(x, model = "rw2"), data = list(y = qset[,set], x = 1:8))

plot(1:8, exp(fitted(g1)), type = "l", ylim = c(0, max(exp(fitted(g1)), qset[,set])))
lines(1:8, exp(i1 $ summary.random $ x $ mean + i1 $ summary.fixed[1,"mean"]))
points(1:8, qset[,set])



}




