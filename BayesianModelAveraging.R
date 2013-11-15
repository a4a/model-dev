

## 
# this script loads the FLa4a library and performs Bayesian Model Averaging
##
options(width=110)

library(FLa4a)

# get data and add some bits
data(ple4)
data(ple4.indices)

# fit two of the simplest models with a small amount of data
# for now we only deal with 1 area, 1 unit.
stock <- ple4
indices <- ple4.indices

f0 <- a4aInternal(~ s(age), list(~1, ~1, ~1), srmodel = ~1, n1model = ~1, stock = stock, indices = indices, fit = "Ext")

AIC(f0)
AIC(a4aFit(f0))
fit <- stock + f0
fit <- stock + a4aFit(f0)


#
# SO, now we will write a model for a simple MCMC update with fixed proposal
#

#
# this method is a bit cheaty and uses an ADMB sample of the posterior
# to get a sample from the params conditional on the model
# a step between this and the mvrnorm independence sampler is to use a copula approx to the posterior
#

 




if (0) {
#
#  this method uses a random walk step proposal for a single model
#
do.mc1 <- function(nmc = 1000, model, nburn = 0, thinby = 1, RWproposal) 
{

  # initial setup
  bmc <- matrix(NA, length(model @ baseLvlPars), nmc, dimnames = list(names(model @ baseLvlPars)))
  probs <- matrix(NA, 2, nmc, dimnames = list(c("loglik","logprior")))

  # initialisation
  b0 <- RWproposal()
  bmc[,1] <- b0
  probs[,1] <- c(calcLogLik(bmc[,1], model, stock), calcLogPrior(bmc[,1], model))
  eprob <- 0

  for (i in 2:nmc) {

    # copy forward reject move
    bmc[,i] <- bmc[,i-1]
    probs[,i] <- probs[,i-1]

    # Proposal step:
    bnew <- bmc[,i] + RWproposal()
  
    # log posterior and proposal densities
    probsnew <- c(calcLogLik(bnew, model, stock), calcLogPrior(bnew, model))

    # metropolis hastings step
    # new-posterior
    acc.prob <- exp(min(0, sum(probsnew) - (sum(probs[,i]))))

    if (runif(1) < acc.prob) {
      probs[,i] <- probsnew
      bmc[,i] <- bnew
    }
    eprob <- eprob + acc.prob
  }
  
  print(eprob/nmc * 100)
  bmc[,seq(nburn + 1, nmc, by = thinby)]
}



## running and plotting

# fit model first
f0 <- a4aInternal(~ 1, list(~1, ~1, ~1), srmodel = ~1, n1model = ~1, stock = stock, indices = indices, fit = "Ext")

# make proposal
RWproposal <- function() propmodel @ L %*% rnorm(ncol(propmodel @ L))
propmodel <- f0
propmodel @ Sigma <- propmodel @ Sigma * 2.38^2 / ncol(propmodel @ Sigma)
propmodel @ L <- chol(propmodel @ Sigma)

system.time(
  fit <- do.mc1(5000, f0, thinby = 10, nburn = 1000, RWproposal = RWproposal)  # approx 1 sec for every 100 iters
)

which <- grep("Mod:", rownames(fit))
par(mfrow = c(3,4), oma = c(0,0,0,0), mar = c(0,0,0,0) + 0.2)
for (i in which) {
  plot(fit[i,], type = "l", axes = FALSE, ann = FALSE)
  abline(h = f0 @ baseLvlPars[i], col = "red")
  box(col = grey(0.8))
  #hist(fit[i,], nclass = 50, axes = FALSE, ann = FALSE, main = "", prob = TRUE)
  #lines(density(fit[i,]))
}

}

if (0) {
#
#  this method uses a random walk step proposal and jumps between models!
#
do.rjmc1 <- function(nmc = 1000, models, nburn = 0, thinby = 1, RWproposal, propRJ) 
{

  # initial setup
  mods <- seq(models)
  npars <- sapply(models, function(x) x @ fitSumm[1])
  bmc <- matrix(NA, max(npars), nmc)
  probs <- matrix(NA, 2, nmc, dimnames = list(c("loglik","logprior")))

  # initialisation
  m0 <- sample(mods, 1)
  b0 <- RWproposal(m0) # random start point centered at zero
  bmc[1:npars[m0],1] <- b0
  probs[,1] <- c(calcLogLik(bmc[1:npars[m0],1], models[[m0]], stock), calcLogPrior(bmc[1:npars[m0],1], models[[m0]]))
  eprob <- 0

  for (i in 2:nmc) {

    # copy forward reject move
    bmc[,i] <- bmc[,i-1]
    probs[,i] <- probs[,i-1]

    # Proposal step:
    mnew <- sample(setdiff(mods, m0), 1)
    bnew <- bmc[1:npars[mnew],i] + RWproposal(mnew)
  
    # log posterior and proposal densities
    probsnew <- c(calcLogLik(bnew, models[[mnew]], stock), calcLogPrior(bnew, models[[mnew]]))

    # metropolis hastings step
    # new-posterior
    acc.prob <- exp(min(0, sum(probsnew) - (sum(probs[,i]))))

    if (runif(1) < acc.prob) {
      probs[,i] <- probsnew
      bmc[,i] <- bnew
    }
    eprob <- eprob + acc.prob
  }
  
  print(eprob/nmc * 100)
  bmc[,seq(nburn + 1, nmc, by = thinby)]
}

# fit models
f0 <- a4aInternal(~ factor(replace(age, age>8, 8)), list(~1, ~1, ~1), srmodel = ~1, n1model = ~1, stock = stock, indices = indices, fit = "Ext")
f1 <- a4aInternal(~ factor(replace(age, age>9, 9)), list(~1, ~1, ~1), srmodel = ~1, n1model = ~1, stock = stock, indices = indices, fit = "Ext")

# make proposal the posterior
models <- list(f0, f1)

RWproposal <- function(m) propmodel[[m]] @ L %*% rnorm(ncol(propmodel[[m]] @ L))
propmodel <- lapply(models, function(x) {
                x @ Sigma <- x @ Sigma * 2.38^2 / ncol(x @ Sigma)
                x @ L <- chol(x @ Sigma)
                x})

system.time(
  fit <- do.mc1(5000, f0, thinby = 10, prop = proposal)  # approx 1 sec for every 100 iters
)

which <- grep("Mod:", rownames(fit))
par(mfrow = c(3,4), oma = c(0,0,0,0), mar = c(0,0,0,0) + 0.2)
for (i in which) {
  plot(fit[i,], type = "l", axes = FALSE, ann = FALSE)
  abline(h = f0 @ baseLvlPars[i], col = "red")
  box(col = grey(0.8))
  #hist(fit[i,], nclass = 50, axes = FALSE, ann = FALSE, main = "", prob = TRUE)
  #lines(density(fit[i,]))
}

}






