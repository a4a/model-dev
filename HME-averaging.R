
## 
# this script loads the FLa4a library and writes data to a directory "test"
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

f0 <- a4aInternal(stock = stock, indices = indices, fit = "Ext")
AIC(f0)
AIC(a4aFit(f0))
fit <- stock + f0
fit <- stock + a4aFit(f0)

# simulate from model - only for Ext type...
b <- mvrnorm(3, f0)
b <- cbind(b, f0 @ baseLvlPars)

## calculate likelihood
calcLogLik(b, f0, stock)
-1 * f0 @ fitSumm[2]
llik0 <- calcLogLik(mvrnorm(3, f10), f10, stock)

# so that seems to work.... Not lets try the HME!

# we need samples from the posterior:

f1 <- a4aInternal(~ te(age, year, k = c(4, 25)), 
          list(~ s(age, k = 4), ~ s(age, k = 4), ~ age), 
          vmodel = list(~ s(age, k = 3), ~1, ~1, ~1),
          stock = stock, indices = indices, fit = "MCMC")
f10 <- a4aInternal(~ te(age, year, k = c(4, 25)), 
          list(~ s(age, k = 4), ~ s(age, k = 4), ~ age), 
          vmodel = list(~ s(age, k = 3), ~1, ~1, ~1),
          stock = stock, indices = indices, fit = "Ext")

llik1 <- calcLogLik(t(f1 @ pars), f10, stock)

f2 <- a4aInternal(~ te(age, year, k = c(4, 25)), 
          list(~ te(age, year, k = c(4, 3)), ~ te(age, year, k = c(4, 3)), ~ s(year, k = 3, by = age)),
          vmodel = list(~ s(age, k = 3), ~1, ~1, ~1),
          stock = stock, indices = indices, fit = "MCMC")
f20 <- a4aInternal(~ te(age, year, k = c(4, 25)), 
          list(~ te(age, year, k = c(4, 3)), ~ te(age, year, k = c(4, 3)), ~ s(year, k = 3, by = age)),
          vmodel = list(~ s(age, k = 3), ~1, ~1, ~1),
          stock = stock, indices = indices, fit = "Ext")

llik2 <- calcLogLik(t(f2 @ pars), f20, stock)
         
          
f3 <- a4aInternal(~ te(age, year, k = c(4, 25)), 
          list(~ te(age, year, k = c(4, 3)), ~ te(age, year, k = c(4, 3)), ~ s(year, k = 3, by = age)),
          vmodel = list(~ s(age, k = 3), ~1, ~1, ~1),
          srmodel = ~ bevholt(CV = 0.5),
          stock = stock, indices = indices, fit = "MCMC")
f30 <- a4aInternal(~ te(age, year, k = c(4, 25)), 
          list(~ te(age, year, k = c(4, 3)), ~ te(age, year, k = c(4, 3)), ~ s(year, k = 3, by = age)),
          vmodel = list(~ s(age, k = 3), ~1, ~1, ~1),
          srmodel = ~ bevholt(CV = 0.5),
          stock = stock, indices = indices, fit = "Ext")

llik3 <- calcLogLik(t(f3 @ pars), f30, stock)

# calculate HME for each chain
ps <- 1/c(mean(exp(-1*llik1)), mean(exp(-1*llik2)), mean(exp(-1*llik3)))
ps <- ps / sum(ps)
round(ps, 2)

# make a list of models
mod <- sample(1:3, 1000, weights = ps, replace = TRUE)

# and so on!!




