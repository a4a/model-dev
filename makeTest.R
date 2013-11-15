
#install_github("FLash", username = "flr")
#install_github("FLBRP", username = "flr")

## 
# this script loads the FLa4a library and writes data to a directory "test"
##
options(width=110)
#options(width=160)

library(FLa4a)

# get data and add some bits
data(ple4)
data(ple4.indices)


# fit two of the simplest models with a small amount of data
# for now we only deal with 1 area, 1 unit.
stock <- ple4[,paste(1995:2008)]
indices <- FLIndices(ple4.indices[[1]][,paste(1995:2008)])

f1 <- a4aInternal(~ 1, stock = stock, indices = indices, fit = "MP")
f1 <- a4aInternal(~ 1, stock = stock, indices = indices, fit = "assessment")

f2 <- a4aInternal(~ year, stock = stock, indices = indices, fit = "assessment")


# choose models to average over
mods <- list(f1, f2)

# get design matrices
X <- lapply(mods, getX)

#
obs <- c(c(catch.n(stock)), unlist(lapply(indices, function(x) c(index(x)))))

# do a quick independent weighted sample from the posterior
f1mc <- a4aInternal(~ 1, stock = stock, indices = indices, fit = "MCMC")
f2mc <- a4aInternal(~ year, stock = stock, indices = indices, fit = "MCMC")




# propose some params from models
b.sim <- lapply(mods, mvrnorm, n = 1)

# get prediction of catch and indices
preds <- lapply(1:length(mods), function(i) makePrediction(b.sim[[i]], X[[i]], stock, mods[[i]]))
predVars <- lapply(1:length(mods), function(i) makeVarPrediction(b.sim[[i]], X[[i]], stock, mods[[i]]))


# calculate likelihood
llik <- 

#



# iterations
niters <- 2
nseasons <- 4
endyear <- 2008

# try a forecast?
if (endyear > 2008) {
  require(plyr)
  require(FLBRP)
  stock <- fwdWindow(ple4, FLBRP(ple4), end = endyear)
} else {
  stock <- ple4
}

# add in iters and seasons
stock <- qapply(qapply(stock, propagate, niters), 
           function(x, n) {
             dms <- dimnames(x)
             dms $ season <- paste(1:n)
             xp <- rep(c(aperm(x, c(1,2,3,5,6,4))), n)
             dim(xp) <- c(dim(x)[c(1,2,3,5,6)], n)
             out <- aperm(xp, c(1,2,3,6,4,5))
             dimnames(out) <- dms
             FLQuant(out)
           }, n = nseasons)



indices <- ple4.indices  # note no iters in indices

# quick test
out <- a4a(srmodel = ~ geomean(CV = 0.3), stock = stock, indices = indices, fit = "MP")

fbar(ple4 + out[,,,,,1])
fbar(ple4)




exp(yearMeans(log(stock.n(out)[1,-1]))) # first recruitment is not in SR relationship!!
stock.n(out)[1,paste(endyear)]




stock <- ple4

indices <- ple4.indices
#index.var(indices[[1]])[] <- 1
#index.var(indices[[2]])[] <- 1
#index.var(indices[[3]])[] <- 2

varslt <- catch.n(stock)
varslt[] <- 1
catch.n(stock) <- FLQuantDistr(catch.n(stock), varslt)
    
covar <- list(btype = replace(stock.n(stock), TRUE, 1:10), 
              ytype = replace(catch(stock), TRUE, 1:52), 
              atype = replace(stock.n(stock)[1:7,1], TRUE, c(1:7)))

# set up model

fmodel <- ~ s(age, k = 4) + btype
qmodel <- ~ fleet + s(age, k = 3, by = breakpts(year, 2000))
rmodel <- ~ bevholt(amodel = ~ s(year, k = 3), CV = -1) + year + ytype
vmodel <- ~ fleet

# do not run only write data to file.
# a4aInternal(fmodel, qmodel, rmodel, vmodel, stock = stock, indices = indices, covar = covar, wkdir = "test", fit = FALSE)

options(width=110)
#options(width=160)

library(FLa4a)

# get data and add some bits
data(ple4)
data(ple4.indices)

endyear <- 2008

if (endyear > 2008) {
  require(plyr)
  require(FLBRP)
  stock <- fwdWindow(ple4, FLBRP(ple4), end = endyear)
} else {
  stock <- ple4
}



indices <- ple4.indices

covar <- list(a1 = replace(stock.n(stock)[,1], TRUE, c(1, rep(0, 9))))

fmodel <- ~ te(age, year, k = c(4, 20)) + s(year, k = 5, by = as.numeric(age == 1))
qmodel <- list(~ te(age, year, k = c(4,3)), ~ te(age, year, k = c(4, 3)), ~ s(year, k = 3, by = age))
srmodel <- ~ geomean(CV = 0.5)
n1model <- ~ factor(age)
vmodel <- list(~ s(age, k = 3), ~1, ~1, ~ 1)
wkdir <- "test"
fit <- "MP"
verbose <- FALSE


out <- a4aInternal(fmodel, qmodel, srmodel, n1model, vmodel, stock = stock, indices = indices, wkdir = wkdir, fit = fit)

fitstk <- stock + out
plot(FLStocks(fitstk, ple4))


#ssb(stock + a4aFit(out))
#ssb(stock + out)
#ssb(propagate(stock, 100) + out)

fitstk <- propagate(stock, 1000) + out

dev.off()
lattice.options(default.theme = geta4aLatticeOptions())

plotIters(fbar(fitstk))

plotIters(stock.n(fitstk))




library(proftools)

Rprof(tmp <- tempfile())
out <- a4aInternal(fmodel, qmodel, srmodel, n1model, vmodel, stock = stock, indices = indices, wkdir = wkdir, fit = fit)
Rprof()
flatProfile(readProfileData(tmp))
unlink(tmp)





mat <- cov2cor(stkmodel(pars(out)) @ vcov)
mat[lower.tri(mat, diag = TRUE)] <- 0
image(Matrix(mat), main = "Stock param correlations")



stf <- attr(out, "stuff")




plot(stock.n(out)[1,])

image(Matrix(drop(harvest(out)@.Data)))
plot(log(drop(stock.n(out)@.Data)[,1]))


rec <- c(stock.n(out)[1,-1])
ssb <-  colSums(drop((stock.n(out) * stock.wt(stock) * mat(stock))@.Data))[1:length(rec) + 1]
plot(ssb, rec)
sd(log(ssb) - log(rec))



if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

    N <- 50                                 # length of the time series
    series <- cumsum(rnorm(N))
#    
# Here is the code that collects bootstrap values of
# the auto-correlation estimate:
    k <- 6                                  # size of moving blocks

x <- sapply(1:999, function(irep) {                   # the bootstrap loop
      series.bt <- rep(NA,N)                # local vector for a bootstrap replication
      for(i in 1:ceiling(N/k)) {            # fill the vector with random blocks
        endpoint <- sample(k:N, size=1)     # by randomly sampling endpoints
        series.bt[(i-1)*k+1:k] <- series[endpoint-(k:1)+1] # and copying blocks
      }
      series.bt[1:N]           # trim overflow when k doesn't divide N
    })
    


