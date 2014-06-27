# testing index variance for simulation

library(FLa4a)
data(ple4)
data(ple4.indices)

# number of iters
nits <- 25
# fit the model
fit <- a4aSCA(ple4, ple4.indices[1]) 
# update the stock data
stk <- ple4 + fit
# simulate controlling the random seed
fits <- simulate(fit, nits, 1234)
# update stock and index data now with iters
stks <- ple4 + fits 

# Option 1.1) implemented at the moment
idxs <- ple4.indices[1]
idx1.1 <- index(fits)[[1]]
index(idxs[[1]]) <- idx1.1 
# rerun assessment
sfit1 <- a4aSCA(stks, idxs)

# Option 1.2)
fits2 <- fit
fits2@pars <- simulate(pars(fit), nits, 1234)
preds <- predict(fits2)

sfrac <- mean(range(ple4.indices[[1]])[c("startf", "endf")])
# fraction of total mortality up to that moment
Z <- (m(ple4) + harvest(fit))*sfrac
lst <- dimnames(fit@index[[1]])
# survivors
lst$x <- stock.n(fit)*exp(-Z)
stkn <- do.call("trim", lst)
stkn <- propagate(stkn, nits)

# index 
idx1.2 <- stkn*preds$qmodel[[1]]
index(idxs[[1]]) <- idx1.2
sfit2 <- a4aSCA(stks, idxs)

# Option 2)
index(idxs[[1]]) <- index(fit)[[1]]
vp <- predict(fit)$vmodel[[2]]
idx2 <- exp(rnorm(nits, log(index(idxs[[1]])), sqrt(vp)))
index(idxs[[1]]) <- idx2 
sfit3 <- a4aSCA(stks, idxs)


