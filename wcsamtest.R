library(FLa4a)
library(parallel)
ncores <- 15
nits <- ncores*40
options(mc.cores=ncores)
data(ple4)
data(ple4.indices)

fit <- sca(ple4, ple4.indices[1], fit="assessment") 
stk <- ple4 + fit

#====================================================================
# option 1: using pars to simulate
#====================================================================
fits <- simulate(fit, nits)
stks <- ple4 + fits 
idxs <- ple4.indices[1]
index(idxs[[1]]) <- index(fits)[[1]]
lst <- mclapply(split(1:nits, 1:nits), function(x){
	out <- try(sca(iter(stks, x), FLIndices(iter(idxs[[1]], x)))) 
	if(is(out, "try-error")) NULL else out
})

lst <- lst[!unlist(lapply(lst, is.null))]

stks1 <- stks
for(i in 1:length(lst)){
	iter(catch.n(stks1), i) <- catch.n(lst[[i]])
	iter(stock.n(stks1), i) <- stock.n(lst[[i]])
	iter(harvest(stks1), i) <- harvest(lst[[i]])
} 
catch(stks1) <- computeCatch(stks1) 
stock(stks1) <- computeStock(stks1) 

# plot
plot(FLStocks(orig=stk, sim=stks, fitsim=stks1), auto.key=list(lines=TRUE, points=FALSE, columns=3), main="simulate")

#====================================================================
# option 2: adding lognormal error with genFLQuant
#====================================================================
cv <- 0.2
index(idxs[[1]]) <- genFLQuant(index(fit)[[1]], cv=cv)
lst <- mclapply(split(1:nits, 1:nits), function(x){
	out <- try(sca(iter(stks, x), FLIndices(iter(idxs[[1]], x)))) 
	if(is(out, "try-error")) NULL else out
})

lst <- lst[!unlist(lapply(lst, is.null))]

stks2 <- stks
for(i in 1:length(lst)){
	iter(catch.n(stks2), i) <- catch.n(lst[[i]])
	iter(stock.n(stks2), i) <- stock.n(lst[[i]])
	iter(harvest(stks2), i) <- harvest(lst[[i]])
} 
catch(stks2) <- computeCatch(stks2) 
stock(stks2) <- computeStock(stks2) 

# plot
plot(FLStocks(orig=stk, sim=stks, fitsim=stks2), auto.key=list(lines=TRUE, points=FALSE, columns=3), main="genFLQuant")

#====================================================================
# option 3: adding random noise
#====================================================================
cv <- 0.2
index(idxs[[1]]) <- rnorm(nits, index(fit)[[1]], 0.2*index(fit)[[1]])
lst <- mclapply(split(1:nits, 1:nits), function(x){
	out <- try(sca(iter(stks, x), FLIndices(iter(idxs[[1]], x)))) 
	if(is(out, "try-error")) NULL else out
})

lst <- lst[!unlist(lapply(lst, is.null))]

stks3 <- stks
for(i in 1:length(lst)){
	iter(catch.n(stks3), i) <- catch.n(lst[[i]])
	iter(stock.n(stks3), i) <- stock.n(lst[[i]])
	iter(harvest(stks3), i) <- harvest(lst[[i]])
} 
catch(stks3) <- computeCatch(stks3) 
stock(stks3) <- computeStock(stks3) 

# plot
plot(FLStocks(orig=stk, sim=stks, fitsim=stks3), auto.key=list(lines=TRUE, points=FALSE, columns=3), main="rnorm")

#====================================================================
# all together
#====================================================================
plot(FLStocks(orig=stk, simpars=stks1, simlog=stks2, simnorm=stks3), auto.key=list(lines=TRUE, points=FALSE, columns=4))


