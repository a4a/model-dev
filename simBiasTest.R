library(FLa4a)
data(ple4)
data(ple4.indices)
err <- 0.05
nits <- 100

fit <-  a4aSCA(ple4, ple4.indices, qmodel=list(~s(age, k=4), ~s(age, k=4), ~s(age, k=3)))
obj <- simulate(fit, nits)

stk.fit <- stock.n(fit)
stk.sim <- stock.n(obj)
stk.rat <- iterMedians(stk.sim)/stk.fit
f.fit <- harvest(fit)
f.sim <- harvest(obj)
f.rat <- iterMedians(f.sim)/f.fit
idx.fit <- index(fit)[[1]]
idx.sim <- index(obj)[[1]]
idx.rat <- iterMedians(idx.sim)/idx.fit
idx2.fit <- index(fit)[[2]]
idx2.sim <- index(obj)[[2]]
idx2.rat <- iterMedians(idx2.sim)/idx2.fit
idx3.fit <- index(fit)[[3]]
idx3.sim <- index(obj)[[3]]
idx3.rat <- iterMedians(idx3.sim)/idx3.fit

max(stk.rat)-min(stk.rat) < err
max(f.rat)-min(f.rat) < err
max(idx.rat)-min(idx.rat) < err
max(idx2.rat)-min(idx2.rat) < err
max(idx3.rat)-min(idx3.rat) < err

xyplot(data~year|age, groups=qname, data=FLQuants(iterMedians(stk.sim), stk.fit), type="l", scales=list(y=list(relation="free")))
xyplot(data~year|age, groups=qname, data=FLQuants(iterMedians(f.sim), f.fit), type="l", scales=list(y=list(relation="free")))
xyplot(data~year|age, groups=qname, data=FLQuants(iterMedians(idx.sim), idx.fit), type="l", scales=list(y=list(relation="free")))
xyplot(data~year|age, groups=qname, data=FLQuants(iterMedians(idx2.sim), idx2.fit), type="l", scales=list(y=list(relation="free")))
xyplot(data~year|age, groups=qname, data=FLQuants(iterMedians(idx3.sim), idx3.fit), type="l", scales=list(y=list(relation="free")))


# WCSAM
fit <-  a4aSCA(ple4, ple4.indices[1], qmodel=list(~s(age, k=4)))
stk <- ple4 + fit
fits <- simulate(fit, nits, 1234)
stks <- propagate(ple4, nits) + fits 
idxs <- ple4.indices[1]
index(idxs[[1]]) <- index(fits)[[1]]

lst <- lapply(split(1:nits, 1:nits), function(x){
	out <- try(a4aSCA(iter(stks, x), FLIndices(iter(idxs[[1]], x)), fit="MP")) 
	if(is(out, "try-error")) NULL else out
})

stks2 <- stks
for(i in 1:nits){
	iter(catch.n(stks2), i) <- catch.n(lst[[i]])
	iter(stock.n(stks2), i) <- stock.n(lst[[i]])
	iter(harvest(stks2), i) <- harvest(lst[[i]])
} 
catch(stks2) <- computeCatch(stks2) 
stock(stks2) <- computeStock(stks2) 

plot(FLStocks(orig=stk, sim=stks, fitsim=stks2))

