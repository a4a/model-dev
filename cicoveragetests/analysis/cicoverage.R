library(FLa4a)

nit <- 250
data(ple4)
data(ple4.index)
fit <- sca(ple4, FLIndices(ple4.index), fit="assessment")

idx <- ple4.index
stk <- ple4 + fit
flq <- fbar(stk)[rep(1,3)]
flq[2] <- ssb(stk)
flq[3] <- rec(stk)

index(idx) <- rlnorm(nit, log(index(idx)), abs(log(index(idx))*0.4))
catch.n(stk) <- rlnorm(nit, log(catch.n(stk)), abs(log(catch.n(stk))*0.2))

flqres <- catch.n(stk)[1:3]
flqres[] <- NA
dimnames(flqres)[[1]] <- c("fbar", "ssb", "rec") 

for(x in 1:nit){
	cat(x, ", ")
	stk0 <- iter(stk, x)
	fit0 <- sca(stk0, FLIndices(iter(idx, x)), fit="assessment")
	stk0 <- stk0 + simulate(fit0, 250)
	flq0 <- fbar(stk0)[rep(1,3)]
	flq0[2] <- ssb(stk0)
	flq0[3] <- rec(stk0)
	flqres[,,,,,x] <- apply(flq0, c(1,2), quantile, 0.05) < flq & apply(flq0, c(1,2), quantile, 0.95) > flq 
}

res <- iterMeans(flqres)[drop=TRUE]


