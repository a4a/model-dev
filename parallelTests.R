nits <- 100
fit <- a4aSCA(ple4, ple4.indices[1]) 
fits <- simulate(fit, nits)
stks <- ple4 + fits 
idxs <- ple4.indices[1]
index(idxs[[1]]) <- index(fits)[[1]]
library(parallel)
options(mc.cores=4)
lst <- mclapply(split(1:nits, 1:nits), function(x){
	cat(".")
	sca(iter(stks, x), FLIndices(iter(idxs[[1]], x))) 
})

stks2 <- stks
for(i in 1:nits){
	iter(catch.n(stks2), i) <- catch.n(lst[[i]])
	iter(stock.n(stks2), i) <- stock.n(lst[[i]])
	iter(harvest(stks2), i) <- harvest(lst[[i]])
} 
catch(stks2) <- computeCatch(stks2) 
stock(stks2) <- computeStock(stks2) 

plot(FLStocks(orig=ple4, sim=stks, fitsim=stks2), auto.key=list(columns=3))

setMethod("FLQuantPoint", "FLQuant", function (object, ..., units = "NA", lowq=0.25, uppq=0.75) 
    {
        res <- new("FLQuantPoint", FLQuant(NA, dimnames = c(dimnames(object)[1:5], 
            iter = list(c("mean", "median", "var", "uppq", "lowq"))), 
            units = units))
        res[, , , , , "mean"] <- apply(object, 1:5, mean, na.rm = TRUE)
        res[, , , , , "median"] <- apply(object, 1:5, median, 
            na.rm = TRUE)
        res[, , , , , "var"] <- apply(object, 1:5, var, NULL, 
            na.rm = TRUE)
        res[, , , , , "lowq"] <- quantile(object, lowq, na.rm = TRUE)
        res[, , , , , "uppq"] <- quantile(object, uppq, na.rm = TRUE)
        return(res)
    })

flq <- FLQuantPoint(stock.n(stks), lowq=0.05, uppq=0.95)


