library(FLa4a)
library(parallel)
data(ple4)
data(ple4.index)

fmodel <- ~factor(age) + factor(year)
qmodel <- list(~factor(age) + year)

fit. <- a4a(stock=window(ple4, 1990, 2004), qmodel = qmodel, fmodel=fmodel, indices=FLIndices(window(ple4.index, 1990, 2004)), fit ="assessment", wkdir="test")

attach("../stksSim.rdata")
options(mc.cores=4)
fmodel <- ~factor(age) + factor(year)
qmodel <- list(~factor(age))

fits01 <- mclapply(stks01, function(x) try(a4a(stock=x$stock, qmodel = qmodel, fmodel=fmodel, indices=FLIndices(x$index), fit ="assessment")))
fits02 <- mclapply(stks02, function(x) try(a4a(stock=x$stock, qmodel = qmodel, fmodel=fmodel, indices=FLIndices(x$index), fit ="assessment")))
fits03 <- mclapply(stks03, function(x) try(a4a(stock=x$stock, qmodel = qmodel, fmodel=fmodel, indices=FLIndices(x$index), fit ="assessment")))
fits04 <- mclapply(stks04, function(x) try(a4a(stock=x$stock, qmodel = qmodel, fmodel=fmodel, indices=FLIndices(x$index), fit ="assessment")))
fits05 <- mclapply(stks05, function(x) try(a4a(stock=x$stock, qmodel = qmodel, fmodel=fmodel, indices=FLIndices(x$index), fit ="assessment")))

save(fits01, fits02, fits03, fits04, fits05, file="fits00")
rm(fits01, fits02, fits03, fits04, fits05)

attach("fits00NoPhase")

nloglNoPhase <- c(
	lapply(fits01, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits02, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits03, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits04, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits05, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",])
)

attach("fits00")

nloglPhase <- c(
	lapply(fits01, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits02, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits03, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits04, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits05, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",])
)


nloglCompare1 <- data.frame(stk=names(stks01), nlogl=c(unlist(nloglNoPhase), unlist(nloglPhase)), phase=rep(c(TRUE, FALSE), rep(length(stks01)*5,2)), expl=rep(c(1:5), rep(length(stks01), 5)))

nloglCompare2 <- data.frame(stk=names(stks01), nlogl=c(unlist(nloglNoPhase), unlist(nloglPhase)), phase=rep(c(TRUE, FALSE), rep(length(stks01)*5,2)), expl=rep(c(1:5), rep(length(stks01), 5)))

data.frame(no=unlist(nloglNoPhase), yes2=unlist(nloglPhase), yes1=subset(nloglCompare, phase==TRUE)$nlogl)

