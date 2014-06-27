library(FLa4a)
library(parallel)
data(ple4)
data(ple4.index)

fmodel <- ~factor(age) + factor(year)
qmodel <- list(~factor(age) + year)

fit. <- a4a(stock=window(ple4, 1990, 2004), qmodel = qmodel, fmodel=fmodel, indices=FLIndices(window(ple4.index, 1990, 2004)), fit ="assessment", wkdir="test")


# these fits were ran with distinct a4a executables
attach("../stksSim.rdata")
options(mc.cores=18)
fmodel <- ~ factor(age) + factor(year)
qmodel <- list(~ factor(age))

fits01 <- mclapply(stks01, function(x) try(a4a(stock=x$stock, qmodel = qmodel, fmodel=fmodel, indices=FLIndices(x$index), fit ="assessment")))
fits02 <- mclapply(stks02, function(x) try(a4a(stock=x$stock, qmodel = qmodel, fmodel=fmodel, indices=FLIndices(x$index), fit ="assessment")))
fits03 <- mclapply(stks03, function(x) try(a4a(stock=x$stock, qmodel = qmodel, fmodel=fmodel, indices=FLIndices(x$index), fit ="assessment")))
fits04 <- mclapply(stks04, function(x) try(a4a(stock=x$stock, qmodel = qmodel, fmodel=fmodel, indices=FLIndices(x$index), fit ="assessment")))
fits05 <- mclapply(stks05, function(x) try(a4a(stock=x$stock, qmodel = qmodel, fmodel=fmodel, indices=FLIndices(x$index), fit ="assessment")))

save(fits01, fits02, fits03, fits04, fits05, file="fits00f2r3")
rm(fits01, fits02, fits03, fits04, fits05)

# analysis

attach("fits00nophase")

nloglNoPhase <- c(
	lapply(fits01, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits02, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits03, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits04, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits05, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",])
)
detach()

attach("fits00f2")

nloglf2 <- c(
	lapply(fits01, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits02, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits03, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits04, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits05, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",])
)
detach()

attach("fits00q2")

nloglq2 <- c(
	lapply(fits01, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits02, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits03, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits04, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits05, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",])
)
detach()

attach("fits00Ny2")

nloglNy2 <- c(
	lapply(fits01, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits02, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits03, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits04, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits05, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",])
)
detach()

attach("fits00v2")

nloglv2 <- c(
	lapply(fits01, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits02, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits03, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits04, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits05, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",])
)
detach()

attach("fits00R2")

nloglR2 <- c(
	lapply(fits01, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits02, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits03, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits04, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits05, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",])
)
detach()

attach("fits00f2q2")

nloglf2q2 <- c(
	lapply(fits01, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits02, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits03, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits04, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits05, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",])
)
detach()

attach("fits00f2v2")

nloglf2v2 <- c(
	lapply(fits01, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits02, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits03, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits04, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits05, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",])
)
detach()

attach("fits00f2R2")

nloglf2R2 <- c(
	lapply(fits01, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits02, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits03, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits04, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits05, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",])
)
detach()

attach("fits00f2Ny2")

nloglf2Ny2 <- c(
	lapply(fits01, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits02, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits03, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits04, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits05, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",])
)
detach()

attach("fits00v1")

nloglv1 <- c(
	lapply(fits01, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits02, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits03, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits04, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits05, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",])
)
detach()

attach("fits00v1f3")

nloglv1f3 <- c(
	lapply(fits01, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits02, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits03, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits04, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits05, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",])
)
detach()

attach("fits00f2r3")

nloglf2r3 <- c(
	lapply(fits01, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits02, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits03, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits04, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",]),
	lapply(fits05, function(x) if(is(x, "try-error")) return(NA) else fitSumm(x)["nlogl",])
)
detach()


nloglCompare <- data.frame(stk=names(stks01), noPhase=unlist(nloglNoPhase), f2=unlist(nloglf2), q2=unlist(nloglq2),  R2=unlist(nloglR2), Ny2=unlist(nloglNy2), v2=unlist(nloglv2), f2q2=unlist(nloglf2q2), f2v2=unlist(nloglf2v2), f2R2=unlist(nloglf2R2), f2Ny2=unlist(nloglf2Ny2), v1=unlist(nloglv1), v1f3=unlist(nloglv1f3), f2r3=unlist(nloglf2r3), expl=rep(c(1:5), rep(length(stks01), 5)))

compareWithNophase <- rbind(
	ok=apply(is.na(subset(nloglCompare, !is.na(noPhase))), 2, sum)[-c(1,ncol(nloglCompare))],
	fail=apply(is.na(subset(nloglCompare, is.na(noPhase))), 2, sum)[-c(1,ncol(nloglCompare))], 
	likhigher=apply(nloglCompare[,-1]>nloglCompare[,rep(2,ncol(nloglCompare)-1)], 2, sum, na.rm=T)[-(ncol(nloglCompare)-1)], 
	liklower=apply(nloglCompare[,-1]<nloglCompare[,rep(2,ncol(nloglCompare)-1)], 2, sum, na.rm=T)[-(ncol(nloglCompare)-1)]
)	








