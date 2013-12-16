library(FLa4a)
data(ple4)
data(ple4.index)

fmodel <- ~factor(age) + factor(year)
qmodel <- list(~factor(age) + year)

fit. <- a4a(stock=window(ple4, 1990, 2004), qmodel = qmodel, fmodel=fmodel, indices=FLIndices(window(ple4.index, 1990, 2004)), fit ="assessment", wkdir="test")


