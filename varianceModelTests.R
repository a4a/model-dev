library(FLa4a)
data(ple4)
data(ple4.index)

stk <- window(ple4, 1990, 2004)
idx <-window(ple4.index, 1990, 2004) 

fmodel <- ~factor(age) + factor(year)
qmodel <- list(~ s(age, k=4))

fit. <- a4a(stock=stk, qmodel = qmodel, fmodel=fmodel, indices=FLIndices(idx), fit ="assessment")
xyplot(data~year|age, data=log(index(window(ple4.index, 1990, 2004))/fit.@index[[1]]), panel=function(x,y,...){panel.abline(h=0, col.line="gray80");panel.xyplot(x,y,...)}) 

varslt <- catch.n(stk)
varslt[] <- 1
catch.n(stk) <- FLQuantDistr(catch.n(stk), varslt)
index.var(idx)[] <- 0.01
fit.. <- a4a(stock=stk, qmodel = qmodel, fmodel=fmodel, indices=FLIndices(idx), fit ="assessment")
xyplot(data~year|age, data=log(index(window(ple4.index, 1990, 2004))/fit..@index[[1]]), panel=function(x,y,...){panel.abline(h=0, col.line="gray80");panel.xyplot(x,y,...)}) 

