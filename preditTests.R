library(FLa4a)
data(ple4)
data(ple4.index)

fmodel <- ~factor(age) + factor(year)
qmodel <- list(~s(age, k=4))

fit. <- a4aSCA(stock=ple4, qmodel = qmodel, fmodel=fmodel, indices=FLIndices(ple4.index), fit ="assessment", wkdir="test")

# catchability
rng <- fit.@pars@qmodel[[1]]@range
df0 <- expand.grid(age=rng["min"]:rng["max"], year=rng["minyear"]:rng["maxyear"])
opts <- options(contrasts = c(unordered = "contr.sum", ordered = "contr.poly"))
Xr <- getX(fit.@pars@qmodel[[1]]@Mod, df0)
options(opts)
pars <- c(fit.@pars@qmodel[[1]]@params)
q.fit <- index(fit.)[[1]]
q.fit[] <- exp(Xr %*% pars + fit.@pars@qmodel[[1]]@centering - fit.@pars@stkmodel@centering)
xyplot(data~age, groups=year, data=q.fit)

q.fit2 <- q.fit
rng <- range(ple4.index)
flq <- stock.n(fit.)*exp(-mean(rng[c("startf","endf")])*(m(fit.)+harvest(fit.)))
flq <- flq[ac(rng["min"]:rng["max"]), ac(rng["minyear"]:rng["maxyear"])]
q.fit2[] <- index(fit.)[[1]]/flq

q.res <- log(index(ple4.index)/index(fit.)[[1]])
xyplot(data~year|age, data=log(fit.@index[[1]]/index(ple4.index)), panel=function(x,y,...){panel.abline(h=0, col.line="gray80");panel.xyplot(x,y,...)})

xyplot(data~year|age, data=q.res, panel=function(x,y,...){panel.abline(h=0, col.line="gray80");panel.xyplot(x,y,...)})

