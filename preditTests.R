library(FLa4a)
data(ple4)
data(ple4.index)

fmodel <- ~factor(age) + factor(year)
qmodel <- list(~factor(age) + year)

fit. <- a4a(stock=ple4, qmodel = qmodel, fmodel=fmodel, indices=FLIndices(ple4.index), fit ="assessment")

# catchability
rng <- fit.@pars@qmodel[[1]]@range
df0 <- expand.grid(age=rng["min"]:rng["max"], year=rng["minyear"]:rng["maxyear"])
opts <- options(contrasts = c(unordered = "contr.sum", ordered = "contr.poly"))
Xr <- getX(fit.@pars@qmodel[[1]]@Mod, df0)
options(opts)
pars <- c(fit.@pars@qmodel[[1]]@params)
q.fit <- index(fit.)[[1]]
q.fit[] <- exp(Xr %*% pars + fit.@pars@qmodel[[1]]@centering)
xyplot(data~age, groups=year, data=q.fit)



