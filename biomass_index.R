

# test
library(FLa4a)


data(ple4)
data(ple4.indices)

fit <- sca(stock = ple4, indices = ple4.indices)

plot(ple4 + fit)



## testing

#library(devtools)
#install_github("a4a", "colinpmillar", subdir = "packages/FLa4a")
#system("chmod a+x /home/millaco/R/x86_64-pc-linux-gnu-library/3.0/FLa4a/bin/linux/a4a")

options(width=140)

library(FLa4a)

data(ple4)

stock <- ple4
harvest.spwn(stock) <- .5
m.spwn(stock) <- .5
indices <- FLIndices(bio = FLIndex(index = exp(-10) * tsb(stock), name = "bio"))
range(indices[[1]])[c("startf","endf")] <- 0.5


mod <- sca(fmodel = ~ s(age, k = 3) + s(year, k = 10), qmodel = list( ~1), stock = stock, indices = indices)

yearTotals(exp(index(mod)[[1]]) * stock.wt(stock))[1,]

index(indices[[1]])




fmodel  = ~ s(age, k = 3) + factor(year)
qmodel  = lapply(seq(length(indices)), function(i) ~ 1)
srmodel = ~ factor(year)
n1model = ~ factor(age)
vmodel  = lapply(seq(length(indices) + 1), function(i) ~ 1)

covar = NULL

verbose = FALSE
fit = "assessment"
niters = 1000
wkdir = NULL


mod <- a4a(fmodel = ~ s(age, k = 3) + s(year, k = 10), stock = stock, indices = indices)

plot(log(c(index(mod)[[1]] / index(indices[[1]]))))



library(FLa4a)
data(ple4)

bioidx <- FLIndex(FLQuant(NA, dimnames=list(age="all", year=range(ple4)["minyear"]:range(ple4)["maxyear"])), name="bioidx")
index(bioidx) <- stock(ple4)*0.001
index(bioidx) <- index(bioidx)*exp(rnorm(index(bioidx), sd=0.01))
range(bioidx)[c("startf","endf")] <- c(0,0)
fit <- sca(ple4, FLIndices(bioidx), fmodel=~factor(age) + s(year, k=25), qmodel=list(~1), fit="assessment")
plot(residuals(fit, ple4, FLIndices(bioidx)))
retro <- ra(ple4, FLIndices(bioidx), 5, fmodel=~factor(age) + s(year, k=25), qmodel=list(~1))
plot(retro)
sim <- simulate(fit, 1000, stock=ple4)
plot(ple4+sim)




# see the difference in parameters precision
vcov(qmodel(pars(fit))[[1]])
diag(vcov(stkmodel(pars(fit)))[,,1])

varslt <- catch.n(ple4)
varslt[] <- 0.4
catch.n(ple4) <- FLQuantDistr(catch.n(ple4), varslt)
# variance of observed indices
varslt <- index(bioidx)
varslt[] <- 0.1
index.var(bioidx) <- varslt

fit <- sca(ple4, FLIndices(bioidx), qmodel=list(~1), fit="assessment")

# see the difference in parameters precision
vcov(qmodel(pars(fit))[[1]])
diag(vcov(stkmodel(pars(fit)))[,,1])

varslt <- catch.n(ple4)
varslt[] <- 0.4
catch.n(ple4) <- FLQuantDistr(catch.n(ple4), varslt)
# variance of observed indices
varslt <- index(bioidx)
varslt[] <- 0.01
index.var(bioidx) <- varslt

fit <- sca(ple4, FLIndices(bioidx), fmodel=~s(age, k=15) + s(year, k=5), qmodel=list(~s(year, k=5)), fit="assessment")


fit <- a4aSCA(ple4, FLIndices(bioidx), qmodel=list(~s(year, k=5)), fmodel=~s(age, k=4) + s(year, k=5), vmodel=list(~1, ~1), fit="assessment")
data(ple4) # this is just to bypass a missing feature I just found
plot(residuals(fit, ple4, FLIndices(bioidx)))


fit <- a4aSCA(ple4, FLIndices(bioidx), qmodel=list(~factor(replace(year, year>1990, 1990))), fmodel=~s(age, k=4) + s(year, k=25), vmodel=list(~1, ~1), fit="assessment")
data(ple4) # this is just to bypass a missing feature I just found
plot(residuals(fit, ple4, FLIndices(bioidx)))


