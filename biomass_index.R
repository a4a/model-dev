

# test
library(FLa4a)


data(ple4)
data(ple4.indices)

fit <- a4a(stock = ple4, indices = ple4.indices)

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


mod <- a4aInternal(fmodel = ~ s(age, k = 3) + s(year, k = 10), qmodel = list( ~1), stock = stock, indices = indices)

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




