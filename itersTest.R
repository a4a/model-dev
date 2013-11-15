library(FLa4a)
data(ple4)
data(ple4.indices)

# data
#-------------------------------------------------------
stock <- propagate(ple4, 100)
indices <- FLIndices(bts=propagate(ple4.indices[["BTS-Isis"]], 2))

# checking
#-------------------------------------------------------
dms <- do.call(rbind.data.frame, c(list(catch = c(dims(stock), startf = NA, endf = NA)), lapply(indices, dims)))

# only allow 1 season for surveys
if (any(dms $ season[-1] > 1)) stop("only one season per survey - please split into seperate surveys.")

# the check is failing 
# SORTED OUT
if (!identical(sort(unique(dms $ iter)), sort(unique(c(1L, max(dms $ iter)))))) stop("incosistent number of iterations in stock and indices")

# model
#-------------------------------------------------------

fmodel <- ~factor(age) + factor(year)
qmodel <- list(~factor(age))
rmodel <- ~factor(year)

# stk with iters idx without
fit1 <- a4a(fmodel = fmodel, qmodel = qmodel, srmodel = rmodel, stock=stock, indices=ple4.indices["BTS-Isis"])

# stk without iters idx with iters
fit2 <- a4a(fmodel = fmodel, qmodel = qmodel, srmodel = rmodel, stock=ple4, indices=indices)

# stk and idx with iters
fit3 <- a4a(fmodel = fmodel, qmodel = qmodel, srmodel = rmodel, stock=stock, indices=indices)
