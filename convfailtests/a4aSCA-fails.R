# EXAMPLES OF FITS THAT FAIL
# NOTE THAT SETTING wkdir IN THE CALL WILL CREATE A FOLDER WHERE ALL THE ADMB FILES ARE DROPPED

library(FLa4a)
data(ple4)
data(ple4.indices)
data(ple4.index)
load("tests.rdata")

# plaice in 4
# wkdir argument will drop the ADMB files in your disk
fit1a <- a4aSCA(ple4, FLIndices(ple4.index), fmodel=~factor(age)+ factor(year), qmodel=list(~s(age, k=3)+ s(year, k=20)), wkdir="ple4test")
fit1b <- a4aSCA(ple4, FLIndices(ple4.index), fmodel=~factor(age)+ factor(year), qmodel=list(~s(age, k=3)+ s(year, k=20)), wkdir="ple4test", center=FALSE)
fit1c <- a4aSCA(ple4, FLIndices(ple4.index), fmodel=~factor(age)+ factor(year), qmodel=list(~s(age, k=3)+ s(year, k=20)), wkdir="ple4test", center=1)

# horse mackerel NS
fmodel <- ~ factor(age) + factor(year)
qmodel <- list(~s(year, k=5), ~s(year, k=4))
fit2a <- a4aSCA(hom, hom.idx, fmodel=fmodel, qmodel=qmodel, wkdir="nshomtest")
fit2b <- a4aSCA(hom, hom.idx, fmodel=fmodel, qmodel=qmodel, wkdir="nshomtest", center=FALSE)
fit2c <- a4aSCA(hom, hom.idx, fmodel=fmodel, qmodel=qmodel, wkdir="nshomtest", center=1)

# southern hake
fmodel <- ~factor(age) + factor(year)
qmodel <- list(~s(age, k = 3), ~s(age, k = 3), ~s(age, k = 3))
srmodel <- ~bevholt(CV = 0.1)
fit3a <- a4aSCA(stock=hke, indices = hke.idx, fmodel = fmodel, qmodel = qmodel, srmodel = srmodel, wkdir="shketest")
fit3b <- a4aSCA(stock=hke, indices = hke.idx, fmodel = fmodel, qmodel = qmodel, srmodel = srmodel, wkdir="shketest", center=FALSE)
fit3c <- a4aSCA(stock=hke, indices = hke.idx, fmodel = fmodel, qmodel = qmodel, srmodel = srmodel, wkdir="shketest", center=1)

# sims
fit4 <- a4aSCA(ple4, FLIndices(ple4.index), fmodel=~factor(age)+ factor(year), qmodel=list(~s(age, k=3)), vmodel=list(~1, ~1))
set.seed(1234)
stk <- ple4 + simulate(fit4, 10)
fit4a <- a4aSCA(stk, FLIndices(ple4.index), fmodel=~factor(age)+ factor(year), qmodel=list(~s(age, k=3)), vmodel=list(~1, ~1))
fit4b <- a4aSCA(stk, FLIndices(ple4.index), fmodel=~factor(age)+ factor(year), qmodel=list(~s(age, k=3)), vmodel=list(~1, ~1), center=FALSE)
fit4c <- a4aSCA(stk, FLIndices(ple4.index), fmodel=~factor(age)+ factor(year), qmodel=list(~s(age, k=3)), vmodel=list(~1, ~1), center=1)

# sims

fit5 <- a4aSCA(hke, hke.idx, fmodel=~factor(age)+ factor(year), vmodel=list(~1, ~1, ~1, ~1), n1model=~s(age, k=3))
set.seed(1234)
stk <- hke + simulate(fit5, 10)
fit5a <- a4aSCA(stk, hke.idx, fmodel=~factor(age)+ factor(year), vmodel=list(~1, ~1, ~1, ~1), n1model=~s(age, k=3))
fit5b <- a4aSCA(stk, hke.idx, fmodel=~factor(age)+ factor(year), vmodel=list(~1, ~1, ~1, ~1), n1model=~s(age, k=3), center=FALSE)
fit5c <- a4aSCA(stk, hke.idx, fmodel=~factor(age)+ factor(year), vmodel=list(~1, ~1, ~1, ~1), n1model=~s(age, k=3), center=1)

