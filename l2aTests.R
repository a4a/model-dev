mm <- matrix(NA, ncol=3, nrow=3)
diag(mm) <- c(50, 0.001,0.001)
mm[upper.tri(mm)] <- mm[lower.tri(mm)] <- c(0.1,0.01,0.00004)
vbObj <- a4aGr(grMod=~linf*(1-exp(-k*(t-t0))), grInvMod=~t0-1/k*log(1-len/linf), params=FLPar(linf=58.5, k=0.086, t0=0.001, units=c("cm","ano-1","ano")), vcov=mm, distr="norm")
data(rfLen)
cth <- catch.n(rfLen.stk) 
# both with iter=1
cthA1 <- l2a(cth, vbObj)
# both with iter=n
cthA2 <- l2a(propagate(cth,10), mvrnorm(10, vbObj))
# mod: iter=n, data: iter=1
cthA3 <- l2a(cth, mvrnorm(10, vbObj))
# mod: iter=1, data: iter=n
cthA4 <- l2a(propagate(cth,10), vbObj)

# converting a stock object
stk <- l2a(rfLen.stk, vbObj)
stkn <- stock.n(stk)
cthn <- catch.n(stk)

# aggregate areas, seasons
dnm <- dimnames(stkn)
gr <- expand.grid(dnm[4:5])

stknt <- stkn[,,,gr[1,"season"],gr[1,"area"],]
cthnt <- cthn[,,,gr[1,"season"],gr[1,"area"],]

for(i in 2:nrow(gr)){
	stknt <- stknt + stkn[,,,gr[i,"season"],gr[i,"area"],]
	cthnt <- cthnt + cthn[,,,gr[i,"season"],gr[i,"area"],]
}

dimnames(stknt)[4:5] <- dimnames(cthnt)[4:5] <- dimnames(FLQuant())[4:5]

# compute Z
dm <- dim(stknt)
am <- dm[1]
ym <- dm[2]
ni1j1 <- stknt[-1,-1]
nij <- stknt[-am,-ym]
zij <- log(nij/ni1j1)

mij <- m(stk)[-am,-ym]
fij <- zij-mij




cthn <- catch.n(stk)









cthn <- harvest(stk)/z(stk)*stock.n(stk)*(1-exp(-z(stk)))
sum(round(cthn/catch.n(stk), 3)>1)


rfAge.stk <- l2a(rfLen.stk, mvrnorm(10, vbObj))
rfAge.stk <- l2a(propagate(rfLen.stk, 10), vbObj)




cthn <- l2a(catch.n(rfLen.stk), vbObj)
stkn <- l2a(stock.n(rfLen.stk), vbObj)

mrat <- l2a(m(rfLen.stk)/(m(rfLen.stk)+rfLen.stk@harvest), vbObj, "mean")

stkn <- stock.n(stk)

