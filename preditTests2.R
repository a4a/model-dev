vmodel <- list(~1, ~1)
qmodel=list(~s(age, k=4))

fit <-  a4aSCA(ple4, FLIndices(ple4.index), fit="assessment", vmodel=vmodel, qmodel=qmodel)

qmod <- qmodel(pars(fit))[[1]]
qmod@centering
stkmod <- stkmodel(pars(fit))
qmod@centering <- qmod@centering-stkmod@centering
predict(qmod)

sfrac <- mean(range(ple4.index)[c("startf", "endf")])
Z <- (m(ple4) + harvest(fit2))*sfrac
lst <- dimnames(fit2@index[[1]])
lst$x <- stock.n(fit2)*exp(-Z)
stkn <- do.call("trim", lst)
qhat <- index(fit2)[[1]]/stkn
all.equal(c(qhat), c(predict(qmod2)))


vmodel <- list(~1, ~1)
qmodel=list(~s(age, k=4))
fit2 <- FLa4a:::a4aInternal(stock=ple4, indices=FLIndices(ple4.index), fit="assessment", vmodel=vmodel, qmodel=qmodel, center=FALSE)

qmod2 <- qmodel(pars(fit2))[[1]]
qmod2@centering
stkmod2 <- stkmodel(pars(fit2))
qmod2@centering <- qmod2@centering-stkmod2@centering
predict(qmod2)


all.equal(predict(qmod2), predict(qmod))

sfrac <- mean(range(ple4.index)[c("startf", "endf")])
Z <- (m(ple4) + harvest(fit2))*sfrac
lst <- dimnames(fit2@index[[1]])
lst$x <- stock.n(fit2)*exp(-Z)
stkn <- do.call("trim", lst)
qhat <- index(fit2)[[1]]/stkn
all.equal(c(qhat), c(predict(qmod2)))




vmodel <- list(~1, ~1)
qmodel=list(~s(age, k=4))
fit2 <- FLa4a:::a4aInternal(stock=ple4, indices=FLIndices(ple4.index), fit="assessment", vmodel=vmodel, qmodel=qmodel, center=FALSE)

qmod2 <- qmodel(pars(fit2))[[1]]

sfrac <- mean(range(ple4.index)[c("startf", "endf")])
Z <- (m(ple4) + harvest(fit2))*sfrac
lst <- dimnames(fit2@index[[1]])
lst$x <- stock.n(fit2)*exp(-Z)
stkn <- do.call("trim", lst)
qhat <- index(fit2)[[1]]/stkn
all.equal(c(qhat), c(predict(qmod2)))

