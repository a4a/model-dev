library(FLa4a)
data(ple4)
data(ple4.index)
biofull <- 0.001*stock(ple4)
biofull <- FLIndexBiomass(index=biofull)
range(biofull)[c("startf","endf")] <- c(0,0)

fit0 <- sca(ple4, FLIndices(a=ple4.index, b=biofull), qmodel=list(~s(age, k=4), ~1))
fit1 <- sca(ple4, FLIndices(b=biofull, a=ple4.index), qmodel=list(~1, ~s(age, k=4)))


bio24 <- 0.001*quantSums(stock.n(ple4)[2:4]*stock.wt(ple4)[2:4])
bio24 <- FLIndexBiomass(index=bio24)
range(bio24)[c("startf","endf")] <- c(0,0)
range(bio24)[c("min","max")] <- c(2,4)

fit1 <- sca(ple4, FLIndices(ple4.index), fit="assessment")
fit2 <- sca(ple4, FLIndices(biofull), fit="assessment")
fit3 <- sca(ple4, FLIndices(bio24), fit="assessment")
fit4 <- sca(ple4, FLIndices(ple4.index, biofull, bio24), qmodel=list(~s(age, k=4), ~1, ~1))

plot(FLStocks(orig=ple4, f1=ple4+fit1, f2=ple4+fit2, f3=ple4+fit3, f4=ple4+fit4), key=TRUE)

predict(fit3)




