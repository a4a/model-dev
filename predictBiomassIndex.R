library(FLa4a)
data(ple4)
data(ple4.index)

# create a biomass index that reflects the full age range
# note that in this case one doesn't need to set the age range

biofull <- 0.001*stock(ple4)
biofull <- FLIndexBiomass(index=biofull)
range(biofull)[c("startf","endf")] <- c(0,0)

# now create a biomass index that reflects the biomass of ages 2-4
# note that in this case one sets the age range to c(2,4)

bio24 <- 0.001*quantSums(stock.n(ple4)[2:4]*stock.wt(ple4)[2:4])
bio24 <- FLIndexBiomass(index=bio24)
range(bio24)[c("startf","endf")] <- c(0,0)
range(bio24)[c("min","max")] <- c(2,4)

# the fits
# note that the qmodel must be defined for the year effect only

fit1 <- sca(ple4, FLIndices(ple4.index))
fit2 <- sca(ple4, FLIndices(biofull), qmodel=list(~1))
fit3 <- sca(ple4, FLIndices(bio24), qmodel=list(~1), fit="assessment")
fit4 <- sca(ple4, FLIndices(ple4.index, biofull, bio24), qmodel=list(~s(age, k=4), ~1, ~1))

# everybody likes a plot

plot(FLStocks(orig=ple4, f1=ple4+fit1, f2=ple4+fit2, f3=ple4+fit3, f4=ple4+fit4), key=TRUE) 


obj <- predict(fit3)
obj$qmodel

