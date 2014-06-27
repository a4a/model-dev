library(FLa4a)
data(ple4)
data(ple4.indices)
flqs <- FLQuants(catch.n=catch.n(fit), stock.n=stock.n(fit), f=harvest(fit))
wireframe(data ~ age + year| qname, data = as.data.frame(flqs), drape = TRUE, screen = list(x = -90, y=-45))


