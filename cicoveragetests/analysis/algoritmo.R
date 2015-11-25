library(FLa4a)
load("EJ_stock.Rdata")

iter <- 2

stk_pg <- setPlusGroup(iter(shake_gis_tri,iter), pg_tri[iter])
age_smoother <- max(3,floor(dim(stock.n(stk_pg))[1]/2)+1) # has to be at least 3
fmodel_temp <-as.formula(paste("~te(age, year, k = c(",age_smoother,",15))",sep=""))

idxs <- FLIndices(lapply(indices_tri, "iter", iter))

fit <- sca(stk_pg, idxs, fmodel=fmodel_temp, qmodel=qmodel3, srmodel=rmodel1, fit="assessment")

stk <- stk_pg + simulate(fit, 2)

refit <- sca(iter(stk, 2), idxs, fmodel=fmodel_temp, qmodel=qmodel3, srmodel=rmodel1, fit="assessment")


