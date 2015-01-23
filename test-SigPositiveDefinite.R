library(FLa4a)
library(R2admb)

load("sim_fail.Rdata")
fit_try <- a4aSCA(
      stock=test_stk,
    indices = ids,
    fmodel = fmodel_temp,
    qmodel = qmodel_temp,
    srmodel = srmodel_temp,
    fit="assessment", wkdir="test")

vcov(vmodel(pars(fit_try))[[1]])
vcov2 <- read_pars("test-1/a4a")$vcov
vcov2[grep("vpar", dimnames(vcov2)[[2]]), grep("vpar", dimnames(vcov2)[[2]])]

# hessian is not symetric
dim(getADMBHessian("test")$hes)
dim(getADMBCovariance("test")$cov)

all.equal(c(hess[upper.tri(hess)]),c(hess[lower.tri(hess)]))

# what about the cov matrix ?
wkdir <- "test"
filename <- file(paste0(wkdir, "/admodel.cov"), "rb")
on.exit(close(filename))
num.pars <- readBin(filename, "integer", 1)
cov <- readBin(filename, "numeric", num.pars^2)
cov <- matrix(cov, ncol = num.pars, nrow = num.pars)
all.equal(c(cov[upper.tri(cov)]),c(cov[lower.tri(cov)]))

# with ple4
data(ple4)
data(ple4.indices)
wkdir <- "ple4test"
fit <- a4aSCA(ple4, ple4.indices, wkdir=wkdir)

hess <- getADMBHessian(wkdir)$hes
all.equal(c(hess[upper.tri(hess)]),c(hess[lower.tri(hess)]))

filename <- file(paste0(wkdir, "/admodel.cov"), "rb")
on.exit(close(filename))
num.pars <- readBin(filename, "integer", 1)
cov <- readBin(filename, "numeric", num.pars^2)
cov <- matrix(cov, ncol = num.pars, nrow = num.pars)
all.equal(c(cov[upper.tri(cov)]),c(cov[lower.tri(cov)]))

vcov2 <- read_pars("ple4test/a4a")$vcov
vcov2[grep("vpar", dimnames(vcov2)[[2]]), grep("vpar", dimnames(vcov2)[[2]])]
vcov(vmodel(pars(fit))[[1]])




