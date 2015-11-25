##################################################################################
##################################################################################

# code for implementing a maximum likelihood version of the CBBM in R
# leire ibaibarriaga, 29/01/2013

##################################################################################
##################################################################################

# the final Bayesian run in the new stock annex is in "C:/use/ices/wghansa/wghansa2013/mio/cbbm_newdata/runs_varestim"

##################################################################################
##################################################################################

# read data file

bio.dat <- read.table("realdata_87-13_jun2013.txt", header=T)

# compute age 1 biomass proportion for depm, acoustics and commercial catches

bio.dat$BP.depm <- bio.dat$Bage1.depm/bio.dat$Btot.depm
bio.dat$BP.ac <- bio.dat$Bage1.ac/bio.dat$Btot.ac
bio.dat$CPsem1 <- bio.dat$Cage1sem1/bio.dat$Ctotsem1
bio.dat$CPsem2 <- bio.dat$Cage1sem2/bio.dat$Ctotsem2

bio.dat$CPsem1[is.nan(bio.dat$CPsem1)] <- NA
bio.dat$CPsem2[is.nan(bio.dat$CPsem2)] <- NA

# number of years

Y <- dim(bio.dat)[1]

# when there are no catches we do not have an observation for the age1 proportion

bio.dat$Cage1sem1[bio.dat$Ctotsem1==0] <- NA
bio.dat$Cage1sem2[bio.dat$Ctotsem2==0] <- NA

################################################################################
################################################################################

source("functions_loglik.r")

################################################################################
################################################################################

# general initial values

logqdepm.ini1 <- 0
logqac.ini1 <- 0
logqrobs.ini1 <- -4
logkrobs.ini1 <- 0.5
logpsidepm.ini1 <- log(7)
logpsiac.ini1 <- log(7)
logpsirobs.ini1 <- log(5)
xidepm.ini1 <- 4
xiac.ini1 <- 4
xicatch.ini1 <- 4
logB0.ini1 <- 10
logR.ini1 <- rnorm(Y, 11, 1) # because when all the values are the same cf() can give problems 
logsage1sem1.ini1 <- log(0.5)
logsage1sem2.ini1 <- log(1.5)
logfsem1.ini1 <- rep(log(0.5), Y)
logfsem2.ini1 <- rep(log(0.5), Y)
logG1.ini1 <- 0.5
logG2.ini1 <- 0.5
logpsig.ini1 <- log(30)

# initials for loglik.f 

ini1 <- c(logqdepm.ini1,
         logqac.ini1,
         logqrobs.ini1,
         logkrobs.ini1,
         logpsidepm.ini1,
         logpsiac.ini1,
         logpsirobs.ini1,
         xidepm.ini1,
         xiac.ini1,
         xicatch.ini1,
         logB0.ini1,
         logR.ini1,
         logsage1sem1.ini1,
         logsage1sem2.ini1,
         logfsem1.ini1,
         logfsem2.ini1,
         logG1.ini1,
         logG2.ini1,
         logpsig.ini1)

# initials for loglikcf.f (closed form)

inicf1 <- c(logpsidepm.ini1,
          logpsiac.ini1,
          xidepm.ini1,
          xiac.ini1,
          xicatch.ini1,
          logB0.ini1,
          logR.ini1,
          logsage1sem1.ini1,
          logsage1sem2.ini1,
          logfsem1.ini1,
          logfsem2.ini1,
          logG1.ini1,
          logG2.ini1)

# initial values from the posterior median from Bayesian CBBM

logqdepm.ini2 <- -0.47095
logqac.ini2 <- 0.27295
logqrobs.ini2 <- -3.8615
logkrobs.ini2 <- 0.40575
logpsidepm.ini2 <- log(6.963)
logpsiac.ini2 <- log(7.5305)
logpsirobs.ini2 <- log(4.5655)
xidepm.ini2 <- 4.0895
xiac.ini2 <- 3.5875
xicatch.ini2 <- 2.825
logB0.ini2 <- 9.983
logR.ini2 <- c(9.69, 10.38, 9.15, 11.14, 10.05, 11.43, 11.08, 10.67, 10.81, 10.89,
                   10.78, 11.44, 10.68, 11.41, 11.22, 9.475, 9.885, 10.32, 8.28, 9.85,
                   10.03, 9.12, 9.23, 10.73, 11.52, 10.57, 10.35)
logsage1sem1.ini2 <- log(0.48245)
logsage1sem2.ini2 <- log(1.322)
logfsem1.ini2 <- log(c(1.19, 0.98, 0.91, 1.18, 1.11, 1.11,0.81, 1.09, 1.36, 1.08,
                    0.53, 0.41, 0.51, 0.70, 0.65, 0.55, 0.38, 0.84, 0.16, 0.22,
                    0.01, 0.001, 0.001, 0.41, 0.32, 0.22, 0.33))
logfsem2.ini2 <- log(c(0.31, 0.31, 0.16, 0.58, 0.24, 0.29, 0.47, 0.5, 0.27, 0.52,
                    0.39, 0.39, 0.36, 0.33, 0.43, 0.46, 0.56, 0.52, 0.001, 0.008,
                    0.001, 0.001, 0.001, 0.16, 0.06, 0.16, 0.001))
logG1.ini2 <- -0.61
logG2.ini2 <- -1.42
logpsig.ini2 <- log(26.06)


# initials for loglik.f

ini2 <- c(logqdepm.ini2,
         logqac.ini2,
         logqrobs.ini2,
         logkrobs.ini2,
         logpsidepm.ini2,
         logpsiac.ini2,
         logpsirobs.ini2,
         xidepm.ini2,
         xiac.ini2,
         xicatch.ini2,
         logB0.ini2,
         logR.ini2,
         logsage1sem1.ini2,
         logsage1sem2.ini2,
         logfsem1.ini2,
         logfsem2.ini2,
         logG1.ini2,
         logG2.ini2,
         logpsig.ini2)

# initials for loglikcf.f (closed form)

inicf2 <- c(logpsidepm.ini2,
          logpsiac.ini2,
          xidepm.ini2,
          xiac.ini2,
          xicatch.ini2,
          logB0.ini2,
          logR.ini2,
          logsage1sem1.ini2,
          logsage1sem2.ini2,
          logfsem1.ini2,
          logfsem2.ini2,
          logG1.ini2,
          logG2.ini2)

inicf3 <- c(logpsidepm.ini2,
            logpsiac.ini2,
            xidepm.ini2,
            xiac.ini2,
            xicatch.ini2,
            logB0.ini2,
            logR.ini1,
            logsage1sem1.ini2,
            logsage1sem2.ini2,
            logfsem1.ini1,
            logfsem2.ini1,
            logG1.ini2,
            logG2.ini2)


################################################################################
################################################################################

# try different optimization methods for the closed form

print(Sys.time())
cf1.nelder <- optim(inicf1, loglikcf.f, bio.dat=bio.dat, method="Nelder-Mead", control=list(fnscale=-1, maxit=50000, reltol=1e-8))
print(Sys.time()) # warnings & not converged

print(Sys.time())
cf2.nelder <- optim(inicf2, loglikcf.f, bio.dat=bio.dat, method="Nelder-Mead", control=list(fnscale=-1, maxit=50000, reltol=1e-8))
print(Sys.time()) # warnings & not converged 

print(Sys.time())
cf3.nelder <- optim(inicf3, loglikcf.f, bio.dat=bio.dat, method="Nelder-Mead", control=list(fnscale=-1, maxit=50000, reltol=1e-8))
print(Sys.time()) # warnings & not converged

# HASTA AQUI DAN Mensajes de aviso perdidos, e.g. 
# In log(rest[match("krobs", names(rest))]) : Se han producido NaNs

print(Sys.time())
cf1.bfgs <- optim(inicf1, loglikcf.f, bio.dat=bio.dat, method="BFGS", control=list(fnscale=-1, maxit=50000, reltol=1e-8))
print(Sys.time()) # 

print(Sys.time())
cf2.bfgs <- optim(inicf2, loglikcf.f, bio.dat=bio.dat, method="BFGS", control=list(fnscale=-1, maxit=50000, reltol=1e-8))
print(Sys.time()) # warnings 

print(Sys.time())
cf3.bfgs <- optim(inicf3, loglikcf.f, bio.dat=bio.dat, method="BFGS", control=list(fnscale=-1, maxit=50000, reltol=1e-8))
print(Sys.time()) # warnings & Error en optim(inicf3, loglikcf.f, bio.dat = bio.dat, method = "BFGS",: non-finite finite-difference value [25] 

# HASTA AQUI warnings: 
#  1: In dbeta(x, shape1, shape2, log) : Se han producido NaNs
# tardaba aprox 5 min

print(Sys.time())
cf1.cg <- optim(inicf1, loglikcf.f, bio.dat=bio.dat, method="CG", control=list(fnscale=-1, maxit=500, reltol=1e-8))
print(Sys.time()) # around 12 min & warnings & not converged

print(Sys.time())
cf2.cg <- optim(inicf2, loglikcf.f, bio.dat=bio.dat, method="CG", control=list(fnscale=-1, maxit=500, reltol=1e-8))
print(Sys.time()) # around 12 min, not converged & warnings: In dbeta(x, shape1, shape2, log) : Se han producido NaNs  

print(Sys.time())
cf3.cg <- optim(inicf3, loglikcf.f, bio.dat=bio.dat, method="CG", control=list(fnscale=-1, maxit=500, reltol=1e-8))
print(Sys.time()) # around 12 min, converged & warnings


print(Sys.time())
cf1.lbfgsb <- optim(inicf1, loglikcf.f, bio.dat=bio.dat, method="L-BFGS-B", control=list(fnscale=-1, maxit=50000))
print(Sys.time()) # 3 min & converged 

print(Sys.time())
cf2.lbfgsb <- optim(inicf2, loglikcf.f, bio.dat=bio.dat, method="L-BFGS-B", control=list(fnscale=-1, maxit=50000))
print(Sys.time()) # 5 min & converged 

print(Sys.time())
cf3.lbfgsb <- optim(inicf3, loglikcf.f, bio.dat=bio.dat, method="L-BFGS-B", control=list(fnscale=-1, maxit=50000))
print(Sys.time()) # warnings & Error en optim(inicf3, loglikcf.f, bio.dat = bio.dat, method = "L-BFGS-B",  : L-BFGS-B necesita valores finitos de 'fn'

print(Sys.time())
cf1.sann <- optim(inicf1, loglikcf.f, bio.dat=bio.dat, method="SANN", control=list(fnscale=-1, maxit=50000, reltol=1e-8))
print(Sys.time()) # 7 min & converged & warnings  In log(rest[match("krobs", names(rest))]) : Se han producido NaNs

print(Sys.time())
cf2.sann <- optim(inicf2, loglikcf.f, bio.dat=bio.dat, method="SANN", control=list(fnscale=-1, maxit=50000, reltol=1e-8))
print(Sys.time()) #  7 min & converged 

print(Sys.time())
cf3.sann <- optim(inicf3, loglikcf.f, bio.dat=bio.dat, method="SANN", control=list(fnscale=-1, maxit=50000, reltol=1e-8))
print(Sys.time()) # 7 min & converged & Warnings: In log(rest[match("krobs", names(rest))]) : Se han producido NaNs

##################################################################################
##################################################################################


##################################################################################
##################################################################################

# function to simulate data

sim.f <- function(param, bio.dat){

  logqdepm <- param[1]
  logqac <- param[2]
  logqrobs <- param[3]
  logkrobs <- param[4]
  logpsidepm <- param[5]
  logpsiac <- param[6]
  logpsirobs <- param[7]
  xidepm <- param[8]
  xiac <- param[9]
  xicatch <- param[10]
  logB0 <- param[11]
  logR <- param[(11+1):(11+Y)]
  logsage1sem1 <- param[11+Y+1]
  logsage1sem2 <- param[11+Y+2]
  logfsem1 <- param[(13+Y+1):(13+Y+Y)]
  logfsem2 <- param[(13+2*Y+1):(13+2*Y+Y)]
  logG1 <- param[13+3*Y+1]
  logG2 <- param[13+3*Y+2]
  logpsig <- param[13+3*Y+3]

  dsurv <- 0.375
  dsem1 <- 0.5
  dsem2 <- 0.5
  M1 <- 0.8
  M2 <- 1.2

  qdepm <- exp(logqdepm)
  qac <- exp(logqac)
  psidepm <- exp(logpsidepm)
  psiac <- exp(logpsiac)
  sage1sem1 <- exp(logsage1sem1)
  sage1sem2 <- exp(logsage1sem2)
  B0 <- exp(logB0)
  R <- exp(logR)
  G1 <- exp(logG1)
  G2 <- exp(logG2)
  fsem1 <- exp(logfsem1)
  fsem2 <- exp(logfsem2)

  fsem1[bio.dat$Ctotsem1==0] <- 0.001 # when there are no catches we set f to zero
  fsem2[bio.dat$Ctotsem2==0] <- 0.001
  logfsem1 <- log(fsem1)
  logfsem2 <- log(fsem2)

# compute population based on the function "C:/use/tesis/analysis/chapter2/model winbugs blackbox/function_calcpop.r"

  mat <- c(fsem1, fsem2, logB0, logR, sage1sem1, sage1sem2, log(M1), log(M2), log(G1), log(G2))
  names(mat) <- c(paste("fsem1[",1:Y,"]",sep=""), paste("fsem2[",1:Y,"]",sep=""), "logB0", paste("logR[",1:Y,"]",sep=""),
                "sage1[1]", "sage1[2]", "logM1", "logM2", "logG1", "logG2")

  pop <- calc.pop(mat, Y=Y, dsurv=dsurv, d=c(dsem1, dsem2))

  B1 <- pop[match(paste("B1[",1:Y,"]",sep=""), names(pop))]
  Btot <- pop[match(paste("Btot[",1:Y,"]",sep=""), names(pop))]
  C1sem1 <- pop[match(paste("C1sem1[",1:Y,"]",sep=""), names(pop))]
  Ctotsem1 <- pop[match(paste("Ctotsem1[",1:Y,"]",sep=""), names(pop))]
  C1sem2 <- pop[match(paste("C1sem2[",1:Y,"]",sep=""), names(pop))]
  Ctotsem2 <- pop[match(paste("Ctotsem2[",1:Y,"]",sep=""), names(pop))]

# create index variables for observation equations (avoiding missing data)

  year.BP.depm <- (1:Y)[!is.na(bio.dat$BP.depm)]
  year.Btot.depm <- (1:Y)[!is.na(bio.dat$Btot.depm)]
  year.BP.ac <- (1:Y)[!is.na(bio.dat$BP.ac)]
  year.Btot.ac <- (1:Y)[!is.na(bio.dat$Btot.ac)]
  year.Robs <- (1:Y)[!is.na(bio.dat$Robs)]
  year.Ctotsem1 <- (1:Y)[!is.na(bio.dat$Ctotsem1) & bio.dat$Ctotsem1>0]
  year.Ctotsem2 <- (1:Y)[!is.na(bio.dat$Ctotsem2) & bio.dat$Ctotsem2>0]
  year.CPsem1 <- (1:Y)[!is.na(bio.dat$CPsem1)]
  year.CPsem2 <- (1:Y)[!is.na(bio.dat$CPsem2)]
  year.G1 <- (1:Y)[!is.na(bio.dat$G1)] 
  year.G2 <- (1:Y)[!is.na(bio.dat$G2)]

  sim.dat <- bio.dat

  idx <- year.BP.depm
  sh1 <- exp(xidepm)*B1[idx]/Btot[idx]
  sh2 <- exp(xidepm)*(1-B1[idx]/Btot[idx])  
  sim.dat$BP.depm[idx] <- rbeta(length(idx), shape1=sh1, shape2=sh2)

  idx <- year.Btot.depm
  ml <- logqdepm + log(Btot[idx])
  sdl <- sqrt((psidepm+bio.dat$psidepm[idx])/(psidepm*bio.dat$psidepm[idx]))
  sim.dat$Btot.depm[idx] <- rlnorm(length(idx), meanlog=ml, sdlog=sdl)

  idx <- year.BP.ac
  sh1 <- exp(xiac)*B1[idx]/Btot[idx]
  sh2 <- exp(xiac)*(1-B1[idx]/Btot[idx])  
  sim.dat$BP.ac[idx] <- rbeta(length(idx), shape1=sh1, shape2=sh2)

  idx <- year.Btot.ac
  ml <- logqac + log(Btot[idx])
  sdl <- sqrt((psiac+bio.dat$psiac[idx])/(psiac*bio.dat$psiac[idx]))
  sim.dat$Btot.ac[idx] <- rlnorm(length(idx), meanlog=ml, sdlog=sdl)

  idx <- year.CPsem1
  sh1 <- exp(xicatch)*C1sem1[idx]/Ctotsem1[idx]
  sh2 <- exp(xicatch)*(1-C1sem1[idx]/Ctotsem1[idx])  
  sim.dat$CPsem1[idx] <- rbeta(length(idx), shape1=sh1, shape2=sh2)

  
  idx <- year.Ctotsem1
  ml <- log(Ctotsem1[idx])
  sdl <- 1/sqrt(400)
  bio.dat$Ctotsem1[idx] <- rlnorm(length(idx), meanlog=ml, sdlog=sdl)

  idx <- year.CPsem2
  sh1 <- exp(xicatch)*C1sem2[idx]/Ctotsem2[idx]
  sh2 <- exp(xicatch)*(1-C1sem2[idx]/Ctotsem2[idx])  
  sim.dat$CPsem2[idx] <- rbeta(length(idx), shape1=sh1, shape2=sh2)

  idx <- year.Ctotsem2
  ml <- log(Ctotsem2[idx])
  sdl <- 1/sqrt(400)
  sim.dat$Ctotsem2[idx] <- rlnorm(length(idx), meanlog=ml, sdlog=sdl)

  idx <- year.Robs
  ml <- logqrobs + exp(logkrobs)*log(R[idx])
  sdl <- 1/sqrt(psirobs)
  sim.dat$Robs[idx] <- rlnorm(length(idx), meanlog=ml, sdlog=sdl)

  idx <- year.G1
  m <- G1
  sd <- 1/sqrt(psig)
  sim.dat$G1[idx] <- rnorm(length(idx), mean=m, sd=sd)

  idx <- year.G2
  m <- G2
  sd <- 1/sqrt(psig)
  sim.dat$G2[idx] <- rnorm(length(idx), mean=m, sd=sd)

  return(sim.dat)
}

##################################################################################
##################################################################################

# simulation exercise1

niter <- 70
conv1.iter <- NULL
par1.iter <- NULL
for (i in 1:niter){
  print(i)
  dat <- sim.f(ini2, bio.dat)
  aux <- try(optim(inicf1, loglikcf.f, bio.dat=dat, method="L-BFGS-B", control=list(fnscale=-1, maxit=50000)), silent=T)
  if(class(aux)=="try-error"){
    conv1.iter <- c(conv1.iter, NA)
    par1.iter <- rbind(par1.iter, rep(NA, length(inicf1)))
  }else{
    conv1.iter <- c(conv1.iter, aux$convergence)
    par1.iter <- rbind(par1.iter, aux$par)
  }  
}

# simulation exercise

niter <- 70
conv2.iter <- NULL
par2.iter <- NULL
for (i in 1:niter){
  print(paste(i, Sys.time(), sep=" AT "))
  dat <- sim.f(ini2, bio.dat)
  aux <- try(optim(inicf1, loglikcf.f, bio.dat=dat, method="CG", control=list(fnscale=-1, maxit=500)), silent=T)
  if(class(aux)=="try-error"){
    conv2.iter <- c(conv2.iter, NA)
    par2.iter <- rbind(par2.iter, rep(NA, length(inicf1)))
  }else{
    conv2.iter <- c(conv2.iter, aux$convergence)
    par2.iter <- rbind(par2.iter, aux$par)
  }  
}

# check results visually

sum(conv.iter==0, na.rm=T)/niter

boxplot(par.iter)
lines(inicf2, col=2)

boxplot(par.iter[conv.iter==0,])
lines(inicf2, col=2)
abline(v=c(2.5,5.5,6.5,27+6.5,27+8.5,54+8.5,81+8.5), lty=2, col=3)


inicf2 <- c(logpsidepm.ini2,
            logpsiac.ini2,
            xidepm.ini2,
            xiac.ini2,
            xicatch.ini2,
            logB0.ini2,
            logR.ini2,
            logsage1sem1.ini2,
            logsage1sem2.ini2,
            logfsem1.ini2,
            logfsem2.ini2,
            logG1.ini2,
            logG2.ini2)
            
            
