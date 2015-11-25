################################################################################
################################################################################

# functions to compute the log likelihood to be maximised for the CBBM model
# version with and without closed-form solutions for some of the parameters

# leire ibaibarriaga, preparation of stecf ewg 2014-03

################################################################################
################################################################################

# function to compute the log likelihood, without closed-forms solutions

loglik.f <- function(param, bio.dat){
  Y <- dim(bio.dat)[1] # TO BE CHECKED
  if (length(param)!= (3*Y+16)){
    warning("Dimensions do not match")  
  }
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
  
  psidepm <- exp(logpsidepm)
  psiac <- exp(logpsiac)
  psirobs <- exp(logpsirobs)
  psig <- exp(logpsig)
  sage1sem1 <- exp(logsage1sem1)
  sage1sem2 <- exp(logsage1sem2)
  B0 <- exp(logB0)
  R <- exp(logR)
  G1 <- exp(logG1)
  G2 <- exp(logG2)
  fsem1 <- exp(logfsem1)
  fsem2 <- exp(logfsem2)
  
  fsem1[bio.dat$Ctotsem1==0] <- 0 # when there are no catches we set f to zero
  fsem2[bio.dat$Ctotsem2==0] <- 0
  
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
  
  # compute log likelihood
  
  ll <- 0
  #  print("BP.depm")
  for (i in 1:length(year.BP.depm)){
    idx <- year.BP.depm[i]
    sh1 <- exp(xidepm)*B1[idx]/Btot[idx]
    sh2 <- exp(xidepm)*(1-B1[idx]/Btot[idx])  
    aux <- dbeta(x=bio.dat$BP.depm[idx], shape1=sh1, shape2=sh2, log=T)
    ll <- ll + aux 
  }
  #  print("Btot.depm")
  for (i in 1:length(year.Btot.depm)){
    idx <- year.Btot.depm[i]
    ml <- logqdepm + log(Btot[idx])
    sdl <- sqrt((psidepm+bio.dat$psidepm[idx])/(psidepm*bio.dat$psidepm[idx]))
    aux <- dlnorm(x=bio.dat$Btot.depm[idx], meanlog=ml, sdlog=sdl, log=T)
    ll <- ll + aux 
  }  
  #  print("BP.ac")
  for (i in 1:length(year.BP.ac)){
    idx <- year.BP.ac[i]
    sh1 <- exp(xiac)*B1[idx]/Btot[idx]
    sh2 <- exp(xiac)*(1-B1[idx]/Btot[idx])  
    aux <- dbeta(x=bio.dat$BP.ac[idx], shape1=sh1, shape2=sh2, log=T)
    ll <- ll + aux 
  }
  #  print("Btot.ac")
  for (i in 1:length(year.Btot.ac)){
    idx <- year.Btot.ac[i]
    ml <- logqac + log(Btot[idx])
    sdl <- sqrt((psiac+bio.dat$psiac[idx])/(psiac*bio.dat$psiac[idx]))
    aux <- dlnorm(x=bio.dat$Btot.ac[idx], meanlog=ml, sdlog=sdl, log=T)
    ll <- ll + aux 
  }
  #  print("CPsem1")  
  for (i in 1:length(year.CPsem1)){
    idx <- year.CPsem1[i]
    sh1 <- exp(xicatch)*C1sem1[idx]/Ctotsem1[idx]
    sh2 <- exp(xicatch)*(1-C1sem1[idx]/Ctotsem1[idx])  
    aux <- dbeta(x=bio.dat$CPsem1[idx], shape1=sh1, shape2=sh2, log=T)
    ll <- ll + aux 
  }
  #  print("Ctotsem1")  
  for (i in 1:length(year.Ctotsem1)){
    idx <- year.Ctotsem1[i]
    ml <- log(Ctotsem1[idx])
    sdl <- 1/sqrt(400)
    aux <- dlnorm(x=bio.dat$Ctotsem1[idx], meanlog=ml, sdlog=sdl, log=T)
    ll <- ll + aux 
  }
  #  print("CPsem2")  
  for (i in 1:length(year.CPsem2)){
    idx <- year.CPsem2[i]
    sh1 <- exp(xicatch)*C1sem2[idx]/Ctotsem2[idx]
    sh2 <- exp(xicatch)*(1-C1sem2[idx]/Ctotsem2[idx])  
    aux <- dbeta(x=bio.dat$CPsem2[idx], shape1=sh1, shape2=sh2, log=T)
    ll <- ll + aux 
  }
  #  print("Ctotsem2")  
  for (i in 1:length(year.Ctotsem2)){
    idx <- year.Ctotsem2[i]
    ml <- log(Ctotsem2[idx])
    sdl <- 1/sqrt(400)
    aux <- dlnorm(x=bio.dat$Ctotsem2[idx], meanlog=ml, sdlog=sdl, log=T)
    ll <- ll + aux 
  }
  #  print("Robs")  
  for (i in 1:length(year.Robs)){
    idx <- year.Robs[i]
    ml <- logqrobs + exp(logkrobs)*log(R[idx])
    sdl <- 1/sqrt(psirobs)
    aux <- dlnorm(x=bio.dat$Robs[idx], meanlog=ml, sdlog=sdl, log=T)
    ll <- ll + aux 
  }
  #  print("G1") 
  for (i in 1:length(year.G1)){
    idx <- year.G1[i]
    m <- G1
    sd <- 1/sqrt(psig)
    aux <- dnorm(x=bio.dat$G1[idx], mean=m, sd=sd, log=T)
    ll <- ll + aux 
  }
  #  print("G2") 
  for (i in 1:length(year.G2)){
    idx <- year.G2[i]
    m <- G2
    sd <- 1/sqrt(psig)
    aux <- dnorm(x=bio.dat$G2[idx], mean=m, sd=sd, log=T)
    ll <- ll + aux 
  }
  return (ll)
}

################################################################################
################################################################################

# function to derive closed-form of some parameters given all the rest 

cf <- function(Btot, R, G1, G2, psidepm, psiac, bio.dat){
  Y <- dim(bio.dat)[1] 
  
  # create index variables for observation equations (avoiding missing data)
  
  year.Btot.depm <- (1:Y)[!is.na(bio.dat$Btot.depm)]
  year.Btot.ac <- (1:Y)[!is.na(bio.dat$Btot.ac)]
  year.Robs <- (1:Y)[!is.na(bio.dat$Robs)]
  year.G1 <- (1:Y)[!is.na(bio.dat$G1)] 
  year.G2 <- (1:Y)[!is.na(bio.dat$G2)]
  
  wtdepm <- psidepm*bio.dat$psidepm[year.Btot.depm]/(psidepm+bio.dat$psidepm[year.Btot.depm])
  logqdepm <- sum(wtdepm*(log(bio.dat$Btot.depm[year.Btot.depm])-log(Btot[year.Btot.depm])))/sum(wtdepm)
  wtac <- psiac*bio.dat$psiac[year.Btot.ac]/(psiac+bio.dat$psiac[year.Btot.ac])
  logqac <- sum(wtac*(log(bio.dat$Btot.ac[year.Btot.ac])-log(Btot[year.Btot.ac])))/sum(wtac)
  aa <- - sum(log(bio.dat$Robs[year.Robs]))
  bb <- length(year.Robs)
  cc <- sum(log(R[year.Robs]))
  dd <- sum(log(bio.dat$Robs[year.Robs])*log(R[year.Robs]))
  ee <- - cc
  ff <- - sum(log(R[year.Robs])^2)                                    
  krobs <- (aa*ee-dd*bb)/(ff*bb-cc*ee)  
  logqrobs <- -cc / bb * krobs - aa/bb
  psirobs <- length(year.Robs)/sum( (log(bio.dat$Robs[year.Robs]) - logqrobs - krobs*log(R[year.Robs]))^2 )                 
  psig <- (length(year.G1)+length(year.G2))/(sum((bio.dat$G1[year.G1]-G1)^2)+sum((bio.dat$G2[year.G2]-G2)^2))                  
  out <- c(logqdepm, logqac, krobs, logqrobs, psirobs, psig)
  names(out) <- c("logqdepm", "logqac", "krobs", "logqrobs", "psirobs", "psig")
  return (out)
}

################################################################################
################################################################################

# function to compute the log likelihood, with closed-forms solutions

loglikcf.f <- function(param, bio.dat){
  Y <- dim(bio.dat)[1] # TO BE CHECKED
  if (length(param)!= (3*Y+10)){
    warning("Dimensions do not match")  
  }
  logpsidepm <- param[1]
  logpsiac <- param[2]
  xidepm <- param[3]
  xiac <- param[4]
  xicatch <- param[5]
  logB0 <- param[6]
  logR <- param[(6+1):(6+Y)]
  logsage1sem1 <- param[6+Y+1]
  logsage1sem2 <- param[6+Y+2]
  logfsem1 <- param[(8+Y+1):(8+Y+Y)]
  logfsem2 <- param[(8+2*Y+1):(8+2*Y+Y)]
  logG1 <- param[8+3*Y+1]
  logG2 <- param[8+3*Y+2]
  
  dsurv <- 0.375
  dsem1 <- 0.5
  dsem2 <- 0.5
  M1 <- 0.8
  M2 <- 1.2
  
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
  
  fsem1[bio.dat$Ctotsem1==0] <- 0 # when there are no catches we set f to zero
  fsem2[bio.dat$Ctotsem2==0] <- 0
  
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
  
  rest <- cf(Btot, R, G1, G2, psidepm, psiac, bio.dat)
  logqdepm <- rest[match("logqdepm", names(rest))]
  logqac <- rest[match("logqac", names(rest))]
  logkrobs <- log(rest[match("krobs", names(rest))])
  logqrobs <- rest[match("logqrobs", names(rest))]
  psirobs <- rest[match("psirobs", names(rest))]
  psig <- rest[match("psig", names(rest))]
  
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
  
  # compute log likelihood
  
  ll <- 0
  #  print("BP.depm")
  for (i in 1:length(year.BP.depm)){
    idx <- year.BP.depm[i]
    sh1 <- exp(xidepm)*B1[idx]/Btot[idx]
    sh2 <- exp(xidepm)*(1-B1[idx]/Btot[idx])  
    aux <- dbeta(x=bio.dat$BP.depm[idx], shape1=sh1, shape2=sh2, log=T)
    ll <- ll + aux 
  }
  #  print("Btot.depm")
  for (i in 1:length(year.Btot.depm)){
    idx <- year.Btot.depm[i]
    ml <- logqdepm + log(Btot[idx])
    sdl <- sqrt((psidepm+bio.dat$psidepm[idx])/(psidepm*bio.dat$psidepm[idx]))
    aux <- dlnorm(x=bio.dat$Btot.depm[idx], meanlog=ml, sdlog=sdl, log=T)
    ll <- ll + aux 
  }  
  #  print("BP.ac")
  for (i in 1:length(year.BP.ac)){
    idx <- year.BP.ac[i]
    sh1 <- exp(xiac)*B1[idx]/Btot[idx]
    sh2 <- exp(xiac)*(1-B1[idx]/Btot[idx])  
    aux <- dbeta(x=bio.dat$BP.ac[idx], shape1=sh1, shape2=sh2, log=T)
    ll <- ll + aux 
  }
  #  print("Btot.ac")
  for (i in 1:length(year.Btot.ac)){
    idx <- year.Btot.ac[i]
    ml <- logqac + log(Btot[idx])
    sdl <- sqrt((psiac+bio.dat$psiac[idx])/(psiac*bio.dat$psiac[idx]))
    aux <- dlnorm(x=bio.dat$Btot.ac[idx], meanlog=ml, sdlog=sdl, log=T)
    ll <- ll + aux 
  }
  #  print("CPsem1")  
  for (i in 1:length(year.CPsem1)){
    idx <- year.CPsem1[i]
    sh1 <- exp(xicatch)*C1sem1[idx]/Ctotsem1[idx]
    sh2 <- exp(xicatch)*(1-C1sem1[idx]/Ctotsem1[idx])  
    aux <- dbeta(x=bio.dat$CPsem1[idx], shape1=sh1, shape2=sh2, log=T)
    ll <- ll + aux 
  }
  #  print("Ctotsem1")  
  for (i in 1:length(year.Ctotsem1)){
    idx <- year.Ctotsem1[i]
    ml <- log(Ctotsem1[idx])
    sdl <- 1/sqrt(400)
    aux <- dlnorm(x=bio.dat$Ctotsem1[idx], meanlog=ml, sdlog=sdl, log=T)
    ll <- ll + aux 
  }
  #  print("CPsem2")  
  for (i in 1:length(year.CPsem2)){
    idx <- year.CPsem2[i]
    sh1 <- exp(xicatch)*C1sem2[idx]/Ctotsem2[idx]
    sh2 <- exp(xicatch)*(1-C1sem2[idx]/Ctotsem2[idx])  
    aux <- dbeta(x=bio.dat$CPsem2[idx], shape1=sh1, shape2=sh2, log=T)
    ll <- ll + aux 
  }
  #  print("Ctotsem2")  
  for (i in 1:length(year.Ctotsem2)){
    idx <- year.Ctotsem2[i]
    ml <- log(Ctotsem2[idx])
    sdl <- 1/sqrt(400)
    aux <- dlnorm(x=bio.dat$Ctotsem2[idx], meanlog=ml, sdlog=sdl, log=T)
    ll <- ll + aux 
  }
  #  print("Robs")  
  for (i in 1:length(year.Robs)){
    idx <- year.Robs[i]
    ml <- logqrobs + exp(logkrobs)*log(R[idx])
    sdl <- 1/sqrt(psirobs)
    aux <- dlnorm(x=bio.dat$Robs[idx], meanlog=ml, sdlog=sdl, log=T)
    ll <- ll + aux 
  }
  #  print("G1") 
  for (i in 1:length(year.G1)){
    idx <- year.G1[i]
    m <- G1
    sd <- 1/sqrt(psig)
    aux <- dnorm(x=bio.dat$G1[idx], mean=m, sd=sd, log=T)
    ll <- ll + aux 
  }
  #  print("G2") 
  for (i in 1:length(year.G2)){
    idx <- year.G2[i]
    m <- G2
    sd <- 1/sqrt(psig)
    aux <- dnorm(x=bio.dat$G2[idx], mean=m, sd=sd, log=T)
    ll <- ll + aux 
  }
  return (ll)
}

################################################################################
################################################################################
