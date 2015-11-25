################################################################################

# function to compute biomass given the model estimates

################################################################################

calc.pop <- function(mat, Y, dsurv=0.375, d=c(0.5,0.5)){

  dsem1 <- d[1]
  dsem2 <- d[2]
  
  fsem1 <- mat[match(paste("fsem1[",1:Y,"]",sep=""), names(mat))]
  fsem2 <- mat[match(paste("fsem2[",1:Y,"]",sep=""), names(mat))]
  B0 <- exp(mat[match("logB0", names(mat))])
  R <- exp(mat[match(paste("logR[",1:Y,"]",sep=""), names(mat))])
  sage1sem1 <- mat[match("sage1[1]", names(mat))] 
  sage1sem2 <- mat[match("sage1[2]", names(mat))]
  M1 <- exp(mat[match("logM1", names(mat))]) 
  M2 <- exp(mat[match("logM2", names(mat))]) 
  G1 <- exp(mat[match("logG1", names(mat))]) 
  G2 <- exp(mat[match("logG2", names(mat))]) 
  
# initialise B1, Btot, C1sem1, C2plussem1, Ctotsem1, C1sem2, C2plussem2, Ctotsem2

  B1 <- rep(0, Y)        # at survey time
  B2plus <- rep(0, Y)    # at survey time
  Btot <- rep(0, Y)      # at survey time
  
#  B1popsem1 <- rep(0,Y)
#  B2pluspopsem1 <- rep(0,Y)
#  B1popsem2 <- rep(0,Y)
#  B2pluspopsem2 <- rep(0,Y)
  
  C1sem1 <- rep(0, Y)
  C2plussem1 <- rep(0, Y)
  Ctotsem1 <- rep(0, Y)
  C1sem2 <- rep(0, Y)
  C2plussem2 <- rep(0, Y)
  Ctotsem2 <- rep(0, Y)   

# calculations for the first year

  B1[1] <- R[1] * exp( (G1-M1-fsem1[1]*sage1sem1) * dsurv)
  B2plus[1] <- B0 * exp( (G2-M2-fsem1[1]) * dsurv) # sage2 is fixed at 1
  Btot[1] <- B1[1] + B2plus[1]
  C1sem1[1] <- R[1] * (1-exp((G1-M1-fsem1[1]*sage1sem1)*dsem1)) * (fsem1[1]*sage1sem1) / (M1+fsem1[1]*sage1sem1-G1) 
  C2plussem1[1] <- B0 * (1-exp((G2-M2-fsem1[1])*dsem1)) * (fsem1[1]) / (M2+fsem1[1]-G2)
  C1sem2[1] <- (R[1] * exp( (G1-M1-fsem1[1]*sage1sem1) * dsem1)) * (1-exp((G1-M1-fsem2[1]*sage1sem2)*dsem2)) * (fsem2[1]*sage1sem2) / (M1+fsem2[1]*sage1sem2-G1)
  C2plussem2[1] <- (B0 * exp( (G2-M2-fsem1[1]) * dsem1)) * (1-exp((G2-M2-fsem2[1])*dsem2)) * (fsem2[1]) / (M2+fsem2[1]-G2)
  Ctotsem1[1] <- C1sem1[1] + C2plussem1[1]
  Ctotsem2[1] <- C1sem2[1] + C2plussem2[1]

# compute catches in the first and second semester of the first year

# compute B1 and Btot for the rest years

  for (i in 2:Y){
    B1[i] <- R[i] * exp( (G1-M1-fsem1[i]*sage1sem1) * dsurv)
    aux <- B1[i-1] * exp((G1-M1)*(1-dsurv) - fsem1[i-1]*sage1sem1*(dsem1-dsurv) - fsem2[i-1]*sage1sem2*dsem2) + 
           B2plus[i-1] * exp((G2-M2)*(1-dsurv) - fsem1[i-1]*(dsem1-dsurv) - fsem2[i-1]*dsem2) # sage2=1 ; aux is biomass at age 2+ at the beginning of year i 
    B2plus[i] <- aux * exp( (G2-M2-fsem1[i]) * dsurv)
    Btot[i] <- B1[i] + B2plus[i]
    C1sem1[i] <- R[i] * (1-exp((G1-M1-fsem1[i]*sage1sem1)*dsem1)) * (fsem1[i]*sage1sem1) / (M1+fsem1[i]*sage1sem1-G1) 
    C2plussem1[i] <- aux * (1-exp((G2-M2-fsem1[i])*dsem1)) * (fsem1[i]) / (M2+fsem1[i]-G2)
    C1sem2[i] <- (R[i] * exp( (G1-M1-fsem1[i]*sage1sem1) * dsem1)) * (1-exp((G1-M1-fsem2[i]*sage1sem2)*dsem2)) * (fsem2[i]*sage1sem2) / (M1+fsem2[i]*sage1sem2-G1)
    C2plussem2[i] <- (aux * exp( (G2-M2-fsem1[i]) * dsem1)) * (1-exp((G2-M2-fsem2[i])*dsem2)) * (fsem2[i]) / (M2+fsem2[i]-G2)
    Ctotsem1[i] <- C1sem1[i] + C2plussem1[i]
    Ctotsem2[i] <- C1sem2[i] + C2plussem2[i]

  }
  
  out <- c(B1, Btot, C1sem1, Ctotsem1, C1sem2, Ctotsem2)
  names(out) <- c( paste("B1[",1:Y,"]", sep=""), paste("Btot[",1:Y,"]", sep=""), 
                   paste("C1sem1[",1:Y,"]", sep=""), paste("Ctotsem1[",1:Y,"]", sep=""),
                   paste("C1sem2[",1:Y,"]", sep=""), paste("Ctotsem2[",1:Y,"]", sep=""))
  return(out)
}
