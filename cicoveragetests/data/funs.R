###############################################################################
# EJ(20120413)
# Auxiliary functions for LH generation of datasets
# Parameters: 
#	a1, sr, sl = 50% selectivity age, right variance, left variance
#	s, v = steepness and virgin biomass for S/R
#	M1, M2 = two components of M
#	a50, asym = age of 50% maturity, ???
#	linf, k, t0 = vonBertalanffy growth pars
#	a, b = length~weight pars
#	bg ?? ato95 ??
###############################################################################

iterMedians <- function(x, ...){
	return(apply(x, c(1:5), median, na.rm = FALSE))
}

iterSums <- function(x, ...){
	return(apply(x, c(1:5), sum, na.rm = FALSE))
}

iterCv <- function(object, ...){
	sqrt(iterVars(object))/iterMeans(object)
}

#==============================================================================
# gislasim - cleaned version
#==============================================================================

setGeneric("gislasim", function(linf, ...) standardGeneric("gislasim"))

setMethod("gislasim", signature(linf="numeric"), function (linf, t0 = -0.1, a = 1e-05, b = 3, ato95 = 1, sl = 2, sr = 5000, s = 0.9, v = 1000, asym=1, bg=b, iter=1, k="missing", M1="missing", M2="missing", a50="missing", a1="missing"){
    if(missing(k))  k <- 3.15 * linf^(-0.64)
    if(missing(M1)) M1 <- 0.55 + 1.44 * log(linf) + log(k) 
    if(missing(M2)) M2 <- -1.61
    if(missing(a50)) a50 <- FLBRP:::invVonB(FLPar(linf=linf, t0=t0, k=k), 0.72 * linf^0.93)
    if(missing(a1)) a1 <- a50
    par <- FLPar(linf=linf, k=k, t0 = t0, a = a, b = b, asym=asym, bg=bg, sl=sl, sr=sr, s=s, v=v, M1=M1, M2=M2, ato95 = ato95, a50=a50, a1=a1, iter=iter)
    attributes(par)$units = c("cm", "kg", "1000s")
    return(par)
})

setMethod("gislasim", signature(linf="FLPar"), function (linf){
    # Renaming to avoid confusing the argument with the object.
    # linf here is an FLPar object that can contain several parameters 
    object <- linf
    rm(linf)
    # now the real thing
    v0 <- dimnames(object)$params	    
    if(!("linf" %in% v0)) stop("The function requires linf.")
    par <- FLPar(c(linf=NA, t0 = -0.1, a = 1e-05, b = 3, ato95 = 1, sl = 2, sr = 5000, s = 0.9, v = 1000, asym=1, bg=3, k=NA, M1=NA, M2=NA, a50=NA, a1=NA), iter=ncol(object))
    dimnames(par)$iter <- dimnames(object)$iter 
    par[dimnames(object)$params] <- object
    if(!("bg" %in% v0)) par["bg"] = par["b"]
    if(!("k" %in% v0)) par["k"] = 3.15 * par["linf"]^(-0.64)
    if(!("M1" %in% v0)) par["M1"] = 0.55 + 1.44 * log(par["linf"]) + log(par["k"])
    if(!("M2" %in% v0)) par["M2"] = -1.61
    if(!("a50" %in% v0)) par["a50"] = FLBRP:::invVonB(FLPar(linf=par["linf"], t0=par["t0"], k=par["k"]), c(0.72 * par["linf"]^0.93))
    if(!("a1" %in% v0)) par["a1"] = par["a50"]
    attributes(par)$units = c("cm", "kg", "1000s")
    return(par)
})

#==============================================================================
# genFunctions
#==============================================================================

genObs <- function(stock, qcv=0.1, ccv=0.1) {
	ages <- 1:range(stock)["maxfbar"]
	n <- stock.n(stock)[ac(ages)]
	z <- harvest(stock)[ac(ages)] + m(stock)[ac(ages)]
	logq <- -exp(-exp(0.2 * ages)) - 3 # trawl like catchability
	# observe index in 1st quarter with qcv
	index <- FLIndex(index = n[ac(ages)] * exp(-0.25 * z[ac(ages)]) * exp(logq + rnorm(prod(dim(n[ac(ages)])), 0, qcv)))  # 10% cv
	# OR observe index in 1st quarter with qcv and cov N
	#Sig2 <- cor(t(n[drop = T]))*qcv^2
	#index <- FLIndex(index = n * exp(-0.25 * z) * c(exp(logq + t(mvrnorm(dim(n)[2], rep(0,dim(n)[1]), Sig2)))))
	range(index)[c("startf","endf")] <- 0.25
	# observe catch with ccv
	catch.n(stock) <- catch.n(stock) * exp(rnorm(prod(dim(catch.n(stock))), 0, ccv))
	catch(stock) <- computeCatch(stock)
	list(stock = stock, index = list(index))
}

genIdx <- function(stock, qcv=0.1, ccv=0.1) {
	ages <- 1:range(stock)["maxfbar"]
	n <- stock.n(stock)[ac(ages)]
	z <- harvest(stock)[ac(ages)] + m(stock)[ac(ages)]
	logq <- -exp(-exp(0.2 * ages)) - 3 # trawl like catchability
	# observe index in 1st quarter with qcv
	index <- FLIndex(index = n[ac(ages)] * exp(-0.25 * z[ac(ages)]) * exp(logq + rnorm(prod(dim(n[ac(ages)])), 0, qcv)))  # 10% cv
	# OR observe index in 1st quarter with qcv and cov N
	#Sig2 <- cor(t(n[drop = T]))*qcv^2
	#index <- FLIndex(index = n * exp(-0.25 * z) * c(exp(logq + t(mvrnorm(dim(n)[2], rep(0,dim(n)[1]), Sig2)))))
	range(index)[c("startf","endf")] <- 0.25
	list(stock = stock, index = list(index))
}

#==============================================================================
# catch constraint
#==============================================================================

ccstr <- function(x, lo=0.85, hi=1.15){
	if(x>hi) x <- hi
	if(x<lo) x <- lo
	x
}

#==============================================================================
# hcr projection
#==============================================================================

cat3hcr <- function(data, HCR, iyr = 1, fyr = 25, srviyr= -9, oem.cv = 0.2, iem.cv = 0.2, sr.cv=0.2, cc=FALSE, rndseed=12345, ul=0.975, ll=0.33, upr=1.05, lor=0.75){

	# hcr: cat3 , ic, cat3len, iclen
	# lh
	stk <- data$stock
	lhPars <- attr(stk, "lhPars")
	refpts <- attr(stk, "refpts")
	k <- lhPars["k"] 
	linf <- lhPars["linf"] 
	t0 <- lhPars["t0"] 
	srModel <- "bevholt" # sr hard coded ...
	v <- lhPars["v"]
	s <- lhPars["s"]
	# age 50% retention
	a1 <- lhPars["a1"]
	sl <- lhPars["sl"]
	
	# stock object
	stk <- qapply(stk, function(x){dimnames(x)[[2]] <- -15:0; x})
	range(stk)[4:5] <- c(-15,0)
	stk <- window(stk, stf=list(nyears=fyr-iyr+1))
	
	# index
	idx <- qapply(data$index[[1]], function(x){dimnames(x)[[2]] <- -15:0; x})
	idx <- window(idx, start=srviyr, end=fyr)
	idx.wt <- stock.wt(stk)[dimnames(index(idx))[[1]],1]
	
	# S/R
 	sr <- FLSR(model=srModel)
 	params(sr) <- FLPar(abPars(srModel, c(refpts["virgin","ssb"]/refpts["virgin","rec"]), v=c(v), s=c(s)))
	set.seed(rndseed)
	sr.res <- FLQuant(rlnorm(fyr,0,log(sr.cv^2+1)), dimnames=list(year=iyr:fyr))

	# objects for results
	Cadv <- catch(stk)
	Cadv[,ac(iyr)] <- mean(Cadv[,ac(iyr-c(1:3))])
	idxRatio <- catch(stk)
	idxRatio[] <- NA
	Lsq <- idxRatio

	# loop
	set.seed(rndseed)
	for(i in iyr:(fyr-1)){
		# biomass index for advice
		idxadv <- window(index(idx), end=i-1)
		dnms2 <- dimnames(idxadv)[[2]]
		nobs <- length(dnms2)
		idxadv <- quantSums(idxadv*idx.wt[,rep(1,nobs)])
		if(HCR=="cat3"){
			# mean index last 2 years
			i1 <- mean(idxadv[,ac(i-c(1,2))])
			# mean index previous 3 years
			i2 <- mean(idxadv[,ac(i-c(3:5))])
			if(cc) idxRatio[,ac(i-1)] <- ccstr(i1/i2) else idxRatio[,ac(i-1)] <- i1/i2
		}
		if(HCR=="cat3len"){
			# mean index last 2 years
			i1 <- mean(idxadv[,ac(i-c(1,2))])
			# mean index previous 3 years
			i2 <- mean(idxadv[,ac(i-c(3:5))])
			# Lc
			ages <- as.numeric(dimnames(catch.n(stk))[[1]])
			lenbar <- vonB(FLPar(k=k, linf=linf, t0=t0), ages+0.5)
			#Lc <- lenbar[1]
			Lc <- floor(qnorm(0.25, a1, sl)) + 0.5
			Lc <- vonB(FLPar(k=k, linf=linf, t0=t0), Lc)
			Lmf <- (3*Lc+linf)/4
			#Lsq
			Lsq <- weighted.mean(lenbar, c(catch.n(stk)[,ac(i-1)]))
			if(cc) idxRatio[,ac(i-1)] <- ccstr(i1/i2*Lsq/Lmf) else idxRatio[,ac(i-1)] <- i1/i2*Lsq/Lmf
		}
		# idx <> IC
		if(HCR=="ic"){
			idx1 <- idxadv[,ac(i-1)]
			idx0 <- window(idxadv, end=i-2) 
			idxRatio[,ac(i-1)] <- 1
	 		if(idx1 > yearMeans(idx0) + qnorm(ul)*sqrt(yearVars(idx0)/nobs)) idxRatio[,ac(i-1)] <- upr
	 		if(idx1 < yearMeans(idx0) + qnorm(ll)*sqrt(yearVars(idx0)/nobs)) idxRatio[,ac(i-1)] <- lor
			if(cc) idxRatio[,ac(i-1)] <- ccstr(idxRatio[,ac(i-1)])
		}	
		# idx <> IC * len
		if(HCR=="iclen"){
			idx1 <- idxadv[,ac(i-1)]
			idx0 <- window(idxadv, end=i-2) 
			idxRatio[,ac(i-1)] <- 1
			if(idx1 > yearMeans(idx0) + qnorm(0.975)*sqrt(yearVars(idx0)/nobs)) idxRatio[,ac(i-1)] <- 1.05
			if(idx1 < yearMeans(idx0) + qnorm(0.33)*sqrt(yearVars(idx0)/nobs)) idxRatio[,ac(i-1)] <- 0.75
			# Lc
			ages <- as.numeric(dimnames(catch.n(stk))[[1]])
			lenbar <- vonB(FLPar(k=k, linf=linf, t0=t0), ages+0.5)
			#Lc <- lenbar[1]
			Lc <- floor(qnorm(0.25, a1, sl)) + 0.5
			Lc <- vonB(FLPar(k=k, linf=linf, t0=t0), Lc)
			Lmf <- (3*Lc+linf)/4
			#Lsq
			Lsq <- weighted.mean(lenbar, c(catch.n(stk)[,ac(i-1)]))
			if(cc) idxRatio[,ac(i-1)] <- ccstr(idxRatio[,ac(i-1)]*Lsq/Lmf) else idxRatio[,ac(i-1)] <- idxRatio[,ac(i-1)]*Lsq/Lmf
		}	
		# len
		if(HCR=="len"){
			# Lc
			ages <- as.numeric(dimnames(catch.n(stk))[[1]])
			lenbar <- vonB(FLPar(k=k, linf=linf, t0=t0), ages+0.5)
			#Lc <- lenbar[1]
			Lc <- floor(qnorm(0.25, a1, sl)) + 0.5
			Lc <- vonB(FLPar(k=k, linf=linf, t0=t0), Lc)
			Lmf <- (3*Lc+linf)/4
			#Lsq
			Lsq[,ac(i-1)] <- weighted.mean(lenbar, c(catch.n(stk)[,ac(i-1)]))
			if(cc) idxRatio[,ac(i-1)] <- ccstr(Lsq[,ac(i-1)]/Lmf) else idxRatio[,ac(i-1)] <- Lsq[,ac(i-1)]/Lmf	
		}	
		# cacth for next year
		Cadv[,ac(i+1)] <- Cadv[,ac(i-1)]*idxRatio[,ac(i-1)]*rlnorm(1,0,sqrt(log(iem.cv^2+1)))
		# note: will need rewrite to account for iters
		trg <- fwdControl(data.frame(year=i, quantity='catch', val=c(Cadv[,ac(i)])))
		stk[,ac(i)] <- fwd(stk, ctrl=trg, sr=sr, sr.residuals = sr.res[,ac(i)], sr.residuals.mult = TRUE)[,ac(i)]

		#stk <- fwd(stk, ctrl=trg, sr=sr)
		idx[,ac(i)] <- genIdx(stk, qcv = oem.cv)$index[[1]][, ac(i)]
	}
	# output
	flqs <- FLQuants(ssb=ssb(stk), catch=catch(stk), rec=rec(stk), f=fbar(stk), catch.adv=Cadv, idxRatio=idxRatio, Lsq=Lsq)
	attr(flqs, "lhPars") <- lhPars
	attr(flqs, "refpts") <- refpts
	flqs	
}

#==============================================================================
# extract	
#==============================================================================

getStkInfo <- function(object, scn, hcr){
	res <- lapply(object, function(x){
		lh <- as.data.frame(attr(x, "lhPars"))
		rp <- as.data.frame(attr(x, "refpts")[c("msy", "f0.1", "fmax"), "harvest"])
		data.frame(params=c(ac(lh[,1]), ac(rp[,1])), value=c(lh[,3], rp[,4]), scn=scn, hcr=hcr)	
	})
	
	res <- do.call("rbind", res)
	res$stk <- unlist(lapply(strsplit(rownames(res), "[.]"), "[[", 1))
	rownames(res) <- NULL
	res
	
}
	
getRes <- function(object, scn, hcr){
	res <- lapply(object, as.data.frame)
	res <- do.call("rbind", res)
	res$stk <- unlist(lapply(strsplit(rownames(res), "[.]"), "[[", 1))
	rownames(res) <- NULL
	res$scn <- scn
	res$hcr <- hcr
	res
}

