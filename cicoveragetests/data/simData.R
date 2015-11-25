library(FLBRP)
library(plyr)
library(Hmisc)
library(parallel)
library(reshape)
source("funs.R")
options(mc.cores=20)


# read data
wklifeLst <- read.table("allStockslifeHistoryParam.txt", sep="\t", head=TRUE)
wklifeLst <-transform(wklifeLst, value=as.numeric(as.character(value)))
wklifeLst <- wklifeLst[!is.na(wklifeLst$value),]
# remove "fle-2232", "spr-nsea" and pol-celt the parameters are duplicated and seem inconsistent
wklifeLst <- subset(wklifeLst, stock != "fle-2232")
wklifeLst <- subset(wklifeLst, stock != "spr-nsea")
wklifeLst <- subset(wklifeLst, stock != "pol-celt")

# use a and b from her-nits for her-31
wklifeLst <- rbind(wklifeLst, data.frame(stock="her-31", paramlong=c("a","b"), value=c(0.0036, 3.03871), param=c("a","b")))

set.seed(123)
wklife.brp <- mclapply(split(wklifeLst, wklifeLst$stock), function(x, s=0.75){
    cat(as.character(x$stock)[1], "\n")
    # get parameters
    par <- FLPar(x$value, tolower(x$param))
    if(!("linf" %in% dimnames(par)$params) & "lmax" %in% dimnames(par)$params){
        dimnames(par)$params[dimnames(par)$params == "lmax"] <- "linf"
    } 
    if("linf" %in% dimnames(par)$params){
      # complete with gislasim
      dnms <- dimnames(par)$params
      par <- par[dnms %in% dimnames(gislasim(0))$params]
      par <- gislasim(par)
	  par["s"] <- s
	  # set max age if not available compute from linf
	  if(!("tmax" %in% x$param)){ 
	  	tmax <- floor(invVonB(FLPar(c(linf=par["linf"], k=par["k"], t0=par["t0"])), par["linf"]-1))
		} else {
			tmax <- x[x$param=="tmax", "value"]
		}
		# run lh
		res <- lh(par, range=c(min=1, max=tmax, minfbar=ceiling(tmax/10), maxfbar=floor(tmax/2), plusgroup=tmax))    
		res@desc <- as.character(x$stock[1])
		attr(res, "lhPars") <- par
		attr(res, "refpts") <- refpts(res)
		res
    } else {
      NULL
    }
})
wklife.brp <- wklife.brp[!unlist(mclapply(wklife.brp, is.null))]

wklife.stk <- mclapply(wklife.brp, function(x){
	cat(x@desc, "\n")
	lhPars <- attr(x, "lhPars")
	refpts <- attr(x, "refpts")
	Fc <- refpts["msy","harvest"]*2
	if(!is.na(Fc)){	
		Fmsy <- refpts["msy","harvest"]
		Ftrg <- c(seq(0, Fc, len=19), rep(Fc, 20), seq(Fc, Fmsy, len=10))
		trg <- fwdControl(data.frame(year=c(2:50), quantity=rep('f', 49), val=Ftrg))
		ex.stk <- as(x, "FLStock")[,1:50]
		# S/R
		srModel <- "bevholt" # sr hard coded ...
		v <- lhPars["v"]
		s <- lhPars["s"]
	 	ex.sr <- FLSR(model=srModel)
	 	params(ex.sr) <- FLPar(abPars(srModel, c(refpts["virgin","ssb"]/refpts["virgin","rec"]), v=c(v), s=c(s)))
		#ex.sr <- fmle(as.FLSR(ex.stk, model="bevholt"), control=list(trace=0))
		stk <- fwd(ex.stk, ctrl=trg, sr=ex.sr)
		stk@name <- x@desc
		stk@desc <- paste("simulated data loosely based on", x@desc)
		attr(stk, "lhPars") <- lhPars
		attr(stk, "refpts") <- refpts
		stk  
	} else {
    	NULL
  	}
})
wklife.stk <- wklife.stk[!unlist(mclapply(wklife.stk, is.null))]
wklife.stk <- FLStocks(wklife.stk)

# subset full time series to get just developing period
stks01 <- window(wklife.stk, start=5, end=20)
# generating survey index and adding observation error to catches
set.seed(239246)
stks01 <- mclapply(stks01, genObs, ccv=0.2, qcv=0.2)

# subset full time series to get just stable at high exploitation period
stks03 <- window(wklife.stk, 25,40)
# generating survey index and adding observation error to catches
set.seed(239246)
stks03 <- mclapply(stks03, genObs, ccv=0.2, qcv=0.2)

save(wklife.brp,wklife.stk,stks03,stks01, file="RData.data")



