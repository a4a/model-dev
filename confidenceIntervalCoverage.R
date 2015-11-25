summ <- function(object){
	fb <- fbar(object)
	sb <- ssb(object)
	rb <- rec(object)
	cb <- catch(object)
	dnms <- dimnames(fb)
	dnms[[1]] <- c("catch", "fbar", "ssb", "rec")
	names(dnms)[1] <- "summ"
	cb <- cb[rep(1,4)]
	cb[2] <- fb
	cb[3] <- sb
	cb[4] <- rb
	dimnames(cb) <- dnms
	cb
}

iterQuantile <- function(x, prob=0.5, ...){
	return(apply(x, c(1:5), quantile, prob=prob, na.rm = FALSE))
}

# fit
nit <- 250
stk <- cod1 
ids <- window(cod.tun[1], end=2013)
fmod <- ~ te(age, year, k=c(4,20)) + s(year, k=3, by=as.numeric(age==1))
qmod <- list(~s(age, k=4))
fit <- sca(stk, ids, fmodel=fmod, qmodel=qmod, fit="assessment")
fitsim <- simulate(fit, nit)
stk <- stk + fitsim
index(ids[[1]]) <- stock.n(stk)[ac(1:5),ac(1983:2013)]*(0.001+runif(100, 0.0005,0.0015))
range(ids[[1]])["endf"] <- 0

lst <- list()
length(lst) <- nit

for(i in 1:nit){
	cat(".")
	stk0 <- iter(stk, i)
	ids0 <- ids
	index(ids0[[1]]) <- index(iter(ids[[1]],i))
	lst[[i]] <- sca(stk0, ids0, fmodel=fmod, qmodel=list(~1), fit="assessment")
}

lst1 <- lapply(lst, function(x) try(simulate(x, 1000)))



	cat(".")
	try(
	stk <- cod1 + 
	flq <- summ(stk)
	cod.summ <- summ(cod1)
	cod.summ > iterQuantile(flq, 0.05) & cod.summ < iterQuantile(flq, 0.95)
	)	
})


