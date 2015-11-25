# cohort cross validation
gcv <- function(stk, ids, cohorts=range(stk)[c("minyear")]:range(stk)[c("maxyear")]){
	coh <- as.character(cohorts)
	cv <- vector(length=length(coh))
	names(cv) <- coh
	cth <- catch.n(stk)
	for(i in coh){
		cat(i, ", ", sep="")
		flc <- FLCohort(catch.n(stk))
		flc[,i] <- NA
		catch.n(stk) <- as(flc, "FLQuant")
		fit <- sca(stk, ids)
		cv[i] <- mean(log(FLCohort(cth)[,i]/FLCohort(catch.n(fit))[,i])^2, na.rm=T)
		catch.n(stk) <- cth
	}
	cat("\n")
	cv
}

