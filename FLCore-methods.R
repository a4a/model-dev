ssbn <- function(stk) quantSums(mat(stk)*stock.n(stk))

yearCumsum <- function(x, ...){
	x[] <- apply(x, c(1,3:6), cumsum)
	return(x)
}

yearDiffPerc <- function(x, ...){
	#x[,-1] <- x[,-1]/x[,-ncol(x)]-1
	x[,-1] <- x[,-1]/c(x[,1])-1
	x[,1] <- 0
	return(x)
}

as.table.FLQuants <- function(x){
	x <- mcf(x)
	df0 <- as.data.frame(do.call("cbind", lapply(x, c)))
	row.names(df0) <- dimnames(x)[[2]]
	df0
}

iterQuantiles <- function(x, prob=0.5, ...){
	return(apply(x, c(1:5), quantile, prob=prob, na.rm = FALSE))
}

an <- function(x, ...) as.numeric(x, ...)

cbind.FLQuant <- function(x, y){
	lst <- mcf(list(x, y))
	res <- lst[[1]]
	res[,dimnames(y)[[2]]] <- y
	res
} 


