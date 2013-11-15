
library(FLa4a)

data(ple4)

stk <- list(ple4, ple4, ple4)

FUNCs <- c("ssb", "catch", "fbar", "rec")

getSumm <- function(stk, FUNCs) {
  do.call(rbind,
    lapply(FUNCs, function(x) 
    {
      FUN <- match.fun(x)
      out <- as.data.frame(FUN(stk))
      out $ age <- paste(out $ age)
      cbind(out, fun = x)
    })
  )
}

dat <- 
 do.call(rbind,
  lapply(seq(stk), function(i) {
    out <- getSumm(stk[[i]], FUNCs = FUNCs)
    cbind(out, stk = i)
    })
 )

dat $ cfun <- factor(dat $ fun, levels = c("ssb","fbar","catch","rec"), labels = c("SSB","Mean F","Total Catch","Recruitment"))
xyplot(data ~ year | cfun, group = stk, data = dat,
      scales = list(y = list(relation = "free", rot = 0)),
      col = brewer.pal(length(stk), "Set1"),
      lty = 1, type = c("l","g"), lwd = 2,
      as.table = TRUE,
      ylab = "", xlab = "Year", between = list(x = 0.5, y = 0.5),
      alternating = FALSE)
      
      
