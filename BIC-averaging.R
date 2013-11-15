# BIC model averaging
#

weightedModelAverage <- function(mods, stock, FUN = AIC, nsim = 1000)
{

  FUN <- match.fun(FUN)

  # calculate weights
  ICs <- -1 * sapply(mods, FUN)
  eICs <- exp( 0.5 * (ICs - max(ICs)))
  weights <- eICs / sum(eICs)

  wt.table <- data.frame("weight (perc)" = round(weights * 100, 3))
  rownames(wt.table) <- names(mods)
  
  message("model weights are \n\t", paste(capture.output(wt.table), collapse = "\n\t"))

  # now sample from each model 1000 times and randomly select those to combine
  mod.sim <- sample(seq_along(mods), nsim, prob = weights, replace = TRUE)

  stock.sim <- propagate(stock, nsim)
  for (i in seq_along(mods)) {
    if (sum(mod.sim == i)) stock.sim[,,,,,mod.sim == i] <- stock.sim[,,,,,mod.sim == i] + mods[[i]]
  }
  # DONE !!

  stock.sim
}


## 
# this script loads the FLa4a library and writes data to a directory "test"
##
options(width=110)

library(FLa4a)

# get data and add some bits
data(ple4)
data(ple4.indices)

# fit two of the simplest models with a small amount of data
# for now we only deal with 1 area, 1 unit.
stock <- ple4
indices <- ple4.indices

# we need the covariance for proposals / simulations so use fit = "assessment"
f1 <- a4a(~ te(age, year, k = c(4, 25)), 
          list(~ s(age, k = 4), ~ s(age, k = 4), ~ age),
          stock = stock, indices = indices, fit = "assessment")
          
f2 <- a4a(~ te(age, year, k = c(4, 25)), 
          list(~ te(age, year, k = c(4, 3)), ~ te(age, year, k = c(4, 3)), ~ s(year, k = 3, by = age)),
          stock = stock, indices = indices, fit = "assessment")
 
 
# simulate
pstock <- propagate(stock, 1000)
fit1 <- pstock + f1
fit2 <- pstock + f2

stock.sim <- weightedModelAverage(list(f1, f2), stock, BIC, nsim = 1000)
 

getxy <- function(fit, p = 0.95, n = 100, yr = paste(dims(fit) $ maxyear))
{
  x <- c(fbar(fit)[,yr])
  y <- c(ssb(fit)[,yr])
  
  xy <- kde2d(x,y, n = n)
  area <- diff(xy $ x[1:2]) * diff(xy $ y[1:2])

  p2d <- function(v) (sum((xy $ z * area)[log(xy $ z) > v]) - p)^2
  logz <- optimise(p2d, c(-200, 0)) $ minimum

  xy $ z[log(xy $ z) < logz] <- NA
  xy
}

getxyPoly <- function(fit, p = 0.95, n = 100, yr = paste(dims(fit) $ maxyear))
{
  xy <- getxy(fit, p = p, n = n, yr = yr)
  
  dat <- data.frame(x = rep(xy $ x, n),
                    y = rep(xy $ y, each = n),
                    z = c(xy $ z))
  dat <- dat[!is.na(dat $ z),]
  dat[chull(dat $ x, dat $ y),c("x","y")]                  
}


# plot 1

xy1 <- getxy(fit1)
xy2 <- getxy(fit2)
xyMA <- getxy(stock.sim)

xlim = c(0, 0.5)
ylim = c(0, 5e5)
 

cols <- paste0(brewer.pal(3, "Set1"), "80")

plot(0,0,ylim = ylim, xlim = xlim, ylab = "SSB", xlab = "Fbar", type = "n")
image(xyMA, col = cols[3], add = TRUE)
image(xy1, col = cols[1], add = TRUE)
image(xy2, col = cols[1], add = TRUE)

# plot 2

ps <- lapply(list(stock.sim, fit1, fit2), getxyPoly)

lcols <- paste0(brewer.pal(3, "Set1"), "FF")
fcols <- paste0(brewer.pal(3, "Set1"), c("88", rep("44", 2)))

plot(0,0,ylim = ylim, xlim = xlim, ylab = "SSB", xlab = "Fbar", type = "n")
for (i in seq(ps))
{
  polygon(ps[[i]], col = fcols[i], border = lcols[i], lwd = 2)
}



#  polygon(dat[chull(dat $ x, dat $ y),], col = "blue")
#  dat <- dat[log(dat $ z) > logz,]
#  polygon(dat[chull(dat $ x, dat $ y),], col = "red")


