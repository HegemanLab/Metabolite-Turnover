## load packages, need to be installed first (but only once)
## to install isotopelabeling
## > install.packages("isotopelabeling", repos="http://www.stat.umn.edu/~rend0020/repos")
## for xcms see
## http://bioconductor.org/packages/devel/bioc/html/xcms.html

# Package set up if needed. Note ProteinTurnover requires R 3.2 to be installed
install.packages("ProteinTurnover", repos="http://www.stat.umn.edu/~rend0020/repos")

library(ProteinTurnover)
library(xcms)

# If needed
setwd("C:/Users/Lab/Desktop/Coding_Bits/Metabolite_turnover")

# Reads in the files that will need to be operated on. 
file <- "S2_1_10D.mzML"
# file <- files[1,1]


# Will need to figure out how to read in names (if needed) from CVS. Should just have to use them as col names
mz <- c("normal" = 90.0550, "isodope" = 91.0520)
mz.l <- 91.0520

mz.tol <- .001
rt <- c(11, 11.34)*60

# NOTE mz is a list of mzs
readEIC <- function(file, mz, mz.tol, rt) {
  # Create raw object
  xraw <- xcmsRaw(file)
  
  # For each mz make a raw EIC and return that list.
  eic <- lapply(mz, function(x){
    rawEIC(xraw, mzrange = c(x - mz.tol, x + mz.tol), rtrange = rt)
  })
  
  dat <- do.call(rbind, lapply(seq_along(eic), function(idx) { # Assigns indices to eic data
    out <- data.frame(eic[[idx]])
    out$channel <- idx
    out
  }))
  
  # adds in names 
  dat$name <- factor(names(mz)[dat$channel], levels=names(mz)) 
  structure(dat, class=c("EIC", class(dat)))
  
}

out <- readEIC(file, mz, mz.tol, rt)

# Only non used function, figure out what it do ###################
getAb <- function(out) {
  cbind(out[[1]]$prop[,c(1,2)], do.call(cbind, lapply(out, function(x) x$prop$Ab)))
}

# only uses single mz value
dofile <- function(file, mz, mz.tol, rt) {
  eic <- readEIC(file, mz, mz.tol, rt)
  relAb <- relAbFromCounts(eic$intensity, eic$channel, eic$scan, norm_channel=1)
  prop <- getProp(relAb, mz)
  list(eic=eic, relAb=relAb, prop=prop)
}

# Intermediate function
getProp <- function(out, mz) {
  prop <- out$data
  ret <- cbind(Name=factor(names(mz)[prop$Channel], levels=names(mz)), prop)
  ret
}

split_data <- split(out, f = out$name == 'normal')

# xyplot appears to no longer be a thing ALSO NOT CLEAR WHAT THIS IS SUPPOSED TO PLOT AT THIS POINT
plot.EIC <- function(x) {
  plot(intensity ~ scan | name, data=x, type="l", as.table=TRUE)
}

# Random useful bits
plot(split_data$'TRUE'$scan, split_data$'TRUE'$intensity)
rmodel <- lm(split_data$'TRUE'$intensity ~ split_data$'TRUE'$scan)
abline(rmodel)
summary(rmodel)

plot.EIC(current_data)

# To get r squared values
summary(regr)$r.squared







