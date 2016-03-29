## load packages, need to be installed first (but only once)
## to install isotopelabeling
## for xcms see
## http://bioconductor.org/packages/devel/bioc/html/xcms.html

# Package set up if needed. Note ProteinTurnover requires R 3.2 to be installed
# install.packages("ProteinTurnover", repos="http://www.stat.umn.edu/~rend0020/repos")

library(ProteinTurnover)
library(xcms)

# If needed, this sets the wd to where your data is located
setwd("C:/Users/Lab/Desktop/Coding_Bits/Turnover/data")

# Gets input from template
inputValues <- read.csv("input_template.csv")

# Generate files list in same folder as data by navagating to that folder and entering the following command in the command prompt:
# dir *.mzXML /b > files.txt

# Gets names of mzXML files that will be analyzed
fileList <- read.table("files.txt", header = FALSE) 

# reformats data as dataframe with each amino name as the row identifier
values <- as.data.frame(inputValues)
values <- values[,-1]
rownames(values) <- input[,1]



############## START HERE. Figure out how to duplicate data read in from new template. Then just need to develop loop and add code for saving output.
# Will need to figure out how to read in names (if needed) from CVS. Should just have to use them as col names
mz <- c("normal" = 90.0550, "isotope" = 91.0520)


# Sets up input data. This will eventually be read in from a csv to make it much smoother and automated. 
mz.tol <- .005
rt <- c(5.10, 5.55)*60 
rt.min <- rep(5.10, length(mz))*60
rt.max <- rep(5.55, length(mz))*60
mz.min <- mz - mz.tol
mz.max <- mz + mz.tol

tester_raw <- xcmsRaw(file)

out <- readEIC(file, mz, mz.tol, rt)

divided <- split(out, out$channel == 1)

A_i <- divided$T$intensity
AN_i <- divided$F$intensity

plot(A_i, AN_i)
reg <- lm(AN_i ~ A_i)
regSlope <- reg$coefficients[[2]]
abline(reg)
summary(reg)
regSummary <- summary(reg)
rSquared <- regSummary$r.squared



# plotEIC ranges must be min to max, not max to min
a_plot <- plotEIC(tester_raw, mzrange = c(90.0545, 90.0555), rtrange = c(5.10, 5.55)*60)
labeled_plot <- plotEIC(tester_raw, mzrange = c(91.047, 91.057), rtrange = c(5.10, 5.55)*60)

a_eic <- rawEIC(tester_raw, mzrange = c(90.0545, 90.0555), rtrange = c(5.10, 5.55)*60)
labeled_eic <- rawEIC(tester_raw, mzrange = c(91.047, 91.057), rtrange = c(5.10, 5.55)*60)

# Possible EIC optiosn
plot.new()
lines(a_eic$scan, a_eic$intensity)
show(a_eic)

args(plotEIC)


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







