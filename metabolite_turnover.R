######### START ON LINE 72 ######################### - also need to add comments and save to git


## load packages, need to be installed first (but only once)
## to install isotopelabeling
## for xcms see
## http://bioconductor.org/packages/devel/bioc/html/xcms.html

# Package set up if needed. Note ProteinTurnover requires R 3.2 to be installed
# install.packages("ProteinTurnover", repos="http://www.stat.umn.edu/~rend0020/repos")

library(ProteinTurnover)
library(xcms)

# If needed, this sets the wd to where your data is located
input_path <- "C:/Users/Lab/Desktop/Coding_Bits/Turnover/data"

# Where you want the files to go
output_path <- "C:/Users/Lab/Desktop/"

setwd(input_path)

# Generate files list in same folder as data by navagating to that folder and entering the following command in the command prompt:
# dir *.mzXML /b > files.txt

# Needed function that does the official reading in of eic data.
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
  #dat$name <- factor(names(mz)[dat$channel], levels=names(mz)) #### Commented this out so it would run. Could try to tack on names somewhere
  structure(dat, class=c("EIC", class(dat)))
}

## doublecheck order for rt
setup_mzs <- function(x){
  name <- x[["name"]]
  mz <- x[["mz"]]
  i <- 6
  for(iso in 1:as.integer(x[["number.of.labels"]])){
    name <- c(name, x[[i]])
    mz <- c(mz, x[[i+1]])
    i <- i + 1
  }
  out <- as.data.frame(name, stringsAsFactors = FALSE)
  out$mz <- mz
  out
}

# Gets input from template
inputValues <- read.csv("input_template.csv", stringsAsFactors = FALSE)


# Gets names of mzXML files that will be analyzed
fileList <- read.table("files.txt", header = FALSE, stringsAsFactors = FALSE) 

# readEIC input data setup
#### Keep this constant for all analysis?
mz.tol <- .005


finalOut <- lapply(fileList$V1, function(f){
  fileOutput <- data.frame()
  for(record in 1:length(inputValues[, 1])){
    n_mz <- setup_mzs(inputValues[record, ])
    
    rt.max <- inputValues[record, "rt.max"]
    rt.min <- inputValues[record, "rt.min"]
    rt <- c(rt.min, rt.max)*60
    out <- readEIC(f, n_mz[["mz"]], mz.tol, rt) 
    
    ### ugly way to assign names. May be best to look into new way of doing this but works for now
    for(j in 1:length(out$scan)){
      if(out[j,"channel"] == 1){
        out$name[j] <- inputValues[record, "name"]
      }
      else{
        out$name[j] <- inputValues[record, (2 + 2*out$channel[j])]
      }
    }

    # gets base amino values (non labeled)
    b <- out$channel== 1
    base <- out[b, ]

    # Loop through for non bases
    for(c in 2:max(out$channel)){
      channelMatch <- out$channel == c
      labeled <- out[channelMatch, ]
      reg <- lm(labeled$intensity ~ base$intensity) ## will need to do some pretty labeling if we want these plots
      regSlope <- reg$coefficients[[2]]
      rSquared <- summary(reg)$r.squared
      
      ### this is where I could build in an if to generate a plot... if I had one!
      setwd(output_path)

      pdf(paste(substr(f, 1, nchar(f)-6), as.character(base[record, "name"]), "vs", as.character(labeled[record, "name"]), ".pdf", sep = ""))
      plot(base$intensity, labeled$intensity, xlab =  as.character(base[record, "name"]), ylab = as.character(labeled[record, "name"]), main = "Unlabeled vs Labeled")
      abline(reg)
      dev.off()
      
      setwd(input_path)
      ###
      fileOutput <- rbind(fileOutput, data.frame(f, base$name[[1]], labeled$name[[1]], rSquared, regSlope))
    }
    
  }
  colnames(fileOutput) <- c("file_name", "base", "label", "r2", "slope")
  fileOutput
})


setwd(output_path)


for(i in finalOut){
  file_name <- as.character(i[1, "file_name"])
  print(file_name)
  write.csv(i, file = paste(substr(file_name, 1, nchar(file_name) - 6), ".csv", sep = ""), row.names = FALSE)
}


# Sets up input data. This will eventually be read in from a csv to make it much smoother and automated. 
mz.tol <- .005
rt <- c(5.05, 5.37)*60 
mz.min <- mz - mz.tol
mz.max <- mz + mz.tol
file <- "DW_Dried-pos.mzXML"
# tol and rt window determined empirically. Currently generating R2 of .927 for Dried file.
# but only .6732 for Fresh

out <- readEIC(file, mz, mz.tol, rt)

divided <- split(out, out$channel == 1)

A_i <- divided$T$intensity
AN_i <- divided$F$intensity

plot(A_i, AN_i)
reg <- lm(AN_i ~ A_i)
regSlope <- reg$coefficients[[2]]
abline(reg)



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







