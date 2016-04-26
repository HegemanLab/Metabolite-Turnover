## load packages, need to be installed first (but only once)
## to install isotopelabeling
## for xcms see
## http://bioconductor.org/packages/devel/bioc/html/xcms.html

# Package set up if needed. Note ProteinTurnover requires R 3.2 to be installed
# install.packages("ProteinTurnover", repos="http://www.stat.umn.edu/~rend0020/repos")
# install.packages("spatstat")

library(ProteinTurnover)
library(xcms)
library(spatstat)
library(lattice)

source("http://bioconductor.org/biocLite.R")
biocLite("xcms")  

# If needed, this sets the wd to where your data is located
# Note: if path copied from windows all \ characters will need to be switched to / characters. Use ctrl+f and replace to do this quickly for long file pathes.
inputPath <- "C:/Users/Hegeman Lab/Desktop/Data/TO_general_data"

# Where you want the files to go
outputPath <- "C:/Users/Hegeman Lab/Desktop/Output/TimeSeries"

setwd(inputPath)

# Needed function that does the official reading in of eic data.
# NOTE mz is a list of mzs
read_EIC <- function(file, mz, mz.tol, rt) {
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
  
  # adds in names. Not being used by this script but from original.
  #dat$name <- factor(names(mz)[dat$channel], levels=names(mz)) #### Commented this out so it would run. Could try to tack on names somewhere
  
  structure(dat, class=c("EIC", class(dat)))
}

## Reads in compound names and masses from a record and outputs them as a dataframe.
setup_mzs <- function(x){
  
  name <- x[["name"]]
  names <- c(name)
  mz <- x[["mz"]]
  label <- x[["label"]]
  mz.tol <- x[["mz.tol"]]
  i <- 8
  for(iso in 1:as.integer(x[["number.of.labels"]])){
    
    names <- c(names, paste0(name, "+", as.character(iso), label))
    mz <- c(mz, x[[i]])
    i <- i + 1
  }
  out <- as.data.frame(names, stringsAsFactors = FALSE)
  out$mz <- mz
  out$mz.tol <- rep(mz.tol, length(mz))
  out
}

# Gets input from template
inputValues <- read.csv("input_template.csv", stringsAsFactors = FALSE)  # Currently all aminos minus the two that didn't have retention times. 
filesInput <- "files_template.csv" # Currently subset of time series data

# Gets names of mzXML files that will be analyzed
setup_files <- function(input_csv_filename){
  fileInfo <- read.csv(input_csv_filename, header = TRUE, stringsAsFactors = FALSE) 
  fileList <- fileInfo$File
  times <- fileInfo$Time
  sets <- fileInfo$Set
  
  list(fileList, times, sets)
}

fileList <- setup_files(filesInput)

# Function for generating output tables and graphs
generate_output <- function(fileList){
  
  outData <- data.frame()
  
  for(record in 1:length(inputValues[, 1])){
    
    # Possibly move this outside the record for loop. I think that'll work and will save a little run time.
    # Unpack data from input
    fi <- fileList[[1]]
    times <- fileList[[2]]
    sets <- fileList[[3]]
    
    current_set <- -1
    
    # Gets mz values and sets up needed variables.
    n_mz <- setup_mzs(inputValues[record, ])
    rt.max <- inputValues[record, "rt.max"]
    rt.min <- inputValues[record, "rt.min"]
    label <- inputValues[record, "label"]
    name <- inputValues[record, "name"]
    mz.tol <- inputValues[record, "mz.tol"]
    set_data <- data.frame()
    
    # RT window to seconds
    rt <- c(rt.min, rt.max)*60
    
    for(j in 1:length(fi)){
      
      #### Note: This implementation means that sets must be grouped together when read in. ie all data from set one must be together in consecutive rows in csv.
      if(sets[j] != current_set){
        # Fence post. If this isn't the first time running through the loop.
        if (!is.empty(set_data)){
          outData <- rbind(outData, set_data)

          
          x <- relAbForTimes(set_data$Count, set_data$Channel, set_data$RT, set_data$time)
          x$data.long$f1 <- factor(x$data.long$TimePoint, labels = unique(x$data.long$TimePoint)) 
          x$data.long$f2 <- factor(x$data.long$Channel, labels = unique(x$data.long$Channel))
          
          # Kinda the EIC plots... plot(x$data.long$RT, x$data.long$Count)
          
          panel.lines <- function(x, y) {
            panel.xyplot(x, y) # show points 
            panel.lmline(x, y)  # show smoothed line 
            
          }
          # eventually add in xlab.top, ylab.right. Look at panel labels or strip labels for formatting. 
          
          setwd(outputPath)
          
          trellis.device(device = "pdf", file = paste0(set_data$name[1], "-", current_set,"-regression.pdf"))
          print(xyplot(data.long$Count ~ data.long$BaseCount | data.long$f2 * data.long$f1, as.table = TRUE, strips = TRUE, data = x, panel = panel.lines, scales=list(alternating=0), main = "Time Series Regressions", xlab = "Relative Abundance", ylab = "Labeling Time"))
          dev.off()
          
          trellis.device(device = "pdf", file = paste0(set_data$name[1], "-", current_set,"-eic.pdf"))
          print(xyplot(data.long$Count ~ data.long$RT | data.long$f2 * data.long$f1, as.table = TRUE, strips = TRUE, data = x, scales=list(alternating=0), main = "EIC Plot", xlab = "RT", ylab = "Labeling Time"))
          dev.off()
          
          setwd(inputPath)
          
          # Seems to generate same plot. First one is coming out a bit wonky.
          # xyplot(Count ~ RT | f2 * f1, as.table = TRUE, strips = TRUE, data = set_data, scales=list(alternating=0), main = "EIC Plot", xlab = "RT", ylab = "Labeling Time")
        }
        
        
        
        
        #plot_time-series(set_data)
        
        current_set = sets[j]
        set_data <- data.frame()
      }
      
      ## This cunction will likely be slow because it will recreate all the xcmsraws. May be a way to do this once as a set up step then just find the correct xraw in a table or something (dynamic programming ish)
      # Get EIC data for input into relAbFromCounts function
      out <- read_EIC(fi[j], n_mz[["mz"]], mz.tol, rt) 
      
      # Gets relative abundance data and appends a time value to it
      relAb <- relAbFromCounts(out$intensity, out$channel, out$scan, norm_channel=1)
      relAb.data <- relAb$data.long
     
      # only for debugging, remove later
      timej <- times[j]
      
      relAb.data$time <- rep(times[j], length(relAb.data[,1]))
      
      ### ugly way to assign names. May be best to look into new way of doing this but works for now
      for(x in 1:length(relAb.data[, 1])){
        if(relAb.data[x, "Channel"] == 1){
          relAb.data$name[x] <- inputValues[record, "name"]
        }
        else{
          relAb.data$name[x] <- paste0(name, "+", as.character(relAb.data$Channel[x] - 1), label)
        }
      }

      
      relAb.data$fileName <- rep(fi[j], length(relAb.data[,1]))
      relAb.data$unlabeledAmino <- rep(name, length(relAb.data[, 1]))
      
      # bind the data from each time together
      set_data <- rbind(set_data, relAb.data)
      
    }
  }
  
  # Returns output
  outData
  
}

# Takes in data generated by generate_output and writes csvs of the regression data for each file for each amino. 
# Also generates individual regression plots for each pairing of unlabeled and labeled amino. 
write_regression_tables <- function(out){
  
  # Iterate over every unique file name in the out data. 
  for (file in unique.default(out$fileName)){
    
    subsetter <- out$fileName == file
    outSubset <- out[subsetter, ]
    
    # iterates over every amino acid from input template
    for(amino in inputValues[, 1]){
      # Create df for output
      
      aminoSubsetter <- outSubset$unlabeledAmino == amino
      aminoSubset <- outSubset[aminoSubsetter, ]
      
      fileOutput <- data.frame()
      
      # gets base amino values (non labeled)
      b <- aminoSubset$Channel== 1
      base <- aminoSubset[b, ]
      
      # Loop through for non bases
      for(c in 2:max(aminoSubset$Channel)){
        channelMatch <- aminoSubset$Channel == c
        labeled <- aminoSubset[channelMatch, ]
        
        
        ########### Should this be RelAb (basecount~count) or count unlabeled vs count labeled, or something else. 
        
        
        reg <- lm(labeled$Count ~ base$Count) 
        regSlope <- reg$coefficients[[2]]
        rSquared <- summary(reg)$r.squared
        
        ### this is where I could build in an if to generate a plot... if I had one!
        pdf(paste(substr(file, 1, nchar(file)-6), "-", amino, "vs", as.character(labeled[1, "name"]), ".pdf", sep = ""))
        plot(base$Count, labeled$Count, xlab =  amino, ylab = as.character(labeled[1, "name"]), main = "Unlabeled vs Labeled")
       
        # Script crashes if not enough points are detected to generate a regression line so try statement used. 
        try(abline(reg))
        dev.off()
        ###
        fileOutput <- rbind(fileOutput, data.frame(file, base$name[[1]], labeled$name[[1]], rSquared, regSlope, stringsAsFactors = FALSE))
      }
      
      colnames(fileOutput) <- c("file_name", "base", "label", "r2", "slope")
      file_name <- paste0(substr(file, 1, nchar(file)-6), "-", amino, ".csv")
      write.csv(fileOutput ,file = file_name, row.names = FALSE)
    }
  }
}

# Tables stored in out. Graphics generated within the generate_output script and also sent to outputPath location.
setwd(inputPath)

out <- generate_output(fileList)

# Send output to the correct place
outputPath <- "C:/Users/Hegeman Lab/Desktop/Output/Tables"
setwd(outputPath)

# Generate output
write_regression_tables(out)


### Extras. Potential useful tools for preparing data. 
# Generate files list in same folder as data by navagating to that folder and entering the following command in the command prompt:
# dir *.mzXML /b > files.txt





