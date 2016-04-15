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
  
  # adds in names. Not being used by this script but from original.
  #dat$name <- factor(names(mz)[dat$channel], levels=names(mz)) #### Commented this out so it would run. Could try to tack on names somewhere
  
  structure(dat, class=c("EIC", class(dat)))
}

## Reads in compound names and masses from a record and outputs them as a dataframe.
setup_mzs <- function(x){

  name <- x[["name"]]
  mz <- x[["mz"]]
  label <- x[["label"]]
  mz.tol <- x[["mz.tol"]]
  i <- 8
  for(iso in 1:as.integer(x[["number.of.labels"]])){
    name <- c(name, paste0(name, "+", as.character(iso), label))
    mz <- c(mz, x[[i]])
    i <- i + 1
  }
  out <- as.data.frame(name, stringsAsFactors = FALSE)
  out$mz <- mz
  out$mz.tol <- rep(mz.tol, length(mz))
  out
}

# Gets input from template
inputValues <- read.csv("input_template.csv", stringsAsFactors = FALSE)
filesInput <- "files_template.csv"

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
generateOutput <- function(fileList){
  
  for(record in 1:length(inputValues[, 1])){
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
    
    # RT window to seconds
    rt <- c(rt.min, rt.max)*60
    
    for(j in 1:length(fi)){
      
      #### Note: This implementation means that sets must be grouped together when read in. ie all data from set one must be together in consecutive rows in csv.
      if(sets[j] != current_set){
        # Generate plots and csv readouts of  here
        
        browser()
        
        # check if current set == -1 or if set_data == null or something
        #write_amino_tables(set_data)
        #plot_time-series(set_data)
        
        browser()
        
        current_set = sets[j]
        set_data <- data.frame()
      }
      
      ## This cunction will likely be slow because it will recreate all the xcmsraws. May be a way to do this once as a set up step then just find the correct xraw in a table or something (dynamic programming ish)
      # Get EIC data for input into relAbFromCounts function
      out <- readEIC(fi[j], n_mz[["mz"]], mz.tol, rt) 
      
      # Gets relative abundance data and appends a time value to it
      relAb <- relAbFromCounts(out$intensity, out$channel, out$scan, norm_channel=1)
      relAb.data <- relAb$data.long
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
      
      browser()
      
      # bind the data from each time together
      set_data <- rbind(set_data, relAb.data)
      
    }
  }
  
 
  
  
  
#   
#   finalOut <- lapply(fi, function(f){
#     fileOutput <- data.frame()
#     for(record in 1:length(inputValues[, 1])){
#       
#       # Gets mz values and sets up needed variables.
#       n_mz <- setup_mzs(inputValues[record, ])
#       rt.max <- inputValues[record, "rt.max"]
#       rt.min <- inputValues[record, "rt.min"]
#       
#       # RT window to seconds
#       rt <- c(rt.min, rt.max)*60
#       
#       browser()
#       
#       ########### Do this loop for each time in each set (nested)
#       
#       
#       # pulls out values needed for relAb from each file
#       out <- readEIC(f, n_mz[["mz"]], mz.tol, rt) 
# 
#       relAb <- relAbFromCounts(out$intensity, out$channel, out$scan, norm_channel=1)
#       relAb.data <- relAb$data.long
#       relAb.data$time <- rep(times[1], length(relAb.data[,1]))
#       
#       
#       ### ugly way to assign names. May be best to look into new way of doing this but works for now
#       for(j in 1:length(out$scan)){
#         if(out[j,"channel"] == 1){
#           out$name[j] <- inputValues[record, "name"]
#         }
#         else{
#           out$name[j] <- inputValues[record, (2 + 2*out$channel[j])]
#         }
#       }
#       
#       # gets base amino values (non labeled)
#       b <- out$channel== 1
#       base <- out[b, ]
#       
#       # Loop through for non bases
#       for(c in 2:max(out$channel)){
#         channelMatch <- out$channel == c
#         labeled <- out[channelMatch, ]
#         reg <- lm(labeled$intensity ~ base$intensity) ## will need to do some pretty labeling if we want these plots
#         regSlope <- reg$coefficients[[2]]
#         rSquared <- summary(reg)$r.squared
#         
#         ### this is where I could build in an if to generate a plot... if I had one!
#         setwd(output_path)
#         
#         pdf(paste(substr(f, 1, nchar(f)-6), as.character(base[record, "name"]), "vs", as.character(labeled[record, "name"]), ".pdf", sep = ""))
#         plot(base$intensity, labeled$intensity, xlab =  as.character(base[record, "name"]), ylab = as.character(labeled[record, "name"]), main = "Unlabeled vs Labeled")
#         abline(reg)
#         dev.off()
#         
#         setwd(input_path)
#         ###
#         fileOutput <- rbind(fileOutput, data.frame(f, base$name[[1]], labeled$name[[1]], rSquared, regSlope))
#       }
#       
#     }
#     colnames(fileOutput) <- c("file_name", "base", "label", "r2", "slope")
#     fileOutput
#   })
}

# Tables stored in out. Graphics generated within the generateOutput script and also sent to output_path location.
out <- generateOutput(fileList)

# Send output to the correct place
setwd(output_path)

# Generate output
for(i in out){
  file_name <- as.character(i[1, "file_name"])
  print(file_name)
  write.csv(i, file = paste(substr(file_name, 1, nchar(file_name) - 6), ".csv", sep = ""), row.names = FALSE)
}


### Extras. Potential useful tools for preparing data. 
# Generate files list in same folder as data by navagating to that folder and entering the following command in the command prompt:
# dir *.mzXML /b > files.txt






