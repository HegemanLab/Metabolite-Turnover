##############################################################################################################
# Loading Libraries - May need to install missing packages first but will only need to do that once.
# Run everything in this section.
##############################################################################################################

# load packages, need to be installed first (but only once)
# to install isotopelabeling
# for xcms see
# http://bioconductor.org/packages/devel/bioc/html/xcms.html

# Package set up if needed. Note ProteinTurnover requires R 3.2 to be installed
# install.packages("ProteinTurnover", repos="http://www.stat.umn.edu/~rend0020/repos")
# install.packages("spatstat")

library(ProteinTurnover)
library(xcms)
library(spatstat)
library(lattice)
library(latticeExtra)



##############################################################################################################
# User Input - Check that all the values are correct and then run this secton.
##############################################################################################################

# If needed, this sets the working directory to where your data is located
# Note: if path copied from windows all \ characters will need to be switched to / characters. Use ctrl+f and replace to do this quickly for long file pathes.
inputPath <- "C:/Users/Hegeman Lab/Desktop/Data/TO_general_data"

### Where you want the output files to go. Note, each of these folders needs to be created before the script is called
# Time-Series Regressions output location
outputPathTime.SeriesRegressions <- "C:/Users/Hegeman Lab/Desktop/Output/TimeSeries/Regressions"

# Time-Series Regressions output location
outputPathTime.SeriesEICs <- "C:/Users/Hegeman Lab/Desktop/Output/TimeSeries/EICs"

# Output location for the regression data on each amino acid/label combination for each file. Stored as .csv files
outputPathRegressionTables <- "C:/Users/Hegeman Lab/Desktop/Output/TimeSeries/Tables"

# Output location for bad amino acid/label regression. This will show regression plot as a pdf. A bad regression is 
# defined as one with an r-squared value less than the designated minimum r-squared value. 
outputPathBadRegressions <- "C:/Users/Hegeman Lab/Desktop/Output/TimeSeries/Bad"
minimumRSquared <- .96

# Gets input from input csv created using input_template
# File name of input csv stored in the inputPath set as indicated above. 
inputValuesPath <- "input_template.csv"
 
# Gets list of files with time and set information. Also based on provided template. 
filesPath <- "files_template.csv" 
# Command line trick to generate files list in same folder as data -
# Navagate to the folder with the mzXML or mzML files of interest and enter the following command:
# dir *.mzXML /b > files.txt



##############################################################################################################
# Generate Needed Functions - This creates all the needed functions. Run everything in this section.
##############################################################################################################

# Needed function that does the official reading in of eic data.
# NOTE mz is a list of mzs
read_EIC <- function(file, mz, mz.tol, rt) {
  # Create raw object
  xraw <- xcmsRaw(file)
  
  # For each mz make a raw EIC and return that list.
  eic <- lapply(mz, function(x){
    rawEIC(xraw, mzrange = c(x - mz.tol, x + mz.tol), rtrange = rt)
  })
  
  # Assigns indices to eic data
  dat <- do.call(rbind, lapply(seq_along(eic), function(idx) { 
    out <- data.frame(eic[[idx]])
    out$channel <- idx
    out
  }))
  
  # Returns correctly structured data
  structure(dat, class=c("EIC", class(dat)))
}

# Function to read in compound names and masses from a record and outputs them as a dataframe.
setup_mzs <- function(record){

  # Gets set up information from the provided record. A record here is a single row from the input_template file. 
  name <- record[["name"]]
  mz <- record[["mz"]]
  label <- record[["label"]]
  mz.tol <- record[["mz.tol"]]
  
  # helper list
  names <- c(name)
  
  # i is used to add in as many masses as needed by the number of labels an input has. 
  i <- 8
  
  # Loops for each label indicated in the input_template. Creates a label name and reads the mz. Adds them to their
  # respective lists for later use. 
  for(iso in 1:as.integer(record[["number.of.labels"]])){
    names <- c(names, paste0(name, "+", as.character(iso), label))
    mz <- c(mz, record[[i]])
    i <- i + 1
  }
  
  # Put data into dataframe. Dataframe includes base amino name, labeled names, mz values, and mz.tol values.
  out <- as.data.frame(names, stringsAsFactors = FALSE)
  out$mz <- mz
  out$mz.tol <- rep(mz.tol, length(mz))
  
  # Returns output
  out
}

# Gets names of mzXML files for analysis from filesTemplate.csv or similarly structured file
setup_files <- function(input_csv_filename){
  
  # Read in data
  fileInfo <- read.csv(input_csv_filename, header = TRUE, stringsAsFactors = FALSE) 
  
  # Split data into separate lists
  fileList <- fileInfo$File
  times <- fileInfo$Time
  sets <- fileInfo$Set
  
  # Return list of lists for later use. 
  list(fileList, times, sets)
}

# Function for generating output tables and graphs
generate_output <- function(fileList){
  
  # Sets up dataframe to outside loop
  outData <- data.frame()
  
  # Unpack data from input
  fi <- fileList[[1]]
  times <- fileList[[2]]
  sets <- fileList[[3]]
  
  # loops over each base amino from the input_template
  for(record in 1:length(inputValues[, 1])){
    
    # Initialize set to make sure first loop executes
    current_set <- -1
    
    # Gets mz values and sets up needed variables.
    n_mz <- setup_mzs(inputValues[record, ])
    rt.max <- inputValues[record, "rt.max"]
    rt.min <- inputValues[record, "rt.min"]
    label <- inputValues[record, "label"]
    name <- inputValues[record, "name"]
    mz.tol <- inputValues[record, "mz.tol"]
    
    # Dataframe to hold data in the loop before being passed out to outData
    set_data <- data.frame()
    
    # RT window to seconds
    rt <- c(rt.min, rt.max)*60
    
    # Iterates over each file
    for(j in 1:length(fi)){
      
      # Checks if the loop has reached a new set or not
      # Note: This implementation means that sets must be grouped together when read in. ie all data from set one must be together in consecutive rows in csv.
      if(sets[j] != current_set){
        
        # If this isn't the first time running through the loop.
        if (!is.empty(set_data)){
          
          # Attach set_data(internal loop) to out_data(sits outside loop)
          outData <- rbind(outData, set_data)
          
          # Gets relative abundance data for time series using protein turnover package
          x <- relAbForTimes(set_data$Count, set_data$Channel, set_data$RT, set_data$time)
          
          # Sets Time and Label as factors for plotting
          x$data.long$f1 <- factor(x$data.long$TimePoint, labels = unique(x$data.long$TimePoint)) 
          x$data.long$f2 <- factor(x$data.long$Channel, labels = unique(x$data.long$Channel))
          
          browser()
          
          # Sets up pannel for graphing a linear regression line in the lattice package xyplot.
          panel.lines <- function(x, y) {
            panel.xyplot(x, y) # show points 
            panel.lmline(x, y)  # show smoothed line 
          }
         
          # Sets directory for the graph to be generated
          setwd(outputPathTime.SeriesRegressions)
          
          # eventually add in xlab.top, ylab.right. Look at panel labels or strip labels for formatting. 
          
          # Plots Time-Series Regression for a set using lattice graphics package
          trellis.device(device = "pdf", file = paste0(set_data$name[1], "-", current_set,"-regression.pdf"))
          print(useOuterStrips(xyplot(data.long$Count ~ data.long$BaseCount | data.long$f2 * data.long$f1,
                                      as.table = TRUE, data = x, panel = panel.lines,
                                      scales=list(alternating=0), 
                                      main = list(label="Time Series Regressions", cex = 1.75),
                                      xlab = list(label = "Relative Abundance", cex = 1.5), 
                                      ylab = list(label="Labeling Time", cex = 1.5),
                                      par.settings = standard.theme(color = FALSE), 
                                      xlab.top = list(label="Channel \n(1 is unlabeled)", cex = 1.5),
                                      ylab.right =list(label="Channel Relative Abundance", cex = 1.5, srt = -90)
                                      ), 
                               strip.left = strip.custom(bg = "white", style=1, horizontal = T), 
                               strip = strip.custom(bg = "white")
                               )
                )
          dev.off()
          
          # Set directory for EICs
          setwd(outputPathTime.SeriesEICs)
          
          # Plots Time-Series EICs for a set using lattice graphics package
          trellis.device(device = "pdf", file = paste0(set_data$name[1], "-", current_set,"-eic.pdf"))
          print(useOuterStrips(xyplot(data.long$Count ~ data.long$RT | data.long$f2 * data.long$f1, as.table = TRUE,
                                      strips = TRUE, data = x, scales=list(alternating=0), 
                                      main = list(label="EIC Plot", cex = 1.75),
                                      xlab = list(label="RT", cex = 1.5),
                                      ylab = list(label="Labeling Time", cex=1.5),
                                      xlab.top = list(label="Channel \n(1 is unlabeled)", cex = 1.5),
                                      ylab.right =list(label="Channel Relative Abundance", cex = 1.5, srt = -90),
                                      par.settings = standard.theme(color = FALSE) 
                                      ),
                               strip.left = strip.custom(bg = "white", style=1, horizontal = T), 
                               strip = strip.custom(bg = "white")
                               )
                )
          dev.off()
          
          # Sets directory back to input location
          setwd(inputPath)
        }
        
        # Now that the set has changed, assign the new set as the current set
        current_set = sets[j]
        
        # Clear out the set_data variable and start fresh with a new data frame
        set_data <- data.frame()
      }
      
      ## This script will likely be slow because it will recreate all the xcmsraws. Future updates may address this.
      # Get EIC data for input into relAbFromCounts function
      out <- read_EIC(fi[j], n_mz[["mz"]], mz.tol, rt) 
      
      # Gets relative abundance data and appends a time value to it
      relAb <- relAbFromCounts(out$intensity, out$channel, out$scan, norm_channel=1)
      relAb.data <- relAb$data.long
      relAb.data$time <- rep(times[j], length(relAb.data[,1]))
      
      # Assigns names. Future updates may do this more elegantly.
      # for each row (record) in relAb.data
      for(x in 1:length(relAb.data[, 1])){
        
        # if it's the base channel, give it the base amino name
        if(relAb.data[x, "Channel"] == 1){
          relAb.data$name[x] <- inputValues[record, "name"]
        }
        
        # else generate a name based on the label
        else{
          relAb.data$name[x] <- paste0(name, "+", as.character(relAb.data$Channel[x] - 1), label)
        }
      }
      
      # Append filenames and base amino to the data for future use
      relAb.data$fileName <- rep(fi[j], length(relAb.data[,1]))
      relAb.data$unlabeledAmino <- rep(name, length(relAb.data[, 1]))
      
      # bind the data from each time together
      set_data <- rbind(set_data, relAb.data)
    }
  }
  
  # Returns output. This represents all aminos in all files
  outData
}

# Takes in data generated by generate_output and writes csvs of the regression data for each file for each amino. 
# Also generates individual regression plots for each pairing of unlabeled and labeled amino below the minimumRSquared. 
write_regression_tables <- function(out, minimumRSquared){
  
  # Iterate over every unique file name in the out data. 
  for (file in unique.default(out$fileName)){
    
    # Get all data pertaining to the current file
    subsetter <- out$fileName == file
    outSubset <- out[subsetter, ]
    
    # iterates over every amino acid from input template
    for(amino in inputValues[, 1]){
      
      # Get all data pertaining to the current amino
      aminoSubsetter <- outSubset$unlabeledAmino == amino
      aminoSubset <- outSubset[aminoSubsetter, ]
      
      # Creates a dataframe to store output
      fileOutput <- data.frame()
      
      # gets base amino values (non labeled)
      b <- aminoSubset$Channel== 1
      base <- aminoSubset[b, ]
      
      # Loop through for non bases
      for(c in 2:max(aminoSubset$Channel)){
        
        # Gets data pertaining to current label
        channelMatch <- aminoSubset$Channel == c
        labeled <- aminoSubset[channelMatch, ]
       
        # generates regression model and stores slope and r-squared values
        reg <- lm(labeled$Count ~ base$Count) 
        regSlope <- reg$coefficients[[2]]
        rSquared <- summary(reg)$r.squared
        
        # checks if regression meets minimum r-squared value. If regression is bad, plots it and saves the plot.
        if(rSquared < minimumRSquared || is.na(rSquared)){
          
          # Sets working directory to correct location
          setwd(outputPathBadRegressions)
          
          # Plot and save regression
          pdf(paste(substr(file, 1, nchar(file)-6), "-", amino, "vs", as.character(labeled[1, "name"]), ".pdf", sep = ""))
          plot(base$Count, labeled$Count, xlab =  amino, ylab = as.character(labeled[1, "name"]), main = "Unlabeled vs Labeled")
          
          # Script crashes if not enough points are detected to generate a regression line so try statement used. 
          try(abline(reg))
          
          # End saving to pdf
          dev.off()
          
          # Resets input to desired path
          setwd(inputPath)
        }
        
        # Add amino vs label data to overall file data
        fileOutput <- rbind(fileOutput, data.frame(file, base$name[[1]], labeled$name[[1]], rSquared, regSlope, stringsAsFactors = FALSE))
      }
      
      # Format output data so it has column names and then write it to a csv
      setwd(outputPathRegressionTables)
      colnames(fileOutput) <- c("file_name", "base", "label", "r2", "slope")
      file_name <- paste0(substr(file, 1, nchar(file)-6), "-", amino, ".csv")
      write.csv(fileOutput ,file = file_name, row.names = FALSE)
      setwd(inputPath)
    }
  }
}



##############################################################################################################
# Run analysis - Need to step through this sequentially. write_regression_tables is optional, only if desired.
##############################################################################################################

# Sets directory to the provided input directory
setwd(inputPath)

# Read in csv from input_template.csv structured file. 
inputValues <- read.csv(inputValuesPath, stringsAsFactors = FALSE)  # Currently all aminos minus the two that didn't have retention times.

# Gets file lists
fileList <- setup_files(filesPath)

# Generate output Time-Series output
out <- generate_output(fileList)

# Generate regression tables and bad regression plots. 
write_regression_tables(out, minimumRSquared)

# Note: Error in int_abline(a = a, b = b, h = h, v = v, untf = untf, ...) : 
#         'a' and 'b' must be finite
# May be thrown. This just means a regression line could not be drawn for some of the bad regressions for some reason
# ie because they did not have enough points or some data is missing. 
