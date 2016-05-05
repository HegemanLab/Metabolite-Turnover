# Metabolite-Turnover

## Introduction/Motivation
This script was developed for use in analyzing plant metabolomics. This particular script was originally utilized for analyzing amino acids. The script utilizes the XCMS package in R developed by Bioconductor. Ultimately this script takes in a table of files (potentially from a time series) and a table of compounds of interest that have been labeled. The script then looks through the files for the compound and each labeled version of the compound and generates time-series regression plots, time-series EIC plots. There is also the option to output individual regression plots for file/compound combinations that produce poor regression results, and also outputs tables for each compound in each file.

## Example output
Each of the following files shows a portion of the output generated from originally running this script. Note the formatting of these output files may have been changed in the script, but the data presented is the same. These example output are not necessarily from the same files or times. 
- bad_regression_example.pdf
- EIC_example.pdf
- Regression_example.pdf
- table_example.csv

## Usage Notes
All lines that require editing are between line 25 and 55 in the script. 
### Lines that require editing:
- 29: set your path to your data. This location should be where your input files are as well as your data files. 
- 33, 36, 39, 43: File paths to folders where various output files will be sent.  
- 44: A minimum R^2 value that determines which individual plots will be output as "bad regressions".
- 48: Name of file holding filled in input data (data regarding compounds of interest and labels).
- 51: Name of .csv file holding data regarding the files to be processed. 

### Workflow
- Fill in Input_Data_Template.csv and files_template.csv.
- Update file name's in script. 
- If needed, create output directories on your computer.
- Update all file pathes and input file names. 
- Install needed packages then run lines 15-19 to activate libraries.
- Run 29-51 to set up input variables.
- Run 64-352 to read in needed functions.
- Run 359 set input directory
- Run 362 to get data from input templates.
- Run 365 to get file data
- Run 368 to run full analysis and generate time-series plots. Depending on the number of files and compounds this may take quite some time.
- (Optional) Run 371 to write the output data to csv files and to output bad regressions. 

## Script and Input Templates
A brief description of each file in this repo and its role. 

### updated_metabolite_turnover.R
R script that does all the analysis. 

### input_template.csv
Input template that holds all input data about the compounds that are being analyzed. Note, the example data in the template should be replaced with your own data and saved. Also, the columns l1.mz, l2.mz ... represent the masses of the first labeled compound, the second labeled compound... etc. Although only 4 columns are shown for labels, these can be extended indefinitely to the right of l4.mz as long as it is reflected in the "number.of.labels" columns. 

### files_template.csv
Template of how information about the files should be entered and example data of how information should be entered. Each set must be listed together although not neccessarily in time order. Example data should be replaced with your own data but column names should remain unchanged. 

## Hegeman Lab - University of Minnesota Twin-Cities
This code was developed for use in the Hegeman Lab at the University of Minnesota Twin-Cities. If you use this script in your research, please don't forget to site us. Additionally, if there are any questions about how to use this code, feel free to contact the lab or the scripts creater directly through GitHub. 
