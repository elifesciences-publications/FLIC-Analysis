# runFlicAnalysis.R 
# Author: Edric Tam 2016/12/18 edric.tam@yale.edu
# Carlson Lab, Department of Molecular, Cellular, and Developmental Biology, Yale University

# just a simple function that reads in the data and then run the analysis. This function is created separately from 
# run.R in order to make the code more flexible and clean. If one desires to batch process files
# or add any addition custom analysis, it can be easily added here inside this function. 

runFlicAnalysis<-function(filename) {
  df<-read.csv(filename)
  makePlotsAndTables(df,filename)
}

