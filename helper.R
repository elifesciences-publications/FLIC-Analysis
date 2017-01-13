# helper.R 
# Author: Edric Tam 2016/12/18 edric.tam@yale.edu
# Carlson Lab, Department of Molecular, Cellular, and Developmental Biology, Yale University

# The helper.R code contains helper functions that are used in the main analysis

# load the library we used for peak detection, pracma
library(pracma)

# this function helps us find the mode of the input data
# reference: http://stackoverflow.com/questions/2547402/is-there-a-built-in-function-for-finding-the-mode
find.mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# this function helps us find candidate feeding events in the licks based on the timeTolerance parameter
# the timeTolerance parameter specifies how big of a gap in the licks are we willing to tolerate to still call 
# those licks 1 feeding event
# if you increae timeTolerance, # of feeding events should decrease, and vice versa
findCandidateEvents<-function(licks, timeTolerance) {
  beginnings<-c(licks[1],licks[which(diff(licks)>timeTolerance)+1])
  endings<-c(licks[which(diff(licks)>timeTolerance)], licks[length(licks)])
  durationsCandidate<-endings - beginnings + 1
  candidateEvents<-data.frame(beginnings, endings, durationsCandidate)
  return(candidateEvents)
}

# this function take the candidateEvents and peakData as inputs and try to return the duration of the 
# feeding events and the time between feeding events as outputs. 
# Basically, the findCandidateEvents function gives us a list of candidate events based on a temporal parameter, 
# and the pracma peak detection gives us a list of peaks in the data based on a peak detection method, which can
# be seen as another list of candidate feeding events. 
# We combine both lists to return the finalized durations of feeding events and time between feeding events. 
findDurationsAndTimeBetween<-function(peakData,candidateEvents){
  
  durations<-rep(NA, nrow(peakData))
  timeBetween<-rep(NA, length(durations)-1)
  realEventsIndex<-rep(NA,length(durations))
  for (i in 1:length(peakData$peakLocation)) {
    peak = peakData$peakLocation[i]
    for (j in 1:length(candidateEvents$beginnings)){
      if(candidateEvents$beginnings[j]<= peak && candidateEvents$endings[j] >= peak){
        durations[i] <- candidateEvents$durations[j]
        realEventsIndex[i]<-j
      }
    }
  }
  
  timeBetween<-candidateEvents$beginnings[realEventsIndex[2:length(realEventsIndex)]]-candidateEvents$endings[realEventsIndex[1:(length(realEventsIndex)-1)]]
  timeBetween[which(timeBetween<0)]<-diff(peakData$peakLocation)[which(timeBetween<0)]
  if(nrow(peakData)==1) {
    timeBetween<-c(NA)
  }
  return(list("durations" = durations,"timeBetween" = timeBetween))
}