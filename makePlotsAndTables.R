# makePlotsAndTables.R
# Author: Edric Tam 2016/12/18 edric.tam@yale.edu
# Carlson Lab, Department of Molecular, Cellular, and Developmental Biology, Yale University
#
# this function does the actual analysis and then output the results to several files:
#
# outputAnalysis.csv contains the values of all the tracked variables for each well
# boxplot.pdf shows boxplots for each tracked variable
# peak_plots.pdf shows raw data traces along with detected peaks
# Statistical Test Results.txt shows the results of the t-test and mann-whitney test 
#
# * Recall that the statistcal tests here are just for the users' benefit to quickly get a rough idea of
# the statistical significance of their data. For publication/other purposes, we recommend using your custom
# statistical analysis that would be suitable for your data and experiments. 
#
# There are several tunable parameters below: 
# 1. timeTolerance (default = 5)
# the time tolerance parameter tell us how far apart each sip/lick/contact have 
# to be in order for us to consider it a separate feeding event
# 2. peakTolerance (default = 70)
# this is a parameter for the pracma findpeaks function's minpeakDistance parameter
# 3. min.peak.height (default = 20)
# this is a parameter for the pracma findpeaks function's minpeakDistance parameter
# 4. baseline.offset.for.contact.threshold (default = 20)
# this is a parameter for determining the threshold of licks. 
# 5. begin.time (default = 1500)
# this is a parameter for determining when to start the analysis
# Usually, depending on the experimenter, it takes several minutes to load the flies into the FLIC wells. 
# By setting the default begin.time as 1500 (5 minutes), we allow 5 minutes for the experimenter to load in the flies 
# before we actually analyze the data. 
# 6. end.time (default = 6000) 
# this is a parameter for determining when to stop the analysis. Depending on your experiment, 
# you might want to analyze the first 15 minutes of data or you might want to to analyze the full length of the 
# data you collected. The default is 18000. 
# Remember, the time unit here is in samples. We used a default sampling rate of 5 samples/sec,
# hence 1500/5 = 300 seconds = 5 minutes and 18000/5 = 3600 = 60 minutes 

makePlotsAndTables<-function(df,filename) {
  
  # create dataframe to hold data that we will write to csv file
  numWell <- 12
  
  # set default tunable parameters
  # set the baseline offset
  # basically, our code computes the baseline by finding the mode of the data
  # we then add an offset to the baseline, and this will be the threshold 
  # for detecting sips/licks/contacts, whichever way you want to call it
  # for example, if the baseline offset is 20, and your data's is computed 
  # to be baseline is 10, then the threshold for sips/licks/contacts would be 30
  
  baseline.offset.for.contact.threshold <- 20
  
  # other tunable parameters as explained above are here
  
  begin.time = 1500
  end.time = 18000
  timeTolerance = 5
  peakTolerance = 70
  min.peak.height = 20
  
  # create array for each tracked variable to hold the analysis results
  
  numLicks <- rep(NA,numWell)
  numEvents<- rep(NA,numWell)
  meanDuration<- rep(NA,numWell)
  medianDuration<- rep(NA,numWell)
  meanTimeBetween<- rep(NA,numWell)
  medianTimeBetween<- rep(NA,numWell)
  meanIntensity<- rep(NA,numWell)
  medianIntensity<- rep(NA,numWell)
  
  
  # create the raw data traces plots overlayed with detected peaks
  pdf('peak_plots.pdf')
  durationList = list()
  # loop over each of the 12 wells
  for (i in 1:12) {
    # print out progress
    print(paste('Working on ', toString(i), 'out of 12 wells', sep= ' '))
    
    # +4 due to the df having time, samples and other columns (totalling 4 columns) before the 12 actual data columns
    datum<-df[begin.time:end.time,i+4]
    
    # set threshold for sip/contact/licks
    data.baseline <- find.mode(datum)
    contactThreshold <- data.baseline + baseline.offset.for.contact.threshold
    
    # find peaks (feeding events) in the data
    # there is a chance that there are no peaks in data, then peakData would be an empty dataframe
    # keep that in mind
    peakData<-findpeaks(datum, nups = 0, ndown = 0, zero = "0", peakpat = NULL,
                        minpeakheight = min.peak.height, minpeakdistance = peakTolerance,
                        threshold = 0, npeaks = 0, sortstr = FALSE)
    
    
    
    if (is.matrix(peakData)) {
      peakData<-as.data.frame(peakData)
      colnames(peakData)<-c('peakIntensity','peakLocation','peakBegin','peakEnd')
      peakData<-peakData[which(!duplicated(peakData$peakLocation)),]
    } else if (is.numeric(peakData)) {
      peakData<-data.frame('peakIntensity' = peakData[1],'peakLocation'= peakData[2],'peakBegin'= peakData[3],'peakEnd'= peakData[4])
    }
    if(!is.null(peakData)) {
      peakData<-peakData[order(peakData$peakLocation),]
    }
    # define licks as data that is above contactThreshold
    if(length(which(datum>contactThreshold)) > 0) {
      licks <- which(datum>contactThreshold)
    } else {
      licks <-c(0)
    }
    if(!is.null(peakData) && (length(which(datum>contactThreshold)) > 0)) {
      candidateEvents<-findCandidateEvents(licks, timeTolerance)
      
      durationsAndTimeBetween<-findDurationsAndTimeBetween(peakData,candidateEvents)
      durations<-unlist(durationsAndTimeBetween$durations)
      timeBetween<-unlist(durationsAndTimeBetween$timeBetween)
    } else {
      candidateEvents<-c(NA)
      durationsAndTimeBetween<-c(NA)
      durations<-c(NA)
      timeBetween<-c(NA)
    }
    # make plot
    par(mar=c(4.5,4.5,2.5,2.5))
    plot(seq_along(as.vector(datum)), as.vector(datum), type="l",main=colnames(df)[i+4], xlab="Time Indices (5/s)", ylab="Peak Intensity")#ylim=c(min(peakIntensity),max(peakIntensity)))
    points(peakData$peakLocation, peakData$peakIntensity, col="blue", pch=20)
    
    
    numLicks[i]<- length(which(datum>contactThreshold))
    if(!is.null(peakData)) {
      numEvents[i]<- nrow(peakData)
      durations[which(durations == 0)]<-NA
      meanDuration[i]<-mean(na.omit(durations)) 
      medianDuration[i]<-median(na.omit(durations))
      meanTimeBetween[i]<-mean(na.omit(timeBetween))
      medianTimeBetween[i]<-median(na.omit(timeBetween))
      meanIntensity[i]<-mean(peakData$peakIntensity) # might put unique in front?
      medianIntensity[i]<-median(peakData$peakIntensity)
    } else {
      numEvents[i] <- 0
      durations[which(durations == 0)]<-NA
      meanDuration[i]<-NA
      medianDuration[i]<-NA
      meanTimeBetween[i]<-NA
      medianTimeBetween[i]<-NA
      meanIntensity[i]<-NA
      medianIntensity[i]<-NA
    }
    
  durationList[[i]] <-durations
  }
  
  dev.off()
  ### create boxplots here ###
  
  A = c(1,3,5,7,9,11)
  B = c(2,4,6,8,10,12)
  pdf('boxplots_statistics.pdf')
  boxplot(numLicks[A],numLicks[B],main='Licks',ylab ='Number of Licks', xlab ='A vs B', names = c('A','B'))
  boxplot(numEvents[A],numEvents[B],main='Bursts',ylab ='Number of Bursts', xlab ='A vs B', names = c('A','B'))
  boxplot(meanDuration[A],meanDuration[B],main='Mean Duration',ylab ='Time(indices 5/s)', xlab ='A vs B', names = c('A','B'))
  boxplot(meanTimeBetween[A],meanTimeBetween[B],main='Mean Time Between Bursts',ylab ='Time (indices 5/s)', xlab ='A vs B', names = c('A','B'))
  boxplot(meanIntensity[A],meanIntensity[B],main='Mean Intensity',ylab ='Intensity', xlab ='A vs B', names = c('A','B'))
  dev.off()
  
  ### output the csv file here ###
  csvOutput <- data.frame(matrix(ncol = 0, nrow = 6))
  csvOutput$PI<-(numLicks[A]-numLicks[B])/(numLicks[A] + numLicks[B])
  csvOutput$EventPI<- (numEvents[A]-numEvents[B])/(numEvents[A] + numEvents[B])
  csvOutput$numLicksA<-numLicks[A]
  csvOutput$numLicksB<-numLicks[B]
  csvOutput$numEventsA<-numEvents[A]
  csvOutput$numEventsB<-numEvents[B]
  csvOutput$meanDurationA<-meanDuration[A]
  csvOutput$meanDurationB<-meanDuration[B]
  csvOutput$meanTimeBetweenA<-meanTimeBetween[A]
  csvOutput$meanTimeBetweenB<-meanTimeBetween[B]
  csvOutput$meanIntensityA<-meanIntensity[A]
  csvOutput$meanIntensityB<-meanIntensity[B]
  csvOutput$medianDurationA<-medianDuration[A]
  csvOutput$medianDurationB<-medianDuration[B]
  csvOutput$medianTimeBetweenA<-medianTimeBetween[A]
  csvOutput$medianTimeBetweenB<-medianTimeBetween[B]
  csvOutput$medianIntensityA<-medianIntensity[A]
  csvOutput$medianIntensityB<-medianIntensity[B]
  
  outputName<-paste(substr(filename,1,nchar(filename)-4), '_output.csv')

  csvOutput

  write.csv(csvOutput, file = outputName)
  
  ### create statistical test results here ###
  
  
  
  t1 =t.test(csvOutput$numLicksA, csvOutput$numLicksB)
  
  t2 = t.test(csvOutput$numEventsA, csvOutput$numEventsB)
  t3 = t.test(csvOutput$meanDurationA, csvOutput$meanDurationB)
  t4 = t.test(csvOutput$meanIntensityA, csvOutput$meanIntensityB)
  
  w1 = wilcox.test(csvOutput$numLicksA, csvOutput$numLicksB)
  w2 = wilcox.test(csvOutput$numEventsA, csvOutput$numEventsB)
  w3 = wilcox.test(csvOutput$meanDurationA, csvOutput$meanDurationB)
  w4 = wilcox.test(csvOutput$meanIntensityA, csvOutput$meanIntensityB)
  
  # time in between is special because there are cases where it might be NA. 
  t5 = t.test(csvOutput$meanTimeBetweenA, csvOutput$meanTimeBetweenB)
  w5 = wilcox.test(csvOutput$meanTimeBetweenA, csvOutput$meanTimeBetweenB)
  sink("Statistical Test Results.txt")
  cat('The statistical test results are below:')
  capture.output(t1,file = "Statistical Test Results.txt", append = TRUE)
  capture.output(w1,file = "Statistical Test Results.txt", append = TRUE )
  capture.output(t2,file = "Statistical Test Results.txt", append = TRUE)
  capture.output(w2,file = "Statistical Test Results.txt", append = TRUE)
  capture.output(t3,file = "Statistical Test Results.txt", append = TRUE)
  capture.output(w3,file = "Statistical Test Results.txt", append = TRUE)
  capture.output(t4,file = "Statistical Test Results.txt", append = TRUE)
  capture.output(w4,file = "Statistical Test Results.txt", append = TRUE)
  capture.output(t5,file = "Statistical Test Results.txt", append = TRUE)
  capture.output(w5,file = "Statistical Test Results.txt", append = TRUE)
  
  sink()
}