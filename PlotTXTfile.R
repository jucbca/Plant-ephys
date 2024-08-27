### Data Management of electrode array

library (ggplot2)
library (tidyr)
library(dplyr)
library(signal)
# Functions
PlotTrace <- function(trace,from=0,to=max(na.omit(trace$Time)) ){
  ## Theme...
  My_Theme =  theme(
    axis.text.x = element_text(color = "black", size = 20, angle = 0, hjust = 0.5, vjust = 0, face = "plain"),
    axis.text.y = element_text(color = "black", size = 20, angle = 0, hjust = 0 , vjust = 0.5, face = "plain"),
    axis.title.x = element_text(color = "black", size = 30, angle = 0, hjust = 0.5, vjust = 0, face = "plain"), #element_blank()
    axis.title.y = element_text(color = "black", size = 30, angle = 0, hjust = 0, vjust = 0.5, face = "plain"),
    plot.title = element_text(color = "black", size = 30, angle = 0, hjust = 0.5, vjust = 0, face = "bold"), # element_blank(),#
    legend.text = element_text(size = 20),
    legend.title = element_text(hjust = 0.1, size = 20),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(), #element_line()
    panel.grid.minor.x = element_blank()) 
  ##
  
  tracePlot <- trace %>% dplyr::filter(`Time` > from & `Time` < to ) %>% 
    ggplot(aes(`Time`, `V`, color = `Electrode`)) +
    geom_line(size = 0.5) +
    guides(color = guide_legend(override.aes = list(size = 3) )) +
    geom_vline(xintercept = stimTimes[ which(stimTimes>from & stimTimes<to) ], size = 0.1 ) + #Make vertical lines the data point immediately after the stimulus pedal was hit. And that are inside the plotting time range.
    ylab ("V") + 
    xlab ("Time (s)") + 
    labs(color = "Electrode") +
    ggtitle(paste(savename,metaData)) +
    My_Theme
  print(tracePlot)
  
  
  

  
}
NormalizeBaseline <- function(trace, normTime = 0 ){
  stimIndex = which(trace$Time == normTime)
  for(i in c(3:ncol(trace)) ){
    trace[,i] = trace[,i]-trace[stimIndex,i]
  } # Normalize the trace w/ the baseline. make all start at zero.
  trace = gather(trace, "Electrode", "V", 3:ncol(trace)) # change structure for ggplot management
  return(trace)
} ## INPUT trace spreaded i.e. each electrode as column. Won't work if data is not at 10Hz samplerate. Output: normalized, gathered trace for plotting.
CalculateAmplitude <- function(trace, stimTime){
  # calculate baseline (1s average)
  stimIndex = which(trace$Time == stimTime)
  amplitudes = c(savename)
  for(i in (1:max(electrodeNums)) ){
    print(i)
    # calculate baseline (1s average)
    BL = mean(  trace[(stimIndex-20):(stimIndex-10),which(names(trace)==i)]  )
    # calculate ss voltage (last 5s)
    ssV = mean(  trace[(stimIndex+50):(stimIndex+90),which(names(trace)==i)]   )
    # substract
    amplitudes = cbind(amplitudes, ssV-BL)
  }
  amplitudes <- as.data.frame(amplitudes) 
  names(amplitudes)=c("plant",(1:max(electrodeNums)))
  return(amplitudes)
}
CheckSampleRate <- function(trace){
  sampleRate = round(60/trace$Time[60],1) # 60 samples in x seconds that ran in the first 60 samples. (Hz)
  if(round(sampleRate) <= 10){
    print("ERROR! sample rate is not 10Hz. Quantification won't work")
  } else {
    print("Good! The sample rate is: ")
    print(sampleRate)
  }
}
filterTrace <- function(trace){
  # Change NA for the value on the previous row
  for( e in pickElectrodes){
    Electrode = e
    NArow = which(is.na(trace[,which(names(trace)==Electrode)]))
    trace[NArow,] = trace[(E=NArow-1),] 
    # smooth with Savitzky-Golay filter. 
    # n is the amount of point befor/after the nth point to be averaged
    # p is the filter order. Dont know what that is.
    trace[,which(names(trace)==Electrode)] = 
      sgolayfilt(na.omit(trace[,which(names(trace)==Electrode)]) , p = 1, n = 9) 
  }
  return(trace)
}

### INPUT ".txt" file. the columns of the electrodes should have the number of the electrode.
### OUTPUT ".csv" file. with 8 columns: Seconds, leaf,V,genotype,Stim,Date,DayTime,id.
list.files()
filename = list.files()[grep(".txt",list.files())] 
savename = substr( filename, 0, nchar(filename)-4 )
savename = strsplit(savename,split = "_")[[1]][1]
trace = read.delim(  filename  ) # import .txt file in the home folder
# trace = read.csv("venus2.csv"  )
CheckSampleRate(trace)

# Save info from E1 column name
metaData = names(trace)[grep("E1", names(trace))]
metaData = strsplit(metaData,"_")[[1]][3]


electrodeNums = NULL## Save the number of the electrode
for(i in  names(trace)[3:ncol(trace)] ){
  electrodeNums = append(electrodeNums, as.numeric( substr(i,2, 2) ))
}
names(trace)[3:ncol(trace)] = electrodeNums # Make the columns' names numeric variables


# Select just few electrodes.
pickElectrodes = c(1,2)
trace = trace[,c(1,2,which(names(trace)%in%pickElectrodes))]

# Filter trace
stimTimes = trace$Time[which(  is.na(trace$Time)  )-1]
# trace = filterTrace(trace)
# cut chunk of trace
# trace = trace[-c(which(trace$Time>400)),]


## SAVE as .csv
write.csv(trace, paste(savename,".csv",sep=""), row.names = FALSE)




#
##
###
#### NORMALIZE TRACE
gTrace = NormalizeBaseline(trace, normTime = 1)

PlotTrace(gTrace)



#
##
###
#### PLOT 
## Pick what channels to plot
# Zoom into a time window? set here limits

PlotTrace(gTrace, from=   97 , to = 101) #dplyr::filter(trace, Electrode == 1 | Electrode == 2 | Electrode == 3 | Electrode == 4) )
ggsave(paste(savename,"train.pdf", sep = ""))

PlotTrace(gTrace, from=100) #dplyr::filter(trace, Electrode == 1 | Electrode == 2 | Electrode == 3 | Electrode == 4) )
ggsave(paste(savename,"-6.pdf", sep = ""))

#
##
###
#### CALCULATE AMPLITUDES
amplitudes =  CalculateAmplitude(trace, stimTime = 25)
AllAmplitudes = rbind(AllAmplitudes, amplitudes )
write.csv(AllAmplitudes,"Amplitudes.csv", row.names = FALSE)
#AllAmplitudes = NULL


### Plot all APs
for (i in stimTimes) {
  PlotTrace(gTrace, from=   i , to = i+3)
  ggsave(paste(savename,"-AP.",i,".pdf", sep = ""))
}



# fix trace
from = 0
to = 15
e = 2
fixed = cbind( 
  dplyr::filter(gTrace, Electrode == e & Time >from & Time < to )%>% select(V),
  dplyr::filter(gTrace, Electrode == e & Time >from & Time < to )%>% select(Time)
)
plot(fixed$Time,fixed$V)

gTrace$V[min(which(gTrace$Time==from & gTrace$Electrode==e)):min(which(gTrace$Time==to& gTrace$Electrode==e))] = gTrace$V[min(which(gTrace$Time==(from-1)& gTrace$Electrode==e))]

write.csv(gTrace, paste("fixed_",savename,".csv",sep = ""), row.names = FALSE)
