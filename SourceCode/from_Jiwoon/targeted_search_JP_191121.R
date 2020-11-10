###### this function goes through each row of m/z in our siderophore mass library and extracts EICs to align with ICP data 

###### Jiwoon edits
#190508: added "id" variable to findICPpeaks function to export peak list from individual samples (RT and intensity of "truepeaks") and save as "ICPpeaks##.csv" 
#190628: added some lines to read all "ICPpeaks##.csv" files and assign them to groups based on similar RT, RT tolerance for defining each group is changed by the variable "diff"
#also added lines to read all "ICPpeaks##.csv" files into one file with all ICP peaks information from all samples - this code generates lots of intermediate files, I didn't spend much time optimizing it so I just delete all the intermediate files when I'm done.
#190725: added a loop to process all ESI and ICP files at once (if you want to work with each individual pair of files, scroll further down after this loop)
#191121: added "background" and "scanfrq" to EICplot function to filter out ESI peaks that are lower than "background" intensity for "scanfrq" number of scans


### Function to align ICPMS and ESIMS data using cyanocobalamin (internal standard) peaks
multiMSalign <- function(mzxcms,ICPdata,timerange){
  #Plot Co in ICPMS data
  ICPComax<-ICPdata[which.max(ICPdata[,"Y_59Co"]),'Time_59Co']
  print(ICPComax)
  ### if the cyanocobalamin peak is not the max. intensity peak, inactivate the two lines above and activate the three lines below
  #ICPdata2 <- ICPdata[750:900,]              #change row numbers depending on your data to include time windows that include the true cyanocobalamine peak
  #ICPComax<-ICPdata2[which.max(ICPdata2[,"Y_59Co"]),'Time_59Co']
  #print(ICPComax)
  
  #Get retention time (in seconds) of Co peak (cobalamin) in Orbitrap data
  B12profile<-rawEIC(mzxcms,c(678.285,678.295))
  B12profile<-do.call(cbind,B12profile)
  B12maxscan<-B12profile[which.max(B12profile[,2]),1]
  OrbiComax<-mzxcms@scantime[B12maxscan]
  print(OrbiComax)
  
  Offset<-ICPComax-OrbiComax
  # Offset <- 118 # define offset in seconds here.
  mzxcms<-mzxcms
  
  #Use this to view data
  par(mfrow=c(3,1))
  plot(ICPdata[,'Time_59Co'],ICPdata[,'Y_59Co'],type='l',xlim=timerange)
  abline(v=ICPComax, lty=3,col='cyan')
  plotEIC(mzxcms,mzrange = c(400, 1500),rtrange=timerange) #can change this to jsut be the cyanocobalamin mass range if necessary
  abline(v=ICPComax, lty=3,col='cyan')
  plotEIC(mzxcms,mzrange = c(678.285, 678.295),rtrange=timerange) 
  abline(v=ICPComax, lty=3,col='cyan')
  ask<-readline("Before alignment. Press enter to view aligned data")
  
  mzxcms@scantime<-mzxcms@scantime+Offset
  
  par(mfrow=c(3,1))
  plot(ICPdata[,'Time_59Co'],ICPdata[,'Y_59Co'],type='l',xlim=timerange)
  abline(v=ICPComax, lty=3,col='cyan')
  plotEIC(mzxcms,mzrange = c(400, 1500),rtrange=timerange)
  abline(v=ICPComax, lty=3,col='cyan')
  plotEIC(mzxcms,mzrange = c(678.285, 678.295),rtrange=timerange) 
  abline(v=ICPComax, lty=3,col='cyan')
  ask<-readline("After alignment. Press enter to continue.")
  
  return(mzxcms)
}

### Find metal peaks in LC-ICPMS data
findICPpeaks<-function(id,ICPdata,element,timerange,minheight){
  ICP<-cbind(ICPdata[,paste('Time_',element,sep="")],ICPdata[,paste('Y_',element,sep="")])
  ICP_base<-ksmooth(ICP[,1],ICP[,2],kernel="normal",bandwidth=120)
  ICP_peak<-ksmooth(ICP[,1],ICP[,2],kernel="normal",bandwidth=10)
  
  peaks<-which(diff(sign(diff(ICP_peak[[2]])))==-2)+1
  truepeaks<-peaks[which((ICP_peak[[2]][peaks]-ICP_base[[2]][peaks])>minheight)]
  
  plot(ICP[,1],ICP[,2],type='l',
       xlim=timerange,
       xlab="Retention Time (min)",
       ylab=paste(element," Intensity (cps)"),
       bty='n',
       cex.axis=1,
       lwd=1,
       cex.lab=1)
  title(paste('  LC-ICPMS peaks:', element))
  abline(v=ICP[truepeaks,1], lty=3,col='red')
  
  ICPtruepeaks <- as.data.frame(ICP[truepeaks,])
  colnames(ICPtruepeaks) <- c("time", "int")
  ICPtruepeaks$index <- id
  write.csv(ICPtruepeaks, paste0("ICPpeaks", id, ".csv"))
  return(ICP[truepeaks,1])
}

### ICP Plotter that automatically scales y axis w/in time range
ICPplotter4<-function(ICPdata,element,metalpeak,timerange){
  ICP<-cbind(ICPdata[,paste('Time_',element,sep="")]/60,ICPdata[,paste('Y_',element,sep="")])
  ICP<-ICP[which(ICP[,1]>timerange[1] & ICP[,1]<timerange[2]),]
  
  plot(ICP[,1],ICP[,2],type='l',
       xlim=timerange,
       xlab="Retention Time (min)",
       ylab=paste(element," Intensity (cps)"),
       bty='n',
       cex.axis=1,
       lwd=1,
       cex.lab=1)
  abline(v=metalpeak, lty=3,col='gray48')
  title(paste('  ',element,'LC-ICPMS'),line=1,adj=0,cex=1.0)
}

### Plot EICs of masses that could be associated with your isotope pattern (e.g. apo form or differently charged form)
EICplot<-function(mzxcms,mass,title,timerange,metalpeak,background,scanfrq){
  
  mzrange<-c(-0.005,0.005)+mass
  EIC<-rawEIC(mzxcms,mzrange)
  times<-mzxcms@scantime
  df<-data.frame(cbind(times,unlist(EIC[[2]])))
  colnames(df) <- c("times", "int")
  df$int[!with(rle(df$int > background), rep(values & lengths >= scanfrq, lengths))] <- 0
  plot(df$times, df$int,
       type='l',
       lwd=1,
       xlim=timerange,
       ylab='Scaled Intensity',
       xlab='Retention Time (min)',
       col='black',
       bty='n',
       cex.axis=1,
       cex.lab=1
  )
  title(paste('  ',title,'LC-ESIMS.  EIC = ',round(mass,digits=3)),line=1,adj=0,cex=1.0)
  abline(v=metalpeak, lty=3,col='gray48')
}

### Function that makes plots
finalplotter<-function(mzxcms,ICPdata,apomass,timerange,metalpeak,siderophorename,isotopefile){
  par(mfrow=c(4,1))
  par(mar=c(4,4,2,2))
  timerange<-timerange/60
  mzxcms@scantime<-mzxcms@scantime/60
  metalpeak<-metalpeak/60
  apodiff<-52.9115
  isotopediff<-1.995328
  ICPplotter4(ICPdata,'56Fe',metalpeak,timerange)
  EICplot(mzxcms,apomass,paste(siderophorename,'    Apo mass:'),timerange,metalpeak,background,scanfrq)
  EICplot(mzxcms,apomass+apodiff,paste(siderophorename,'    56Fe mass:'),timerange,metalpeak,background,scanfrq)
  EICplot(mzxcms,apomass+apodiff-isotopediff,paste(siderophorename,'    54Fe mass:'),timerange,metalpeak,background,scanfrq)
}

#################
library(xcms)

# read siderophore library file
siderophorelib <- 'siderophore_mass list_080918.csv'
siderophore <- read.csv(siderophorelib, header = TRUE, sep = ",", strip.white=TRUE)


# Set parameters
timerange <- c(0,2000) #time range in seconds
background <- 0        #background for EIC (MS1 data)
scanfrq <- 5           #number of scans that detected peaks with masses, setting this to "1" will include all scans



##### if you want to use a loop for processing all your files, use this; if not, skip to below

mzdatafiles <- list.files(recursive=FALSE, full.names=TRUE, pattern="grad2.*.mzXML")      # change "pattern" for your specific file names
ICPdatafiles <- list.files(recursive=FALSE, full.names=TRUE, pattern="grad2_SPE.*.csv")

for (i in 1:length(mzdatafiles)){
  Orbifile <- mzdatafiles[i]
  ICPMSfile <- ICPdatafiles[i]
  savefile <- paste(gsub(".mzXML","",Orbifile),"_targeted.pdf",sep="")
  
  mzxcms <- xcmsRaw(Orbifile,profstep=0.01,profmethod="bin",profparam=list(),includeMSn=FALSE,mslevel=NULL, scanrange=NULL)
  ICPdata<-read.csv(ICPMSfile, header = TRUE, sep = ",", strip.white=TRUE)
  mzxcms<-multiMSalign(mzxcms,ICPdata,timerange)
  peaks_56Fe<-findICPpeaks(i,ICPdata,'56Fe',timerange,300)
  pdf(savefile)
  for(i in 1:nrow(siderophore)){
    finalplotter(mzxcms,ICPdata,siderophore[i,8],timerange,peaks_56Fe,siderophore[i,1],'Feisotope.csv')
  }
  dev.off()
  rm(mzxcms)
}


##### if you want to work with individual ICPMS and ESIMS data pairs, use this

Orbifile <- "./190617_grad2_SPE102.mzXML"
ICPMSfile <- "./190528_grad2_SPE102.csv"
savefile <- paste(gsub(".mzXML","",Orbifile),"_targeted.pdf",sep="")

# Read in Orbitrap file
mzxcms <- xcmsRaw(Orbifile,profstep=0.01,profmethod="bin",profparam=list(),includeMSn=FALSE,mslevel=NULL, scanrange=NULL)
#Read in ICPMS file
ICPdata<-read.csv(ICPMSfile, header = TRUE, sep = ",", strip.white=TRUE)

# First, Align ICPMS and Orbitrap data based on cyanocobalamin
mzxcms<-multiMSalign(mzxcms,ICPdata,timerange)
#To align with a constant offset instead, uncomment the following lines:
#offset <- 60 # define offset in seconds here.
#mzxcms<-mzxcms
#mzxcms@scantime<-mzxcms@scantime - offset

# Second, define LC-ICPMS peaks:
peaks_56Fe<-findICPpeaks(102,ICPdata,'56Fe',timerange,300)

# Third, Make lots of plots

pdf(savefile)
for(i in 1:nrow(siderophore)){
finalplotter(mzxcms,ICPdata,siderophore[i,8],timerange,peaks_56Fe,siderophore[i,1],'Feisotope.csv')
}
dev.off()




#############These functions are useful for playing around:
C<-12.010736
H<-1.0078250321
O<-15.9949146196
N<-14.0030740049
S<-31.9720710015
Fe<-55.93493757
e<-0.00054857990946

apomass<-832.5032



##### to create peak lists
temp <- list.files(recursive=FALSE, full.names=TRUE, pattern="ICPpeaks")
ICPpeakfiles <- do.call(rbind, lapply(temp, function(x) read.csv(x, stringsAsFactors = FALSE)))
ICPpeakfiles <- ICPpeakfiles[order(ICPpeakfiles$time),]
ICPpeakfiles$diff <- c(0, diff(ICPpeakfiles$time) > 10)
ICPpeakfiles$group <- cumsum(ICPpeakfiles$diff) + 1
ICPpeakfiles <- ICPpeakfiles[order(ICPpeakfiles$index),]
write.csv(ICPpeakfiles, "allICPpeaks.csv")



##### to merge all ICPpeaks file into one file

# read all ICP peakarea files and add sample specific information as index to each file, these are the intermediate files that need to be deleted after they're merged
temp <- list.files(recursive=FALSE, full.names=TRUE, pattern="ICPpeaks")
for (i in 1:length(temp)){
  ICPMSfile <- temp[i]
  savefile <-  paste(gsub("_peakarea.csv","",ICPMSfile), "_index", sep="")
  ICPdata <- read.csv(ICPMSfile, header = TRUE, sep = ",", strip.white=TRUE)
  ICPdata$index <- savefile
  write.csv(ICPdata, paste0(savefile, ".csv"))
}

# read all ICP peakarea files WITH INDEX and merge them into one file
temp2 <- list.files(recursive=FALSE, full.names=TRUE, pattern="index")
ICPpeakfiles <- do.call(rbind, lapply(temp2, function(x) read.csv(x, stringsAsFactors = FALSE)))
write.csv(ICPpeakfiles, "allICPpeakarea.csv")
