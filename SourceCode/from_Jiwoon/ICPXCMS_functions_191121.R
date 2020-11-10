

###########Function to find isotope matches
isotopehunter8<-function(mzxcms,isotopefile,scantime){
  
  starttime<-Sys.time()
  
  if(scantime[1]=="all"){
    startscan<-1
    endscan<-length(mzxcms@scantime)
  }else{
    startscan<-which(mzxcms@scantime>scantime[1])[1]
    endscan<-tail(which(mzxcms@scantime<scantime[2]),1)
  }
  
  readpattern<-read.csv(isotopefile, header = TRUE, sep = ",")
  
  pattern<-readpattern[(readpattern[,6]=='Y'),]
  uppermass<-pattern[-1,2]+pattern[-1,4]-pattern[1,2]
  lowermass<-pattern[-1,2]-pattern[-1,4]-pattern[1,2]
  upperratio<-pattern[-1,3]*pattern[-1,5]/pattern[1,3]
  lowerratio<-pattern[-1,3]/pattern[-1,5]/pattern[1,3]
  nisotope<-length(uppermass)
  
  columnnames<-c(sapply(pattern[,1],FUN=function(x) c(paste(x,'mass'),paste(x,'intensity'))))
  
  results<-matrix(0,ncol=(4+nisotope*2),nrow=1E6)
  r<-1
  colnames(results) = c('scan','time', columnnames)
  
  
  mzint<-20
  mzbuffer<-uppermass[length(uppermass)]

  for(i in startscan:endscan) {
    tscan<-getScan(mzxcms,i)
    tscan<-tscan[which(tscan[,2]>0),]
    
    print(paste(i,'of',endscan))
    
    first<-tscan[1,1]
    last<-first+mzint
    
    while(last<tscan[nrow(tscan),1]){
      
      scan<-tscan[which(tscan[,1]>first & tscan[,1]<last),]
      
      if(length(scan)>2){
        final<-length(which(scan[,1]<(last-mzbuffer)))
        for(j in 1:final){
          x<-scan[j,1]
          y<-scan[j,2]
          isotopes<-matrix(NA,ncol=nisotope)
          k<-1
          while(k<nisotope+1){
            isotopes[k]<-(which(
              scan[,1]>(x+lowermass[k]) & 
                scan[,1]<(x+uppermass[k]) & 
                scan[,2]>(y*lowerratio[k]) & 
                scan[,2]<(y*upperratio[k])
            )[1])
            if(is.na(isotopes[k])){k<-(nisotope+2)
            } else {k<-k+1
            }
          }
          if(k==(nisotope+1)){
            results[r,]<-c(i,mzxcms@scantime[i],x,y,c(t(scan[isotopes,])))
            r<-r+1
          }
          
        }
      }
      first<-last-mzbuffer
      last<-last+mzint
    }
    scan<-tscan[which(tscan[,1]>first & tscan[,1]<last),]
    if(length(scan)>2){
      for(j in 1:nrow(scan)){
        x<-scan[j,1]
        y<-scan[j,2]
        isotopes<-matrix(NA,ncol=nisotope)
        k<-1
        while(k<nisotope+1){
          isotopes[k]<-(which(
            scan[,1]>(x+lowermass[k]) & 
              scan[,1]<(x+uppermass[k]) & 
              scan[,2]>(y*lowerratio[k]) & 
              scan[,2]<(y*upperratio[k])
          )[1])
          if(is.na(isotopes[k])){k<-(nisotope+2)
          } else {k<-k+1
          }
        }
        if(k==(nisotope+1)){
          results[r,]<-c(i,mzxcms@scantime[i],x,y,c(t(scan[isotopes,])))
          r<-r+1
        }
      }
    }
  }
  if('O' %in% readpattern[,6]){
    optionalpattern<-readpattern[(readpattern[,6]=='O')]   
    
    for(i in 1:nrow(results)){
      
    }
    
  }
  
  
  print(Sys.time()-starttime)
  
  return(results[1:(r-1),])
}

########Simplify list:
Removenoise4<-function(result,background,twidth){
  
  #remove elements below x counts (x=background)
  results2<-result[(result[,4]>background),]
  
  #round data for binning
  results3<-round(results2, digits=2)
  
  #remove elements that don't appear twice or more in x chromatographic seconds (x=twidth) and which don't have a modest peak if there are more than 10 points.
  masslist<-unique(results3[which(duplicated(results3[,3])),3])
  masslist2<-c()

  #remove elements that appear in x or more 30 second intervals (x = noise).
  for(i in 1:length(masslist)){
    scantime<-0
    intense<-NULL
    intense2<-NULL
    times<-results3[which(results3[,3]==masslist[i]),2]
    timediff<-times[2:length(times)]-times[1:length(times)-1]
    
if(min(timediff)<twidth){
    while(scantime<max(result[,2])){
      resultint<-results3[which(results3[,3]==masslist[i] & results3[,2]>scantime & results3[,2]<(scantime+30)),4]
      resultint2<-results3[which(results3[,3]==masslist[i] & results3[,2]>scantime & results3[,2]<(scantime+30)),6]
      if(length(resultint)>0){
        intense<-c(intense,max(resultint))
        intense2<-c(intense2,max(resultint2))
      }
      scantime<-(scantime+30)
    }
    intense<-sort(intense)
    intense2<sort(intense2)

    if(length(intense)<7.5){masslist2<-c(masslist2,masslist[i])}
    if(length(intense)>7.5){
      quartile<-as.integer(length(intense)*3/4)
      threshold1<-mean(intense[1:quartile]+3*sd(intense[1:quartile]))
      threshold2<-mean(intense2[1:quartile]+3*sd(intense[1:quartile]))
      if(max(intense)>threshold1 & max(intense2)>threshold2){masslist2<-c(masslist2,masslist[i])}
      }
    }
}
  
  
  ##Now return only the result with the maximum intensity for each element:
  
  result4<-matrix(data=0,nrow=length(masslist2),ncol=ncol(result))
  colnames(result4) = colnames(result)
  
  for(i in 1:length(masslist2)){
    resultint<-results2[(results3[,3]==masslist2[i]),]
    result4[i,]<-resultint[which.max(resultint[,4]),]
  }
  
  return(result4)
}

### Automatically find peaks in LC-ICPMS data
findICPpeaks<-function(ICPdata,element,timerange,minheight){
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
  
  
  return(ICP[truepeaks,1])
}
  
##########Plotting functions
#Isotope plotter for final report. Makes one plot with EIC's of every isotope found, scaled to the 'same' intensity.
isotopeplotter3<-function(mzxcms,isotopefile,lowmass,timerange,maxscan,metalpeak){
  pattern<-read.csv(isotopefile, header = TRUE, sep = ",")
  
  plotrange<-pattern[,2]-pattern[1,2]+lowmass
  isotoperatio<-pattern[,3]/pattern[1,3]
  namerange<-as.vector(pattern[,1])
  
  colors<-c('blue4','darkorange2','burlywood4','black','red')
  usedcolors<-colors[1:ncol(pattern)]
  
  maxy<-0
  times<-mzxcms@scantime
  
  for(i in 1:nrow(pattern)){
    mzrange<-c(-0.002,0.002)+plotrange[i]
    EIC<-rawEIC(mzxcms,mzrange)
    EIC1<-unlist(EIC[2])/isotoperatio[i]
    newmax<-max(EIC1[which(times>timerange[1]&times<timerange[2])])
    maxy<-max(c(maxy,newmax))
  }
  
  
  mzrange<-c(-0.002,0.002)+lowmass
  EIC<-rawEIC(mzxcms,mzrange)
  plot(times,EIC[[2]],
       type='l',
       lwd=1,
       xlim=timerange,
       ylim=c(0,maxy*1.1),
       ylab='Scaled Intensity',
       xlab='Retention Time (min)',
       col=usedcolors[1],
       bty='n',
       cex.axis=1,
       cex.lab=1
  )
  
  for(i in 2:ncol(pattern)){
    mzrange<-c(-0.002,0.002)+plotrange[i]
    EIC<-rawEIC(mzxcms,mzrange)
    lines(times,EIC[[2]]/isotoperatio[i],col=usedcolors[i],lwd=1)
  }
  
  abline(v=maxscan, lty=2,col='red')
  abline(v=metalpeak, lty=3,col='gray48')
  title(paste('  Isotope peaks LC-ESIMS.  EIC = ',toString(round(plotrange,digits=3))))
  legend('topright',bty="n",namerange,lwd=2,col=usedcolors,cex=0.7,horiz=TRUE)
  
}

#Plot other masses that could be associated with your isotope pattern (e.g. apo form or differently charged form)
EICplot<-function(mzxcms,mass,title,timerange,metalpeak,maxpeak){

  mzrange<-c(-0.002,0.002)+mass
  EIC<-rawEIC(mzxcms,mzrange)
  times<-mzxcms@scantime
  plot(times,EIC[[2]],
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
  abline(v=maxpeak, lty=2,col='red')
}

#ICPMS chromatogram plotter
#ICP Plotter that automatically scales y axis w/in time range
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
  mtext(metalpeak, side=3)
  title(paste('  ',element,'LC-ICPMS'),line=1,adj=0,cex=1.0)
}

#Final plot function for generating report for any pattern. Includes MS spectra as last panel. X axis in minutes.
finalplotter_MS2<-function(mzxcms,lowmass,results,isotopefile,timerange,ICPdata,element,metalpeak,i,apodiff){
  par(mfrow=c(4,1))
  par(mar=c(4,4,2,2))
  timerange<-timerange/60
  mzxcms@scantime<-mzxcms@scantime/60
  metalpeak<-metalpeak/60
  results[,2]<-results[,2]/60
  
  ICPplotter4(ICPdata,element,metalpeak,timerange)
  isotopeplotter3(mzxcms,isotopefile,lowmass,timerange,results[i,2],metalpeak)
  EICplot(mzxcms,lowmass-apodiff,paste('    Apo mass: m/z =',lowmass-apodiff),timerange,metalpeak,results[i,2])
  MSplotter2(mzxcms,results,i,isotopefile)
}

###########Plot MS spectrum
MSplotter2<-function(mzxcms,results,i,isotopefile){
  
  #determine correct isotope pattern
  readpattern<-read.csv(isotopefile, header = TRUE, sep = ",")
  pattern<-readpattern[(readpattern[,6]=='Y'),]
  theorymass<-pattern[,2]-pattern[1,2]+results[i,3]
  theoryratio<-pattern[,3]/pattern[1,3]*results[i,4]
  
  isotopes<-matrix(results[i,-1:-2],byrow=TRUE,ncol=2)
  massrange<-c(mean(isotopes[,1])-10,mean(isotopes[,1])+10)
  scanrange<-getScan(mzxcms,results[i,1],mzrange=massrange)
  #scanrange<-scan[which(scan[,1]>mzrange[1] & scan[,1]<mzrange[2]),]
  plot(scanrange,type='h',ylim=c(0,max(scanrange[,2])))
  lines(theorymass,theoryratio,type='h',col='cornsilk3',lwd=3)
  lines(isotopes,type='h',col='red')
  title(paste('MS: ',toString(round(results[i,2],digits=2)),' min'))
  
}

####### Generate a final report for any metal
ICPXCMSreport2<-function(savefile,mzxcms,results,isotopefile,timerange,ICPdata,element,metalpeak){
  
  masslist<-results[,3]
  
  pdf(paste(savefile,'_',element,'.pdf'))
  for(i in 1:length(masslist)){
    finalplotter_MS(mzxcms,masslist[i],results,isotopefile,timerange,ICPdata,element,metalpeak,i)
  }
  dev.off()
  return()
}

####### Generate a final report for any metal plus an 'Apo' plot
ICPXCMSreport3<-function(savefile,mzxcms,results,isotopefile,timerange,ICPdata,element,metalpeak,apodiff){
  
  masslist<-results[,3]
  
  pdf(paste(savefile,'_',element,'.pdf'))
  for(i in 1:length(masslist)){
    finalplotter_MS2(mzxcms,masslist[i],results,isotopefile,timerange,ICPdata,element,metalpeak,i,apodiff)
  }
  dev.off()
  return()
}
############## Get MS2 spectrum
ms2spectra<-function(xs,mass,time,labelmin){
  masslow<-mass-0.2
  masshigh<-mass+0.2
  timelow<-time-15
  timehigh<-time+15
  scanlist<-which(xs@msnPrecursorMz>masslow & xs@msnPrecursorMz<masshigh & xs@scantime>timelow & xs@scantime<timehigh)
  scanlist2<-scanlist[which.max(xs@msnPrecursorIntensity[scanlist])]
  #   print(scanlist)
  #   print(xs@msnPrecursorIntensity[scanlist])
  #   print(scanlist2)
  #   print(xs@msnPrecursorMz[scanlist2[1]])
  #   print(xs@scantime[scanlist2[1]])
  #   print(xs@msnLevel[scanlist2[1]])
  if(scanlist2!=0){
    scanrange<-getScan(xs,scanlist2[1])
    maxscans<-scanrange[which(scanrange[,2]>labelmin),]
    #    print(maxscans)
    plot(scanrange,type='h',xlab='m/z',ylab='Intensity',ylim=c(0,1.1*max(scanrange[,2])),xlim=c(150,620))
    text(maxscans,paste(round(maxscans[,1],digits=3)),cex=1,pos=3)
    title(paste('  ms2:',round(xs@msnPrecursorMz[scanlist2[1]],digits=2)),line=-1,adj=0,cex=1.0)
  }else{print('No MS2 spectra available')}
}
  
 