####### This code finds peaks on ICPMS data and calculates each peak area


### Function to find peaks in LC-ICPMS data
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

###### calculate peak areas

## import ICPMS data and separate columns with Fe data
ICPdata<-read.csv("2019_05_29_Bundy_spe_Sample_4.csv", header = TRUE, sep = ",", strip.white=TRUE)
ICPdataFe<-read.csv("2019_05_29_Bundy_spe_Sample_4.csv", header = TRUE, sep = ",", strip.white=TRUE)[ , c(12, 30)]
ICPdataFe.t <- t(ICPdataFe$Y_56Fe)  #(need to transpose matrix to use the baseline function)

## using 'als' method as in Boiteau et al. 2016
## smaller "lambda" will make the baseline follow outline of the spectra more closely
## higher "p" will make the baseline closer to outline of the spectra
library(baseline)
bc.als <- baseline(ICPdataFe.t, lambda = 5, p = 0.01, maxit = 20, method='als')
plot(bc.als)

ICPdataFe$Y_56Fe_baseline <- c(getBaseline(bc.als))    #give baseline count
ICPdataFe$Y_56Fe_corr <- c(getCorrected(bc.als))    #give corrected count = original count - baseline count
## if you want to export your baseline and corrected counts, uncomment the line below
#write.csv(ICPdataFe, '2019_05_29_Bundy_spe_Sample_4_baseline.csv')  

peaks_56Fe <- findICPpeaks(ICPdataFe,'56Fe',timerange = c(300, 2400),300)        #identify peaks
peaks_56Fe_df <-as.data.frame(peaks_56Fe)
colnames(peaks_56Fe_df) <- c("Time")
##subset ICPdataFe to only get rows with retention times where peaks were detected
peaks_56Fe_reduced <- ICPdataFe[ICPdataFe$Time_56Fe %in% peaks_56Fe_df$Time,]

plot(ICPdataFe$Time_56Fe,ICPdataFe$Y_56Fe_corr, type='l')

x <- ICPdataFe$Time_56Fe
y <- ICPdataFe$Y_56Fe_corr

## original linear interpolation function
f <- approxfun(x, y)

## a function zeroing out below-zero part of `f`
g <- function (x) {
  fx <- f(x)
  ifelse(fx > 0, fx, 0)
}

## local minima
x_minima <- x[which(diff(sign(diff(y))) == 2) + 1]

## break points for numerical integration
xx <- c(x[1], x_minima, x[length(x)])

## integrate peak area
intarea <- mapply(function (lwr, upr) integrate(g, lower = lwr, upper = upr)$value,
                  xx[-length(xx)], xx[-1])
Msort <- x_minima[cut(peaks_56Fe_df$Time, x_minima)]
#sapply(peaks_56Fe_df$Time, function(a) x_minima[max(which(x_minima<a))])
match <- which(x_minima %in% Msort) + 1              #match is an index of retention time; from int, find values with indexes from match to get integrated area
area <- data.frame("int_area" = intarea[match])
peaks_56Fe_reduced2 <- cbind(peaks_56Fe_reduced, area)

## export result
write.csv(peaks_56Fe_reduced2,'2019_05_29_Bundy_spe_Sample_4_peakarea.csv')
