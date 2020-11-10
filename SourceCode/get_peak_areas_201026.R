####### This code finds peaks on ICPMS data and calculates each peak area, based on Jiwoon's code
# Right now this only works for the _1.csv for some reason.


# Load appropriate libraries ---
library(tidyverse)
library(baseline)

# Load functions ---
g <- function (x) {
  fx <- f(x)
  ifelse(fx > 0, fx, 0)
}

# Set files to analyze names ----
files_to_analyze <- c("2020_10_26_Heal_AsLipidswB12_1.csv", 
                    "2020_10_26_Heal_AsLipidswB12_2.csv")
element_to_analyze <- "75As" #This code only works for 75As
smp_tags <- c("Smp_") #Only gets peak areas of peaks in the samples

# call up the pivoted ICPMS data and make into one file
ICPdata_all <- read.csv(paste0("RawDat/20201026_icapdata_secondLipidRun/Pivoted/", files_to_analyze[1])) 
#rbind(paste0("RawDat/20201026_icapdata_secondLipidRun/Pivoted/", files_to_analyze[2]))

samples <- ICPdata_all %>% select(sampleID) %>% unique() %>%
  filter(str_detect(sampleID, paste(smp_tags, collapse = '|')))

# call out one sample at a time
# loop around each sample-----
#j=1
for (j in 1:length(samples$sampleID)){
ICPdat_onesamp <- ICPdata_all %>%
  filter(sampleID == samples$sampleID[j])

# call out one element at at time
# loop around each element -----
ICPdat_onesamp_onelement <- ICPdat_onesamp %>%
  filter(element == element_to_analyze)
ICPdat_onesamp_onelement.t <- t(ICPdat_onesamp_onelement$intenstiy) %>% as.matrix()

# check out your baseline and save out the results-----
bc.als <- baseline(ICPdat_onesamp_onelement.t, lambda = 60, p = .01, maxit = 20, method='als') #l6, p0.01, maxit20
pdf(paste0("Intermediates/Baselines/",samples$sampleID[j], "_", element_to_analyze, ".pdf")) 
plot(bc.als)
text(x = 8000, y = 500, label = paste0(samples$sampleID[j], "\n", element_to_analyze))
dev.off()

# Add your baseline spectra to the same df, should save out this information somewhere, but not now-----
ICPdat_onesamp_onelement_wbl <- ICPdat_onesamp_onelement %>%
  mutate(intenstiy_baselined = c(getBaseline(bc.als)),
         intensity_corrected =  c(getCorrected(bc.als)))

# Find detected peaks (adapted from findICPpeaks(), with min intensity good for As)-----
ICP<-ICPdat_onesamp_onelement_wbl %>% select(time, intensity_corrected)
ICP_base<-ksmooth(ICP[,1],ICP[,2],kernel="normal", bandwidth = 50)
ICP_peak<-ksmooth(ICP[,1],ICP[,2],kernel="normal", bandwidth = 10)

peaks<-which(diff(sign(diff(ICP_peak[[2]])))==-2)+1
truepeaks<-peaks[which((ICP_peak[[2]][peaks]-ICP_base[[2]][peaks]) > 40)]

# Save out these plots with sample names on them!-----
pdf(paste0("Intermediates/DetectedPeaks/",samples$sampleID[j], "_", element_to_analyze, ".pdf") ,
    height = 3) 
plot(ICP[,1],ICP[,2],type='l',
     xlab="Retention Time (min)",
     bty='n',
     cex.axis=1,
     lwd=1,
     cex.lab=1)
abline(v=ICP[truepeaks,1], lty=3,col='red')
title(paste0(samples$sampleID[j], "\n", element_to_analyze))
dev.off()


peaks_df <- ICP[truepeaks,1] %>% as.data.frame() %>% rename(peak_times = '.')

ICP_reduced <- ICP %>%
  filter(time %in% peaks_df$peak_times)


plot(ICP$time,ICP$intensity_corrected, type='l')

x <- ICP$time
y <- ICP$intensity_corrected

## original linear interpolation function
f <- approxfun(x, y)

## a function zeroing out below-zero part of `f`



## local minima
x_minima <- x[which(diff(sign(diff(y))) == 2) + 1]

## break points for numerical integration
xx <- c(x[1], x_minima, x[length(x)])

## integrate peak area
intarea <- mapply(function (lwr, upr) integrate(g, lower = lwr, upper = upr)$value,
                  xx[-length(xx)], xx[-1])
Msort <- x_minima[cut(peaks_df$peak_times, x_minima)]
#sapply(peaks_56Fe_df$Time, function(a) x_minima[max(which(x_minima<a))])
match <- which(x_minima %in% Msort) + 1              #match is an index of retention time; from int, find values with indexes from match to get integrated area
area <- data.frame("int_area" = intarea[match])
peaks_df2 <- cbind(peaks_df, area)
}

