# Load appropriate libraries ---
library(baseline)
library(signal)
library(tidyverse)
detach("package:dplyr")
library(dplyr)

# Load functions ---
g <- function (x) {
  fx <- f(x)
  ifelse(fx > 0, fx, 0)
}

# Set files to analyze names ----
files_to_analyze <- c("2020_10_26_Heal_AsLipidswB12_1.csv", 
                      "2020_10_26_Heal_AsLipidswB12_2.csv", 
                    "2020_10_26_Heal_AsLipidswB12_3.csv")
element_to_analyze <- "75As" #This code only works for 75As
smp_tags <- c("_crude", "elute", "KM") #Only gets peak areas of crude samples
blank_tags <- c("blank", "prefilter")
co_rf_dat <- read_csv("Intermediates/CoRF.csv")
co_to_as <- 0.75 #this is the observed ratio in Co / As signal at > 30% B
meta_dat <- read_csv("MetaData/SampleMetaDat.csv") %>%
  select(`Sample ID`, sampleID, lipid_fraction_reconst_vol_mL, L_extracted)


# call up the pivoted ICPMS data and make into one file
ICPdata_all <- read_csv(paste0("RawDat/20201026_icapdata_secondLipidRun/Pivoted/", files_to_analyze[1]))  %>%
  rbind(read_csv(paste0("RawDat/20201026_icapdata_secondLipidRun/Pivoted/", files_to_analyze[2]))) %>%
  rbind(read_csv(paste0("RawDat/20201026_icapdata_secondLipidRun/Pivoted/", files_to_analyze[3])))

samples <- ICPdata_all %>% select(sampleID) %>% unique() %>%
  filter(str_detect(sampleID, paste(smp_tags, collapse = '|')))%>%
  filter(!str_detect(sampleID, paste(blank_tags, collapse = '|')))

peaks_df_all <- list()

# call out one sample at a time
# loop around each sample-----
for (j in 1:length(samples$sampleID)){
ICPdat_onesamp <- ICPdata_all %>%
  filter(sampleID == samples$sampleID[j])

ICPdat_onesamp_onelement <- ICPdat_onesamp %>%
  filter(element == element_to_analyze) %>%
  filter(time < 2500)
ICPdat_onesamp_onelement.t <- t(ICPdat_onesamp_onelement$intenstiy) %>% as.matrix()

# check out your baseline and save out the results-----
bc.als <- baseline(ICPdat_onesamp_onelement.t, 
                   lambda = 6, p = .01, maxit = 20, method='als') #l6, p0.01, maxit20
#plot(bc.als, xlim=c(0, 6500))
#text(x = 6000, y = 500, label = paste0(samples$sampleID[j], "\n", element_to_analyze))

# Add your baseline spectra to the same df, should save out this information somewhere, but not now-----
ICPdat_onesamp_onelement_wbl <- ICPdat_onesamp_onelement %>%
  mutate(intenstiy_baselined = c(getBaseline(bc.als)),
         intensity_corrected =  c(getCorrected(bc.als))) %>%
  mutate(intensity_smoothed = sgolayfilt(intenstiy, p = 15),
         intensity_corrected_smoothed = sgolayfilt(intensity_corrected, p = 15))
#write_csv(ICPdat_onesamp_onelement_wbl, paste0("Intermediates/Baseline_CSVs/", samples$sampleID[j], "_ICPwithbaseline.csv"))

# Find detected peaks (adapted from findICPpeaks(), with min intensity good for As)-----
ICP<-ICPdat_onesamp_onelement_wbl %>% 
  select(time, intensity_corrected_smoothed) %>% 
  as.data.frame()
ICP_base<-ksmooth(ICP$time,ICP$intensity_corrected_smoothed,
                  kernel="normal", bandwidth = 100)
ICP_peak<-ksmooth(ICP$time,ICP$intensity_corrected_smoothed,
                  kernel="normal", bandwidth = 10)

peaks<-which(diff(sign(diff(ICP_peak[[2]])))==-2)+1
truepeaks<-peaks[which((ICP_peak[[2]][peaks]-ICP_base[[2]][peaks]) > 40)]

# Save out these plots with sample names on them!-----
pdf(paste0("Intermediates/DetectedPeaks/",samples$sampleID[j], "_", 
           element_to_analyze, ".pdf") ,
    height = 9) 
par(mfrow=c(2,1))
plot(ICPdat_onesamp_onelement_wbl$time/60,
     ICPdat_onesamp_onelement_wbl$intenstiy,type='l',
     xlab="Retention Time (min)",
     ylab = "Intensity",
     bty='n',
     cex.axis=1,
     lwd=1,
     cex.lab=1,
     xlim=c(5,40))
title(paste0(samples$sampleID[j], "\n", element_to_analyze, " raw ICP signal"))

plot(ICPdat_onesamp_onelement_wbl$time/60,
     ICPdat_onesamp_onelement_wbl$intensity_smoothed,type='l',
     xlab="Retention Time (min)",
     ylab = "Intensity",
     bty='n',
     cex.axis=1,
     lwd=1,
     cex.lab=1,
     xlim=c(5,40))
points(ICPdat_onesamp_onelement_wbl$time/60,
       ICPdat_onesamp_onelement_wbl$intenstiy_baselined,
       type='l',
       col = 'red')
abline(v=ICP[truepeaks,1]/60, lty=3,col='grey')
title(paste0(samples$sampleID[j], "\n", element_to_analyze, " baselined and detected peaks"))
dev.off()

peaks_df <- ICP[truepeaks,1] %>% as.data.frame() %>% rename(peak_times = '.')

ICP_reduced <- ICP %>%
  filter(time %in% peaks_df$peak_times)
x <- ICP$time
y <- ICP$intensity_corrected

## original linear interpolation function
f <- approxfun(x, y)

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
peaks_df_all[[j]] <- cbind(peaks_df, area) %>% mutate(sampID = samples$sampleID[j])
}

peaks_df_all <- do.call(rbind, peaks_df_all) 

peaks_df_all_concentrations <- peaks_df_all%>%
  mutate(nMolAs_in_vial = int_area/co_to_as*(co_rf_dat$concentration[1]/co_rf_dat$area[1])) %>%
  rename(sampleID = sampID) %>%
  left_join(meta_dat, by = "sampleID") %>%
  mutate(pMolAs_enviro = nMolAs_in_vial/10^3*lipid_fraction_reconst_vol_mL/L_extracted*10^3) %>%
  select(-lipid_fraction_reconst_vol_mL, -L_extracted)


write_csv(peaks_df_all_concentrations, "Intermediates/AsLipids_ICPpeakareas_concen.csv")

# This is starting to explore if we can assign the same compounds to each of these, 
peaks_df_all_2 <- peaks_df_all_concentrations %>%
  filter(sampleID %in% meta_dat$sampleID) %>%
  as.data.frame() %>% 
  mutate(present = 1) 
g <- ggplot(data = peaks_df_all_2, aes(x = peak_times, y = sampleID, color = log(pMolAs_enviro), size = log(pMolAs_enviro))) +
              geom_point()
g
