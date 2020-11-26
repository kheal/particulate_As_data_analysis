# Load libraries ---
library(tidyverse)
library(xcms)
library(CluMSID)
library(CluMSIDdata)
library(Rdisop)
library(fuzzyjoin)
options(digits=10)

# Name files ----
sample_matcher_filename <- "MetaData/QE_ICAP_samplematcher.csv"
location_of_QEfiles <- "RawDat/20201026_QE_secondLipidRun"
MS2_spectra_filename <- "MetaData/AsLipidDatabase/MS2_Fragments.csv"
MS1s_library_filename_1 <- "MetaData/AsLipidDatabase/highmass_inclusionlist.csv"
MS1s_library_filename_2 <- "MetaData/AsLipidDatabase/lowmass_inclusionlist.csv"
found_icap_peaks_filename <- "Intermediates/AsLipids_ICPpeakareas.csv"


#Load file of the samples to analyze, load in detected peaks
sample_matcher <- read_csv(sample_matcher_filename) %>%
  filter(str_detect(ICP_sampID, "crude|elute|KM")) %>%
  filter(!str_detect(ICP_sampID, "lank|filter")) 
MS2_library <- read_csv(MS2_spectra_filename) %>%
  mutate(mz = sapply(Formula, 
                     function(x)getMolecule(x)$exactmass))
MS1_library <- read_csv(MS1s_library_filename_1) %>%
  bind_rows(read_csv(MS1s_library_filename_2))
peaks_all <- read_csv(found_icap_peaks_filename)


#Loop through sample_matcher$QE_filename[j] on j ----------
for (j in 1:length(sample_matcher$QE_filename)){
#j=1
QE_file <- paste0(location_of_QEfiles, "/", sample_matcher$QE_filename[j])
savefile <- paste("Intermediates/MS2_matched_only/", gsub(".mzXML","",sample_matcher$QE_filename[j]),"_targeted.pdf",sep="")

# Load up the QE mzXML file
mzxcms_raw <- xcmsRaw(QE_file,profstep=0.01,profmethod="bin",profparam=list(),includeMSn=FALSE,mslevel=NULL, scanrange=NULL)

# Grab ICPMS peaks from already written out file
peaks <- peaks_all %>%
  filter(sampID == sample_matcher$ICP_sampID[j]) %>%
  filter(peak_times > 120) # don't bother looking at the early peaks
detected_peak_times <- peaks$peak_times/60 

# Grab ICP data
ICPdata <- read_csv(paste0("Intermediates/Baseline_CSVs/", 
                           sample_matcher$ICP_sampID[j], "_ICPwithbaseline.csv"))


pdf(savefile, width = 7, height = 10)
# Grab the MS2 data
savedfile <- str_replace(QE_file, ".mzXML", '_mergedMS2.rds')
my_spectra <- readRDS(savedfile)

Matched_MS2 <- list()
for (i in 1:length(MS2_library$Formula))
  {Matched_MS2_small <- findFragment(my_spectra, 
                               getMolecule(MS2_library$Formula[i])$exactmass,
                               tolerance = 3E-06)
Matched_MS2 <- Matched_MS2 %>% append(Matched_MS2_small)
}

Matched_MS2 <- Matched_MS2 %>% unique()
MS1s <- map_chr(Matched_MS2, "precursor") %>% as_tibble() %>%
  mutate(mz = round(as.numeric(value), digits = 5)) %>%
  mutate(rt = as.numeric(map_chr(Matched_MS2, "rt"))) %>% 
  mutate() %>%
  distance_left_join(MS1_library, max_dist = 0.005)%>%
  filter(is.na(LipidClass)) %>%
  mutate(mz = mz.x) %>%
  select(mz, rt) %>%
  arrange(mz) %>%
  filter(rt > min(peaks$peak_times)) %>%
  filter(rt < max(peaks$peak_times))

for (k in 1:length(MS1s$mz)){
  cmp_mass <- MS1s$mz[k]
  cmp_mass_C13 <- MS1s$mz[k]+1.00335
  cmp_mass_S34 <- MS1s$mz[k]+(33.96787-31.97207)
  MS2_time <- MS1s$rt[k]
  
  # Get the MS1 chromats ready to plot
  mzxcms <- mzxcms_raw
  mzrange<-c(-0.002,0.002)+cmp_mass
  EIC<-rawEIC(mzxcms,mzrange)
  times<-mzxcms@scantime/60
  EIC_df<-data.frame(cbind(times,unlist(EIC[[2]])))
  colnames(EIC_df) <- c("times", "int")
  
  EIC_13C<-rawEIC(mzxcms,c(-0.002,0.002)+cmp_mass_C13)
  EIC_13C_df<-data.frame(cbind(times,unlist(EIC_13C[[2]])))
  colnames(EIC_13C_df) <- c("times", "int")
  
  EIC_34S<-rawEIC(mzxcms,c(-0.002,0.002)+cmp_mass_S34)
  EIC_34S_df<-data.frame(cbind(times,unlist(EIC_13C[[2]])))
  colnames(EIC_13C_df) <- c("times", "int")
  
  # Get the MS1 data ready to plot
  MS1_extracted <- getSpec(mzxcms, rtrange = c(MS2_time-5, MS2_time+10))
  MS1_df<-MS1_extracted %>% as.data.frame() %>%
    filter(mz > cmp_mass-4) %>%
    filter(mz < cmp_mass+4)%>%
    arrange(desc(mz)) %>%
    mutate(intensity = ifelse(is.na(intensity), 0, intensity))
  MS1_1mz <- MS1_df %>% 
    filter(mz > cmp_mass-.04) %>%
    filter(mz < cmp_mass+.04) 
  MS1_observed <- weighted.mean(MS1_1mz$mz, MS1_1mz$intensity)
  
  # Get the MS2 data ready to plot
  MS2s_of_mass <- getSpectrum(my_spectra, "precursor", 
                              cmp_mass, mz.tol = 1E-03)
  if (length(MS2s_of_mass) > 1){
  MS2 <- getSpectrum(MS2s_of_mass,"rt", MS2_time, rt.tol = 1)} else{
    MS2 <-  MS2s_of_mass
  }
  MS2_to_plot <- tibble(mz = as.numeric(MS2@spectrum[,1]),
                        intensity = as.numeric(MS2@spectrum[,2]))
  
  MS2_to_plot_withAs <- MS2_to_plot %>%
    distance_left_join(MS2_library %>% #select(mz, Formula) %>% 
                         unique(), max_dist = 0.005) %>%
    filter(!is.na(Formula))
  
  # plot all the results together
  #This makes the ICP MS plot with the detected peaks
  par(mfrow=c(5,1))
  
  plot(ICPdata$time/60,ICPdata$intensity_smoothed,type='l',
       xlab="Retention Time (min)",
       ylab = "Intensity",
       bty='n',
       cex.axis=1,
       lwd=1,
       cex.lab=1,
       xlim=c(8,35))
  points(ICPdata$time/60,
         ICPdata$intenstiy_baselined,
         type='l',
         col = 'red')
  abline(v=detected_peak_times, lty=3,col='gray48')
  title(paste('LC-ICPMS'),line=1,adj=0,cex=1.0)
  
  #This makes the ESI plot for the C12
  plot(EIC_df$times, EIC_df$int,
       type='l',
       lwd=1,
       xlim=c(8,35),
       ylab='Intensity',
       xlab='Retention Time (min)',
       col='black',
       bty='n',
       cex.axis=1,
       cex.lab=1
  )
  title(paste('  ','LC-ESIMS.  EIC = ',round(cmp_mass, digits=4)), line=1,adj=0,cex=1.0)
  abline(v=detected_peak_times, lty=3,col='gray48')
  abline(v=MS2_time/60, lty=1,lwd = 2, col='green')
  
  #This makes the 13C plot
  plot(EIC_13C_df$times, EIC_13C_df$int,
       type='l',
       lwd=1,
       xlim=c(8,35),
       ylab='Intensity',
       xlab='Retention Time (min)',
       col='black',
       bty='n',
       cex.axis=1,
       cex.lab=1
  )
  title(paste('  ','LC-ESIMS of 13C.  EIC = ',round(cmp_mass_C13, digits=4)),
        line=1,adj=0,cex=1.0)
  abline(v=detected_peak_times, lty=3,col='gray48')

  #Plot the MS1
  plot(MS1_df$mz, MS1_df$intensity,
       type = 'h', 
    #   xlim=c(cmp_mass-3,cmp_mass+5),
       ylab='Intensity',
       xlab='m/z',
       col='black',
       bty='n',
       cex.axis=1,
       cex.lab=1)
  abline(v=cmp_mass, lty=1,lwd = 3, col='green')
  title(paste('Calculated exact mass:', round(MS1_observed, digits = 5)),
        line=1,adj=0,cex=1.0)

  
  #Plot the MS2
  plot(MS2_to_plot$mz, MS2_to_plot$intensity,
       type = 'h', 
       xlim=c(90,500),
       ylab='Intensity',
       xlab='m/z',
       col='black',
       bty='n',
       cex.axis=1,
       cex.lab=1)
  points(MS2_to_plot_withAs$mz.x, MS2_to_plot_withAs$intensity,
         type = 'h', col = 'green', lwd = 2) 
  title(paste('MS2 spectra for mass', round(cmp_mass, digits = 4),
              "at", round(MS2_time/60, digits = 1), "minutes"),
        line=1,adj=0,cex=1.0)
  
}
dev.off()

}





