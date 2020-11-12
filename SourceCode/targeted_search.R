###### based on JP targeted_search_JP_191121.R

# Load libraries ---
library(tidyverse)
library(xcms)

# Set parameters
#timerange <- c(0,2000) #time range in seconds
min_height <- 60000        #minimum height to be plotted the max scan needs to be
scanfrq <- 20           #number of scans that detected targeted mass
element <- "75As"

# Name files ----
as_database_filename_lowmass <- "MetaData/AsLipidDatabase/lowmass_inclusionlist.csv"
found_icap_peaks_filename <- "Intermediates/AsLipids_ICPpeakareas.csv"

# Load files ----
as_db_lowmass <- read_csv(as_database_filename_lowmass) %>% arrange(mz)

# Try just one of the samples at the low mass samples
Orbifile <- "RawDat/20201026_QE_secondLipidRun/201028_Smp_ALOHA_elute_lowmass.mzXML"
ICPMSfile <- "RawDat/20201026_icapdata_secondLipidRun/Pivoted/2020_10_26_Heal_AsLipidswB12_1.csv" 
savefile <- paste(gsub(".mzXML","",Orbifile),"_targeted.pdf",sep="") # probably better to save this somewhere special

# Read in Orbitrap file
mzxcms_raw <- xcmsRaw(Orbifile,profstep=0.01,profmethod="bin",profparam=list(),includeMSn=FALSE,mslevel=NULL, scanrange=NULL)
#Read in ICPMS file
ICPdata<-read.csv(ICPMSfile, header = TRUE, sep = ",", strip.white=TRUE)  %>%
  filter(sampleID == "Smp_ALOHA_elute") %>%
  filter(element == "75As")

#To align with a constant offset instead, uncomment the following lines:
#offset <- 0 # define offset in seconds here.
#mzxcms<-mzxcms_raw
#mzxcms@scantime<-mzxcms@scantime - offset

# Second, define LC-ICPMS peaks - grab these from already written out file
peaks <- read_csv(found_icap_peaks_filename) %>%
  filter(sampID == "Smp_ALOHA_crude") %>%
  filter(peak_times > 120)

# Get a list of possible good lipids based on criteria at top of script
possible_lipids <- c()
for(i in 1:nrow(as_db_lowmass)){
  #i=1
  timerange <- c(500, 2500)
  cmp_mass <- as_db_lowmass$mz[i]
  cmp_name <- as_db_lowmass$Lipid_Name[i]
  detected_peak_times <- peaks$peak_times/60 
  timerange_detected_peaks <- c(min(peaks$peak_times)-60, 
                                max(peaks$peak_times)+60)/60
  timerange<-timerange/60
  
  #This checks if theres any signal between the top and bottom times
  mzxcms <- mzxcms_raw
  mzrange<-c(-0.005,0.005)+cmp_mass
  EIC<-rawEIC(mzxcms,mzrange)
  times<-mzxcms@scantime/60
  EIC_df<-data.frame(cbind(times,unlist(EIC[[2]])))
  colnames(EIC_df) <- c("times", "int")
  scans_detected_df <- EIC_df %>% filter(between(times,
                                      min(peaks$peak_times)/60-60,
                                      max(peaks$peak_times)/60+60)) %>%
    filter(int > 0) 
  scans_detected <- length(scans_detected_df$int)
  max_int <- max(EIC_df$int)
  if (scans_detected > scanfrq & max_int > min_height){
  print(paste0("plotting peaks for mass ", i))   
  possible_lipids[i] <- as_db_lowmass$Lipid_Name[i] 
  }
}
possible_lipids <- na.omit(possible_lipids)   
    

# Make plots of all the possibles
as_db_lowmass_culled <- as_db_lowmass %>%
  filter(Lipid_Name %in% possible_lipids)
pdf(savefile)
for(i in 1:length(possible_lipids)){
  timerange <- c(500, 2500)
  cmp_mass <- as_db_lowmass_culled$mz[i]
  cmp_name <- as_db_lowmass_culled$Lipid_Name[i]
  detected_peak_times <- peaks$peak_times/60 
  timerange_detected_peaks <- c(min(peaks$peak_times)-60, 
                                max(peaks$peak_times)+60)/60
  timerange<-timerange/60
  
  #This checks if theres any signal between the top and bottom times
  mzxcms <- mzxcms_raw
  mzrange<-c(-0.005,0.005)+cmp_mass
  EIC<-rawEIC(mzxcms,mzrange)
  times<-mzxcms@scantime/60
  EIC_df<-data.frame(cbind(times,unlist(EIC[[2]])))
  colnames(EIC_df) <- c("times", "int")
  scans_detected_df <- EIC_df %>% filter(between(times,
                                                 min(peaks$peak_times)/60-60,
                                                 max(peaks$peak_times)/60+60)) %>%
    filter(int > 0) 
  scans_detected <- length(scans_detected_df$int)
  max_int <- max(EIC_df$int)
  
  # This is where the plots are made
  par(mfrow=c(2,1))
  
  #This makes the ICP MS plot with the detected peaks
  plot(ICPdata$time/60,ICPdata$intenstiy,type='l',
       xlim=timerange,
       xlab="Retention Time (min)",
       ylab=paste(element," Intensity (cps)"),
       bty='n',
       cex.axis=1,
       lwd=1,
       cex.lab=1)
  abline(v=detected_peak_times, lty=3,col='gray48')
  title(paste('  ',element,'LC-ICPMS'),line=1,adj=0,cex=1.0)

  #This makes the ESI plot
  plot(EIC_df$times, EIC_df$int,
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
  title(paste('  ','LC-ESIMS.  EIC = ',round(cmp_mass, digits=3), cmp_name, "\n Max intensity = ", round(max_int, digits = -1), "Scans = ", scans_detected))
  abline(v=detected_peak_times, lty=3,col='gray48')
    }
dev.off()


