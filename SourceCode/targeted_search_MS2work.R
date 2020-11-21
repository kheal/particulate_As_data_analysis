###### based on JP targeted_search_JP_191121.R

# Load libraries ---
library(tidyverse)
library(xcms)
library(CluMSID)
library(CluMSIDdata)

# Set parameters
#timerange <- c(0,2000) #time range in seconds
min_height <- 50000        #minimum height to be plotted the max scan needs to be; 100000 is good for low mass
scannum <- 20           #number of scans that detected targeted mass
element <- "75As"
scanfrq <- 7  #scans in a row we need to keep before counting it as an actual scan
background <- 5000 #intensity background to keep 

# Name files ----
as_database_filename_lowmass <- "MetaData/AsLipidDatabase/lowmass_inclusionlist.csv"
as_database_filename_highmass <- "MetaData/AsLipidDatabase/highmass_inclusionlist.csv"
sample_matcher_filename <- "MetaData/QE_ICAP_samplematcher.csv"
location_of_QEfiles <- "RawDat/20201026_QE_secondLipidRun"

#Function to getMs2 spectra -----

#Load file of the samples to analyze, load in detected peaks
sample_matcher <- read_csv(sample_matcher_filename) %>%
  filter(str_detect(ICP_sampID, "crude|elute|KM")) %>%
  filter(!str_detect(ICP_sampID, "lank|filter")) 

j=2
QE_file <- paste0(location_of_QEfiles, "/", sample_matcher$QE_filename[j])
print(QE_file)

#This should get us the spectrum of a certain mass
my_spectra <- extractMS2spectra(QE_file,  min_peaks = 1)
my_spectra_merged <- mergeMS2spectra(my_spectra, mz_tolerance = 3E-06, 
                                     rt_tolerance = 1, peaktable = NULL, 
                                     exclude_unmatched = FALSE)

#AsSugars
MS2test <- findFragment(my_spectra, getMolecule("C10H22O7As")$exactmass, tolerance = 3E-06)
#31 found with this mass in ALOHA elute high mass

MS2test <- findFragment(my_spectra, getMolecule("C10H23O10AsP")$exactmass, tolerance = 5E-06) #38 found with this mass in ALOHA elute high mass

MS2test <- findFragment(my_spectra, getMolecule("C2H8OAs")$exactmass, tolerance = 5E-06) 
#19 found with this mass in ALOHA crude low mass

MS2test <- findNL(my_spectra, getMolecule("C2H4As")$exactmass, tolerance = 5E-06) 
#12 found with this mass in ALOHA crude low mass
test <- findNL(my_spectra, 18.0105647) #this only works after you've run mergeMS2spectra, which fills in the neutral loss aspect

my_spectra_merged <- mergeMS2spectra(my_spectra, mz_tolerance = 3E-06,
                                     rt_tolerance = 20, peaktable = NULL, 
                                     exclude_unmatched = FALSE)




ICP_sampID <- sample_matcher$ICP_sampID[j]
if(str_detect(QE_file, "highmass")){
  as_db <- read_csv(as_database_filename_highmass) %>% arrange(mz)
} else {as_db <- read_csv(as_database_filename_lowmass) %>% arrange(mz)}
savefile <- paste("Intermediates/Targeted_search/", gsub(".mzXML","",sample_matcher$QE_filename[j]),"_targeted.pdf",sep="")

# Read in Orbitrap file
mzxcms_raw <- xcmsRaw(QE_file,profstep=0.01,profmethod="bin",profparam=list(),includeMSn=FALSE,mslevel=NULL, scanrange=NULL)

# Grab ICPMS peaks from already written out file
peaks <- peaks_all %>%
  filter(sampID == sample_matcher$ICP_sampID[j]) %>%
  filter(peak_times > 120) # don't bother looking at the early peaks

# Get a list of possible good lipids based on criteria at top of script, plot each and ask if we should keep it as a possibility
possible_lipids <- c()
for(i in 1:nrow(as_db)){
  #i=1
  timerange <- c(500, 2500)
  cmp_mass <- as_db_culled$mz[i]
  cmp_mass_C13 <- as_db_culled$mz_13C[i]
  ratioC13_C12 <- as_db_culled$ratio_12Cto13C[i]
  cmp_name <- as_db_culled$Lipid_Name[i]
  detected_peak_times <- peaks$peak_times/60 
  timerange_detected_peaks <- c(min(peaks$peak_times)-60, 
                                max(peaks$peak_times)+60)/60
  
  mzxcms <- mzxcms_raw
  mzrange<-c(-0.005,0.005)+cmp_mass
  EIC<-rawEIC(mzxcms,mzrange)
  times<-mzxcms@scantime/60
  EIC_df<-data.frame(cbind(times,unlist(EIC[[2]])))
  colnames(EIC_df) <- c("times", "int")
  EIC_13C<-rawEIC(mzxcms,c(-0.005,0.005)+cmp_mass_C13)
  EIC_13C_df<-data.frame(cbind(times,unlist(EIC_13C[[2]])))
  colnames(EIC_13C_df) <- c("times", "int")
  EIC_df$int[!with(rle(EIC_df$int > background),
                   rep(values & lengths >= scanfrq, lengths))] <- 0
  scans_detected_df <- EIC_df %>% filter(between(times,
                                                 min(peaks$peak_times)/60-60,
                                                 max(peaks$peak_times)/60+60)) %>%
    filter(int > 0) 
  scans_detected <- length(scans_detected_df$int)
  max_int <- max(EIC_df$int)
  
  if (scans_detected > scannum & max_int > min_height)
  {  # Grab ICP data
    ICPdata <- read_csv(paste0("Intermediates/Baseline_CSVs/", 
                               sample_matcher$ICP_sampID[j], "_ICPwithbaseline.csv"))
    
    
    # This is where the plots are made
    par(mfrow=c(2,1))
    
    #This makes the ICP MS plot with the detected peaks
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
    title(paste(i,'  ',element,'LC-ICPMS'),line=1,adj=0,cex=1.0)
    
    #This makes the ESI plot
    plot(EIC_df$times, EIC_df$int,
         type='l',
         lwd=1,
         xlim=c(8,35),
         ylab='Intensity (scaled for 13C)',
         xlab='Retention Time (min)',
         col='black',
         bty='n',
         cex.axis=1,
         cex.lab=1
    )
    lines(EIC_13C_df$times, EIC_13C_df$int*ratioC13_C12, 
          type='l',
          col = '#FF000088')
    title(paste('  ','LC-ESIMS.  EIC = ',round(cmp_mass, digits=3), cmp_name, "\n Max intensity = ", round(max_int, digits = -1), "Scans = ", scans_detected))
    abline(v=detected_peak_times, lty=3,col='gray48')
    ask<-readline(prompt="Enter 'y' if this is a possible match: ")
    if(ask=='y'){possible_lipids[i] <- as_db$EmpiricalFormula[i]}
    }
}
possible_lipids <- na.omit(possible_lipids)   
print(paste(length(possible_lipids), "possible hits"))

# Make plots of all the possibles
as_db_culled <- as_db %>%
  filter(EmpiricalFormula %in% possible_lipids)
if (length(possible_lipids) > 0){
  pdf(savefile)
  for(i in 1:length(possible_lipids)){
    timerange <- c(500, 2500)
    cmp_mass <- as_db_culled$mz[i]
    cmp_mass_C13 <- as_db_culled$mz_13C[i]
    ratioC13_C12 <- as_db_culled$ratio_12Cto13C[i]
    cmp_name <- as_db_culled$Lipid_Name[i]
    detected_peak_times <- peaks$peak_times/60 
    timerange_detected_peaks <- c(min(peaks$peak_times)-60, 
                                  max(peaks$peak_times)+60)/60
    
    mzxcms <- mzxcms_raw
    mzrange<-c(-0.005,0.005)+cmp_mass
    EIC<-rawEIC(mzxcms,mzrange)
    times<-mzxcms@scantime/60
    EIC_df<-data.frame(cbind(times,unlist(EIC[[2]])))
    colnames(EIC_df) <- c("times", "int")
    EIC_13C<-rawEIC(mzxcms,c(-0.005,0.005)+cmp_mass_C13)
    EIC_13C_df<-data.frame(cbind(times,unlist(EIC_13C[[2]])))
    colnames(EIC_13C_df) <- c("times", "int")
    EIC_df$int[!with(rle(EIC_df$int > background),
                     rep(values & lengths >= scanfrq, lengths))] <- 0
    scans_detected_df <- EIC_df %>% filter(between(times,
                                                   min(peaks$peak_times)/60-60,
                                                   max(peaks$peak_times)/60+60)) %>%
      filter(int > 0) 
    scans_detected <- length(scans_detected_df$int)
    max_int <- max(EIC_df$int)
    
    
    # Grab ICP data
    ICPdata <- read_csv(paste0("Intermediates/Baseline_CSVs/", 
                               sample_matcher$ICP_sampID[j], "_ICPwithbaseline.csv"))
    
    
    # This is where the plots are made
    par(mfrow=c(2,1))
    
    #This makes the ICP MS plot with the detected peaks
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
    title(paste('  ',element,'LC-ICPMS'),line=1,adj=0,cex=1.0)
    
    #This makes the ESI plot
    plot(EIC_df$times, EIC_df$int,
         type='l',
         lwd=1,
         xlim=c(8,35),
         ylab='Intensity (scaled for 13C)',
         xlab='Retention Time (min)',
         col='black',
         bty='n',
         cex.axis=1,
         cex.lab=1
    )
    lines(EIC_13C_df$times, EIC_13C_df$int*ratioC13_C12, 
          type='l',
          col = '#FF000088')
    title(paste('  ','LC-ESIMS.  EIC = ',round(cmp_mass, digits=3), cmp_name, "\n Max intensity = ", round(max_int, digits = -1), "Scans = ", scans_detected))
    abline(v=detected_peak_times, lty=3,col='gray48')
  }
  dev.off()
}
}

