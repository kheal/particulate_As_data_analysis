### THIS NEEDS TO BE FIXED TO INCORPORATE THE NEW CHANGES - LOOKING ONLY AT RTS of MS2s in appropriate rt range and opening up the MS1 window for the MS2s to 0.1


# Load libraries ---
library(tidyverse)
library(xcms)
library(zoo)
library(fuzzyjoin)
library(CluMSID)
library(Rdisop)

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
found_icap_peaks_filename <- "Intermediates/AsLipids_ICPpeakareas.csv"
icp_files_to_analyze <- c("2020_10_26_Heal_AsLipidswB12_1.csv", 
                      "2020_10_26_Heal_AsLipidswB12_2.csv", 
                      "2020_10_26_Heal_AsLipidswB12_3.csv")
sample_matcher_filename <- "MetaData/QE_ICAP_samplematcher.csv"
location_of_QEfiles <- "RawDat/20201026_QE_secondLipidRun"
files_to_analyze <- c("2020_10_26_Heal_AsLipidswB12_1.csv", 
                      "2020_10_26_Heal_AsLipidswB12_2.csv", 
                      "2020_10_26_Heal_AsLipidswB12_3.csv")
MS2_spectra_filename <- "MetaData/AsLipidDatabase/MS2_Fragments.csv"



#Load file of the samples to analyze, load in detected peaks
sample_matcher <- read_csv(sample_matcher_filename) %>%
  filter(str_detect(ICP_sampID, "crude|elute|KM")) %>%
  filter(!str_detect(ICP_sampID, "lank|filter")) 
peaks_all <- read_csv(found_icap_peaks_filename)
MS2_library <- read_csv(MS2_spectra_filename) %>%
  mutate(mz = sapply(Formula, 
                     function(x)getMolecule(x)$exactmass))


#Loop here!!
for (j in 1:length(sample_matcher$QE_filename)){

#Define appropriate QE file name, ICP filename, savefile name, get correct db
QE_file <- paste0(location_of_QEfiles, "/", sample_matcher$QE_filename[j])
ICP_sampID <- sample_matcher$ICP_sampID[j]
if(str_detect(QE_file, "highmass")){
  as_db <- read_csv(as_database_filename_highmass) %>% arrange(mz)
} else {as_db <- read_csv(as_database_filename_lowmass) %>% arrange(mz)}
savefile <- paste("Intermediates/Targeted_search/", gsub(".mzXML","",sample_matcher$QE_filename[j]),"_targeted.pdf",sep="")

# Read in Orbitrap file
mzxcms_raw <- xcmsRaw(QE_file,profstep=0.01,profmethod="bin",profparam=list(),includeMSn=FALSE,mslevel=NULL, scanrange=NULL)

# Read in MS2 file
savedfile <- str_replace(QE_file, ".mzXML", '_mergedMS2.rds')
my_spectra <- readRDS(savedfile)

# Grab ICPMS peaks from already written out file
peaks <- peaks_all %>%
  filter(sampID == sample_matcher$ICP_sampID[j]) %>%
  filter(peak_times > 120) # don't bother looking at the early peaks

# Get a list of possible good lipids based on criteria at top of script, plot each and ask if we should keep it as a possibility
possible_lipids <- c()

for(i in 1:nrow(as_db)){
  timerange <- c(500, 2500)
  cmp_mass <- as_db$mz[i]
  cmp_mass_C13 <- as_db$mz_13C[i]
  ratioC13_C12 <- as_db$ratio_12Cto13C[i]
  cmp_name <- as_db$Lipid_Name[i]
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
  { possible_lipids[i] <- as_db$EmpiricalFormula[i]}
}
possible_lipids <- na.omit(possible_lipids)   
print(paste(length(possible_lipids), "possible hits"))

# Make plots of all the possibles
as_db_culled <- as_db %>%
  filter(EmpiricalFormula %in% possible_lipids)
pdf(savefile, width = 7, height = 10)
if (length(possible_lipids) > 0){
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
    
    # Grab the MS2 data
    MS2s <- getSpectrum(my_spectra, "precursor", cmp_mass, mz.tol = 1E-03)
    if (!is.null(MS2s)){
    for (k in 1:length(MS2s)){
      if (length(MS2s) >1){
      MS2_time <- MS2s[[k]]@rt/60
      MS2_precursormass <- MS2s[[k]]@precursor
      MS2_to_plot <- tibble(mz = as.numeric(MS2s[[k]]@spectrum[,1]),
        intensity = as.numeric(MS2s[[k]]@spectrum[,2])) %>%
        distance_left_join(MS2_library %>% select(mz, Formula) %>% unique(), max_dist = 0.02)
      MS2_to_plot_withAs <- tibble(mz = as.numeric(MS2s[[k]]@spectrum[,1]),
                            intensity = as.numeric(MS2s[[k]]@spectrum[,2])) %>%
        distance_left_join(MS2_library %>% select(mz, Formula) %>% unique(), max_dist = 0.02) %>%
        filter(!is.na(Formula))} else
        {
        MS2_time <- MS2s@rt/60
        MS2_precursormass <- MS2s@precursor
        MS2_to_plot <- tibble(mz = as.numeric(MS2s@spectrum[,1]),
                              intensity = as.numeric(MS2s@spectrum[,2])) %>%
          distance_left_join(MS2_library %>% select(mz, Formula) %>% unique(), max_dist = 0.02)
        MS2_to_plot_withAs <- tibble(mz = as.numeric(MS2s@spectrum[,1]),
                                     intensity = as.numeric(MS2s@spectrum[,2])) %>%
          distance_left_join(MS2_library %>% select(mz, Formula) %>% unique(), max_dist = 0.02) %>%
          filter(!is.na(Formula))
          }
      
      par(mfrow=c(3,1))
      
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
      abline(v=MS2_time, lty=1,lwd = 2, col='green')
      
      
      plot(MS2_to_plot$mz.x, MS2_to_plot$intensity,
           type = 'h', 
          # xlim=c(8,35),
           ylab='Intensity',
           xlab='m/z',
           col='black',
           bty='n',
           cex.axis=1,
           cex.lab=1)
      points(MS2_to_plot_withAs$mz.x, MS2_to_plot_withAs$intensity,
           type = 'h', col = 'green', lwd = 2) 
      title(paste('MS2 spectra for mass', round(MS2_precursormass, digits = 4),
                  "at", round(MS2_time, digits = 1), "minutes"))
      
    }
     } else {
    
    # This is where the plots are made
    par(mfrow=c(3,1))
    
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
    abline(v=detected_peak_times, lty=3,col='orange')
    
    plot.new()

    }
  }
  dev.off()
}
}

