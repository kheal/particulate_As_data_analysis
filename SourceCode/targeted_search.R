###### based on JP targeted_search_JP_191121.R

# Load libraries ---
library(tidyverse)
library(xcms)

# Set parameters
#timerange <- c(0,2000) #time range in seconds
min_height <- 50000        #minimum height to be plotted the max scan needs to be; 100000 is good for low mass
scannum <- 20           #number of scans that detected targeted mass
element <- "75As"
scanfrq <- 5  #scans in a row we need to keep before counting it as an actual scan
background <- 0 #intensity background to keep 

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


#Load file of the samples to analyze, load in all ICP data
sample_matcher <- read_csv(sample_matcher_filename)
ICPdata_all <- read_csv(paste0("RawDat/20201026_icapdata_secondLipidRun/Pivoted/", files_to_analyze[1]))  %>%
  rbind(read_csv(paste0("RawDat/20201026_icapdata_secondLipidRun/Pivoted/", files_to_analyze[2]))) %>%
  rbind(read_csv(paste0("RawDat/20201026_icapdata_secondLipidRun/Pivoted/", files_to_analyze[3])))
peaks_all <- read_csv(found_icap_peaks_filename)

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

# Subset ICPMS data
ICPdata <- ICPdata_all %>%
  filter(sampleID == sample_matcher$ICP_sampID[j]) %>%
  filter(element == "75As")

# Grab ICPMS peaks from already written out file
peaks <- peaks_all %>%
  filter(sampID == sample_matcher$ICP_sampID[j]) %>%
  filter(peak_times > 120)

# Get a list of possible good lipids based on criteria at top of script
possible_lipids <- c()
for(i in 1:nrow(as_db)){
  #i=1
  timerange <- c(500, 2500)
  cmp_mass <- as_db$mz[i]
  cmp_name <- as_db$Lipid_Name[i]
  detected_peak_times <- peaks$peak_times/60 
  timerange_detected_peaks <- c(min(peaks$peak_times)-60, 
                                max(peaks$peak_times)+60)/60

  #This checks if theres any signal between the top and bottom times
  mzxcms <- mzxcms_raw
  mzrange<-c(-0.005,0.005)+cmp_mass
  EIC<-rawEIC(mzxcms,mzrange)
  times<-mzxcms@scantime/60
  EIC_df<-data.frame(cbind(times,unlist(EIC[[2]])))
  colnames(EIC_df) <- c("times", "int")
  EIC_df$int[!with(rle(EIC_df$int > background),
               rep(values & lengths >= scanfrq, lengths))] <- 0
  scans_detected_df <- EIC_df %>% filter(between(times,
                                      min(peaks$peak_times)/60-60,
                                      max(peaks$peak_times)/60+60)) %>%
    filter(int > 0) 
  scans_detected <- length(scans_detected_df$int)
  max_int <- max(EIC_df$int)
  if (scans_detected > scannum & max_int > min_height){
  possible_lipids[i] <- as_db$Lipid_Name[i] 
  }
}
possible_lipids <- na.omit(possible_lipids)   
print(paste(length(possible_lipids), "possible hits"))
    

# Make plots of all the possibles
as_db_culled <- as_db %>%
  filter(Lipid_Name %in% possible_lipids)

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

  # This is where the plots are made
  par(mfrow=c(2,1))
  
  #This makes the ICP MS plot with the detected peaks
  plot(ICPdata$time/60,ICPdata$intenstiy,type='l',
       xlim=timerange_detected_peaks,
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
       xlim=timerange_detected_peaks,
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

