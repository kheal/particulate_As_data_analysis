# Load libraries ---
library(tidyverse)
library(xcms)
library(CluMSID)
library(CluMSIDdata)

#Name files
sample_matcher_filename <- "MetaData/QE_ICAP_samplematcher.csv"
location_of_QEfiles <- "RawDat/20201026_QE_secondLipidRun"

#Load file of the samples to analyze, load in detected peaks
sample_matcher <- read_csv(sample_matcher_filename) %>%
  filter(str_detect(ICP_sampID, "crude|elute|KM")) %>%
  filter(!str_detect(ICP_sampID, "lank|filter")) 

#Loop through here
# start at j = 12
for (j in 4:length(sample_matcher$QE_filename)){
  QE_file <- paste0(location_of_QEfiles, "/", sample_matcher$QE_filename[j])
  print(QE_file)
  savename <- str_replace(QE_file, ".mzXML", '_mergedMS2.rds')
  my_spectra <- extractMS2spectra(QE_file,  min_peaks = 1)
  print(paste0("original MS2 scans = ", length(my_spectra)))
  start.time <- Sys.time()
  my_spectra_merged <- mergeMS2spectra(my_spectra, mz_tolerance = 5E-06, 
                                       rt_tolerance = 20, peaktable = NULL, 
                                       exclude_unmatched = FALSE)
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
  print(paste0("merged MS2 scans = ", length(my_spectra_merged)))
  saveRDS(my_spectra_merged, savename)
}

#To reload the data, do this :)
#my_data <- readRDS(savename)
