
# Load libraries ---
library(tidyverse)
library(xcms)
library(CluMSID)
library(CluMSIDdata)
library(Rdisop)
library(fuzzyjoin)

# Name files ----
sample_matcher_filename <- "MetaData/QE_ICAP_samplematcher.csv"
location_of_QEfiles <- "RawDat/20201026_QE_secondLipidRun"
MS2_spectra_filename <- "MetaData/AsLipidDatabase/MS2_Fragments.csv"
MS1s_library_filename_1 <- "MetaData/AsLipidDatabase/highmass_inclusionlist.csv"
MS1s_library_filename_2 <- "MetaData/AsLipidDatabase/lowmass_inclusionlist.csv"

#Load file of the samples to analyze, load in detected peaks
sample_matcher <- read_csv(sample_matcher_filename) %>%
  filter(str_detect(ICP_sampID, "crude|elute|KM")) %>%
  filter(!str_detect(ICP_sampID, "lank|filter")) 
MS2_library <- read_csv(MS2_spectra_filename)
MS1_library <- read_csv(MS1s_library_filename_1) %>%
  bind_rows(read_csv(MS1s_library_filename_2))

#Loop through sample_matcher$QE_filename[j] on j
j=25
QE_file <- paste0(location_of_QEfiles, "/", sample_matcher$QE_filename[j])
print(QE_file)
savedfile <- str_replace(QE_file, ".mzXML", '_mergedMS2.rds')
my_spectra <- readRDS(savedfile)

Matched_MS2 <- list()
#Loop through list of MS2_library on i
for (i in 1:length(MS2_library$Formula)){
Matched_MS2_small <- findFragment(my_spectra, 
                               getMolecule(MS2_library$Formula[i])$exactmass,
                               tolerance = 3E-06)
Matched_MS2 <- Matched_MS2 %>% append(Matched_MS2_small)
}

Matched_MS2 <- Matched_MS2 %>% unique()

MS1s <- map_chr(Matched_MS2, "precursor") %>% as_tibble() %>%
  mutate(mz = as.numeric(value)) %>% 
  distance_left_join(MS1_library, max_dist = 0.005)
print(QE_file)




srts <- map_chr(Matched_MS2, "retention time")
MS2 <- map_chr(Matched_MS2, paste("spectrum"))
Matched_MS2[[1]]@spectrum






