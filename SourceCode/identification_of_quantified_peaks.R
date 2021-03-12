
# Load appropriate libraries ---
library(tidyverse)
options(readr.num_columns = 0)
library(fuzzyjoin)

# Load up data---
quan_peaks <- read_csv("Intermediates/AsLipids_ICPpeakareas_concen.csv")
meta_dat <- read_csv("MetaData/SampleMetaDat.csv") %>%
  select(`Sample ID`, sampleID)
ids <- read_csv("MetaData/Targeted_visually_IDd_2.csv") %>% select(lipid_ID:sampleID)


# Clean up the data
quan_peaks2 <- quan_peaks %>%
  filter(sampleID %in% meta_dat$sampleID) %>%
  select(-sampleID) %>% rename(sampleID = `Sample ID`) %>%
  filter(peak_times > 750 & peak_times < 2000) %>%
  mutate(rt_min = peak_times/60)%>%
  difference_left_join(ids, by = "rt_min", max_dist = 0.6, distance_col = "rt_diff" ) %>%
  filter(sampleID.x == sampleID.y) %>%
  select(-sampleID.y) %>%
  rename(sampleID = sampleID.x)

# When there are double matches, return only the closer one in RT
quan_peaks3 <- quan_peaks2 %>%
  arrange(rt_diff) %>%
  group_by(sampleID, lipid_ID) %>%
  slice(1)

# Write out these quantified IDd lipids :)
write_csv(quan_peaks3, "Intermediates/IDd_Lipids.csv")
