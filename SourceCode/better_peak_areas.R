# Load appropriate libraries ---
library(tidyverse)

# Set file names
peaks_to_integrate_filename <- "Intermediates/IDd_Lipids_Quanified.csv"

# Load up files
peaks_to_integrate <- read_csv(peaks_to_integrate_filename) %>% 
  select(sampleID, peak_times)

# Load up ICP data and save it out in a easy to use format for just the crude samples for easy to work with
dat1 <- "Intermediates/Baseline_CSVs/Smp_ALOHA_crude_ICPwithbaseline.csv"
dat2 <- "Intermediates/Baseline_CSVs/Smp_BATS_crude_ICPwithbaseline.csv"
dat3 <- "Intermediates/Baseline_CSVs/Smp_PS1_crude_ICPwithbaseline.csv"
dat4 <- "Intermediates/Baseline_CSVs/Smp_PS2_crude_ICPwithbaseline.csv"
dat5 <- "Intermediates/Baseline_CSVs/Smp_PS3_crude_ICPwithbaseline.csv"

dat_combo <- rbind(read_csv(dat1), read_csv(dat2)) %>%
  rbind(read_csv(dat3)) %>% rbind(read_csv(dat4)) %>% rbind(read_csv(dat5))

write_csv(dat_combo, "Intermediates/ICP_combined.csv")




# Lets try to get better areas for the detected peaks - THIS LOOKS WAY BETTER, DO THIS FOR the detected peaks that we've seen from the targeted list to get better areas.
deteced_peaks <- read_csv("Intermediates/AsLipids_ICPpeakareas_concen.csv") %>%
  select(peak_times, sampleID, int_area)
# Test for ALOHA try integrating this time : 1394.13135 +- 
test_time <- 1394.13135
peak_width <- 25
ICP_dat <- read_csv(files[1])
ICP_dat_to_integrate <- ICP_dat %>% 
  filter(time > (test_time-peak_width) & time < (test_time+peak_width)) %>%
  filter(intensity_corrected_smoothed > 0)
total_area <- AUC(ICP_dat_to_integrate$time, ICP_dat_to_integrate$intensity_corrected_smoothed)
total_area
plot(ICP_dat_to_integrate$time, ICP_dat_to_integrate$intensity_corrected_smoothed, type = "l")
