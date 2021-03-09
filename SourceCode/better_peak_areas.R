# Load appropriate libraries ---
library(tidyverse)
library(DescTools)

# Set file names and constantss
peaks_to_integrate_filename <- "Intermediates/IDd_Lipids.csv"
peaks_to_integrate_filename2 <- "MetaData/Targeted_visually_IDd_2.csv"
meta_dat <- read_csv("MetaData/SampleMetaDat.csv") 
co_rf_dat <- read_csv("Intermediates/CoRF.csv")
co_to_as <- 0.75

# Load up files
# peaks_to_integrate <- read_csv(peaks_to_integrate_filename)  %>%
#   rename(`Sample ID` = sampleID) %>%
#   mutate(rt_min = rt_min.x) %>%
#   select(`Sample ID`, lipid_ID, rt_min) %>%
#   left_join(meta_dat %>% select(`Sample ID`, sampleID), by = "Sample ID") %>%
#   filter(!is.na(rt_min))

peaks_to_integrate <- read_csv(peaks_to_integrate_filename2)  %>%
  rename(`Sample ID` = sampleID) %>%
  select(`Sample ID`, lipid_ID, rt_sec, rt_sec_lower, rt_sec_higher) %>%
  left_join(meta_dat %>% select(`Sample ID`, sampleID), by = "Sample ID") %>%
  filter(!is.na(rt_sec_lower))

# Load up ICP data and save it out in a easy to use format for just the crude samples for easy to work with
dat1 <- "Intermediates/Baseline_CSVs/Smp_ALOHA_crude_ICPwithbaseline.csv"
dat2 <- "Intermediates/Baseline_CSVs/Smp_BATS_crude_ICPwithbaseline.csv"
dat3 <- "Intermediates/Baseline_CSVs/Smp_PS1_crude_ICPwithbaseline.csv"
dat4 <- "Intermediates/Baseline_CSVs/Smp_PS2_crude_ICPwithbaseline.csv"
dat5 <- "Intermediates/Baseline_CSVs/Smp_PS3_crude_ICPwithbaseline.csv"

dat_combo <- rbind(read_csv(dat1), read_csv(dat2)) %>%
  rbind(read_csv(dat3)) %>% rbind(read_csv(dat4)) %>% rbind(read_csv(dat5))

write_csv(dat_combo, "Intermediates/ICP_combined.csv")


# Loop around each item in peaks_to_integrate
areas = c()
for (i in 1:length(peaks_to_integrate$sampleID)){
ICP_to_integrate = dat_combo %>% filter(sampleID == peaks_to_integrate$sampleID[i]) %>% 
  filter(time > (peaks_to_integrate$rt_sec_lower[i]) & 
           time < (peaks_to_integrate$rt_sec_higher[i])) %>%
  filter(intensity_corrected_smoothed > 0)
areas[i] = AUC(ICP_to_integrate$time, ICP_to_integrate$intensity_corrected_smoothed)
#plot(ICP_to_integrate$time, ICP_to_integrate$intensity_corrected_smoothed, type = "l")
}

peaks_to_integrate <- peaks_to_integrate %>%
  mutate(area = areas)

peaks_to_integrate_concentrations <- peaks_to_integrate %>%
  mutate(nMolAs_in_vial = area/co_to_as*(co_rf_dat$concentration[1]/co_rf_dat$area[1])) %>%
  left_join(meta_dat %>% select(sampleID, lipid_fraction_reconst_vol_mL, L_extracted), by = "sampleID") %>%
  mutate(pMolAs_enviro = nMolAs_in_vial/10^3*lipid_fraction_reconst_vol_mL/L_extracted*10^3) %>%
  select(-lipid_fraction_reconst_vol_mL, -L_extracted)


write_csv(peaks_to_integrate_concentrations, "Intermediates/Quantified_Individual_Lipids_betterint.csv")


