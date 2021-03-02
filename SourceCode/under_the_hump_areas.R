####### This code finds peaks on ICPMS data and calculates each peak area, based on Jiwoon's code

# Load appropriate libraries ---
library(baseline)
library(signal)
library(tidyverse)
library(DescTools)

# at > 30% B, Co:As = 0.75
# Get RF of cobalamin standards
co_rf_dat <- read_csv("Intermediates/CoRF.csv")
co_to_as <- 0.75 #this is the observed ratio in Co / As signal at > 30% B
meta_dat <- read_csv("MetaData/SampleMetaDat.csv") %>%
  select(`Sample ID`, sampleID, lipid_fraction_reconst_vol_mL, L_extracted)

# Point to folder with smoothed ICP data
folder <- "Intermediates/Baseline_CSVs/"

# Load up each sample after baselining
files <- list.files(folder, full.names =  TRUE) %>% 
  str_subset(., "crude")

#Loop here
total_area_list <- list()

for (i in 1:length(files)){

ICP_dat <- read_csv(files[i])
sampleID <- ICP_dat$sampleID %>% unique() %>% as.character()
ICP_dat_to_integrate <- ICP_dat %>% filter(time > 1000 & time < 2500)
total_area <- AUC(ICP_dat_to_integrate$time, ICP_dat_to_integrate$intensity_smoothed)
total_area_list[[i]]<- tibble(sampleID = sampleID, area_under_hump = total_area)


}

total_areas <- do.call(bind_rows, total_area_list)

dat2 <- total_areas %>%
  mutate(nMolAs_in_vial = area_under_hump/co_to_as*(co_rf_dat$concentration[1]/co_rf_dat$area[1])) %>%
  left_join(meta_dat, by = "sampleID") %>%
  mutate(pMolAs_enviro = nMolAs_in_vial/10^3*lipid_fraction_reconst_vol_mL/L_extracted*10^3) %>%
  select(-lipid_fraction_reconst_vol_mL, -L_extracted)

# Write out these results
write_csv(dat2, "Intermediates/total_75As_areas_and_concentrations.csv")





