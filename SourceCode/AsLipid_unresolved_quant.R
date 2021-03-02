####### This code finds peaks on ICPMS data and calculates each peak area, based on Jiwoon's code

# Load appropriate libraries ---
library(tidyverse)

# Get RF of cobalamin standards
co_rf_dat <- read_csv("Intermediates/CoRF.csv")
co_to_as <- 0.75 #this is the observed ratio in Co / As signal at > 30% B
meta_dat <- read_csv("MetaData/SampleMetaDat.csv") %>%
  select(`Sample ID`, sampleID, lipid_fraction_reconst_vol_mL, L_extracted)

# Read in data
dat <- read_csv("Intermediates/total_75As_areas.csv")

# Get in vial concentrations!
dat2 <- dat %>%
  mutate(nMolAs_in_vial = area_under_hump/co_to_as*(co_rf_dat$concentration[1]/co_rf_dat$area[1])) %>%
  left_join(meta_dat, by = "sampleID") %>%
  mutate(pMolAs_enviro = nMolAs_in_vial/10^3*lipid_fraction_reconst_vol_mL/L_extracted*10^3) %>%
  select(-lipid_fraction_reconst_vol_mL, -L_extracted)

# Write out these results
write_csv("Intermediates/")
  



