library(tidyverse)

# Set your file names
meta_dat_filename <- "MetaData/SampleMetaDat.csv"

# Load your files
samp_info <- read_csv(meta_dat_filename)

# Make it look nice
samp_info2 <- samp_info %>%
  select(`Sample ID`:L_digested) %>%
  arrange(`Sample ID`) %>%
  rename(Date = Date_Collected,
         Sample = `Sample ID`,
         `Filter set up (\\textmu m)` = Filter_Set_Up,
         `Volume (L)` = L_Filtered) %>%
  mutate(`Volume digested for bulk As (L)` = round(L_digested, digits = 2),
         `Volume extracted for arsenolipids (L)` = round(L_extracted, digits = 0)) %>%
  select(-L_digested, -L_extracted)

# Make QE datafiles look nice
samp_info3 <- samp_info %>%
  select(`Sample ID`, QE_crudeHM_ID:QE_eluteLM_ID) %>%
  rename(Sample = `Sample ID`,
         `Crude extract, high mass HRESIMS data` = QE_crudeHM_ID,
         `Crude extract, low mass HRESIMS data` = QE_crudeLM_ID,
         `Silica-cleaned extract, high mass HRESIMS data` = QE_eluteHM_ID,
         `Silica-cleaned extract, low mass HRESIMS data` = QE_eluteLM_ID)

# Mash them together
samp_info_4 <- samp_info2 %>%
  left_join(samp_info3, by = "Sample")

# Write it out for Supplemental Table 1
write_csv(samp_info_4, "Tables/ManuscriptReady/SuppTables/sample_descriptions_metadata.csv")
