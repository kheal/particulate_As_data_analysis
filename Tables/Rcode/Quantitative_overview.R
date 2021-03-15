# TO DO: breaks about at line 70, pick it up there soon

library(here)
library(tidyverse)
library(xtable)
library(lubridate)

# Name files
meta_dat_file <- "MetaData/SampleMetaDat.csv"
quan_dat_file <- "Intermediates/bulk_quantification.csv"
under_hump_file <- "Intermediates/total_75As_areas_and_concentrations.csv"
individual_file <- "Intermediates/Quantified_Individual_Lipids_betterint.csv"

# Load files
meta_dat <- read_csv(meta_dat_file) %>% rename(Sample = `Sample ID`) %>%
  select(Sample, PC_uMperL_mean, PC_uMperL_stdev, sampleID, Latitude, Longitude)
quan_dat <- read_csv(quan_dat_file)
under_hump_dat <- read_csv(under_hump_file)
indiv_dat <- read_csv(individual_file)

#Mundge sample names
quan_dat2 <- quan_dat %>%
  filter(!str_detect(descript_ID, "0.3")) %>%
  filter(!str_detect(descript_ID, "blank")) %>%
  mutate(Sample = ifelse(str_detect(descript_ID, "ETNP_PS"), 
                              str_extract(descript_ID, "ETNP_PS\\d"), 
                              ifelse(str_detect(descript_ID, "BATS_0.7"),
                                     "BATS", 
                                     ifelse(str_detect(descript_ID, "ALOHA.bulk"),
                                            "ALOHA", NA)))) %>%
  mutate(Sample = Sample %>% str_replace("_", "-")) %>%
  mutate(replicate = str_extract(descript_ID, ".bulk") %>%
           str_remove("bulk")) 

# Get the pAs, pC and pAs/C in correct format
quan_dat3 <- quan_dat2 %>%
  left_join(meta_dat %>% select(Sample, PC_uMperL_mean, PC_uMperL_stdev), by = "Sample") %>%
  group_by(Sample, PC_uMperL_mean, PC_uMperL_stdev) %>%
  summarise(pAs_pM_ave = mean(enviro_nM)*1000,
            pAs_pM_stdev = sd(enviro_nM)*1000) %>%
  mutate(pAs_per_Cx106_ave = pAs_pM_ave/PC_uMperL_mean,
         pAs_per_Cx106_stdev = pAs_per_Cx106_ave*sqrt((pAs_pM_stdev/pAs_pM_ave)^2 + 
                                                        (PC_uMperL_stdev/PC_uMperL_mean)^2 )) %>%
  mutate(`particulate As (pM)` = paste0(signif(pAs_pM_ave, digits = 2), " ± ",signif(pAs_pM_stdev, digits = 2))) %>%
  mutate(`particulate C (μM)` = paste0(signif(PC_uMperL_mean, digits = 2), " ± ",signif(PC_uMperL_stdev, digits = 2))) %>%
  mutate(`particulate As/C (x10^6)` = paste0(signif(pAs_per_Cx106_ave, digits = 2), " ± ",
                                         signif(pAs_per_Cx106_stdev, digits = 2))) 
#  select(-(PC_uMperL_mean:pAs_per_Cx106_stdev))

# Adding in under hump data
quan_dat4 <- quan_dat3 %>%
  left_join(under_hump_dat %>% 
              rename(Sample = `Sample ID`) %>%
              select(Sample, pMolAs_enviro), by = "Sample") %>%
  mutate(pMolAs_enviro = signif(pMolAs_enviro, digits =2)) %>%
  rename(`total As lipids (pM)` = pMolAs_enviro) %>%
  mutate(`particulate As as lipids (%)` = signif(`total As lipids (pM)`/pAs_pM_ave*100, digits = 0))

# Adding all IDd lipids as classes
quan_dat_test <- indiv_dat %>%
  mutate(lipid_type = ifelse(str_detect(lipid_ID, "HC"), "AsHC",
                             ifelse(str_detect(lipid_ID, "AsSugPL"), "AsSugPL",
                                    ifelse(str_detect(lipid_ID, "AsSugPeL"), "AsSugPeL",
                                           ifelse(str_detect(lipid_ID, "mz"), "unknown, but m/z known",
                                                  "unknown")))))%>%
  group_by(lipid_type, sampleID) %>%
  summarise(total = sum(pMolAs_enviro)) %>%
  pivot_wider(id_cols = sampleID, names_from = lipid_type, values_from = total) %>%
  left_join(meta_dat %>% select(Sample, sampleID), by = "sampleID") %>%
  select(-sampleID)

quan_dat5 <- quan_dat4 %>%
  left_join(quan_dat_test, by = "Sample") %>%
  mutate(AsSugPL = signif(AsSugPL, digits =2),
         AsHC = signif(AsHC, digits =2),
         AsSugPeL = signif(AsSugPeL, digits = 2),
         `unknown, but with m/z known` = signif(`unknown, but m/z known`, digits = 2))%>%
  rename(`AsSugPL (pM)`= AsSugPL,
         `AsHC (pM)` = AsHC,
         `AsSugPeL (pM)` = AsSugPeL)
  
# Adding Lat and Long
quan_dat6 <- quan_dat5 %>% ungroup() %>%
  left_join(meta_dat %>% select(Sample, Latitude, Longitude), by = "Sample") %>%
  select(Sample, Latitude, Longitude, everything()) %>%
  select(-(PC_uMperL_mean)) %>%
  select(-(PC_uMperL_stdev:pAs_per_Cx106_stdev)) %>%
  select(-(unknown:`unknown, but with m/z known`))

write_excel_csv(quan_dat6, "Tables/ManuscriptReady/QuantitativeSummary.csv")

