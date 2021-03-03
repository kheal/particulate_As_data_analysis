library(here)
library(tidyverse)
library(xtable)
library(lubridate)

# Name files
meta_dat_file <- "MetaData/SampleMetaDat.csv"
quan_dat_file <- "Intermediates/bulk_quantification.csv"
under_hump_quan <- "Intermediates/total_75As_areas_and_concentrations.csv"
individual_quan <- "Intermediates/IDd_Lipids_Quanified.csv"

# Load files
meta_dat <- read_csv(meta_dat_file) %>% rename(Sample = `Sample ID`) %>%
  select(Sample, PC_uMperL_mean)
quan_dat <- read_csv(quan_dat_file)
under_hump_dat <- read_csv("Intermediates/total_75As_areas_and_concentrations.csv")
indiv_dat <- read_csv("Intermediates/IDd_Lipids_Quanified.csv")

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
  left_join(meta_dat) %>%
  mutate(nMAsperuMC = enviro_nM/PC_uMperL_mean) %>%
  group_by(Sample) %>%
  summarise(`ave` = round(mean(enviro_nM)*1000, digits = 1),
            stdev = round(sd(enviro_nM)*1000,  digits = 1),
            `aveperC` = round(mean(nMAsperuMC)*1000, digits = 1),
            stdevperC = round(sd(nMAsperuMC)*1000,  digits = 1)) %>%
  mutate(`particulate As (pM)` = paste0(ave, "  $\\pm$ ",stdev)) %>%
  mutate(`particulate As/C ($x$10\\textsuperscript{6})` = paste0(aveperC, "  $\\pm$ ", stdevperC)) %>%
  select(-(ave:stdevperC))

# Adding in under hump data
quan_dat4 <- quan_dat3 %>%
  left_join(under_hump_dat %>% 
              rename(Sample = `Sample ID`) %>%
              select(Sample, pMolAs_enviro), by = "Sample") %>%
  mutate(pMolAs_enviro = round(pMolAs_enviro, digits =1)) %>%
  rename(`total lipid As (pM)` = pMolAs_enviro)

# Adding all IDd lipids as classes
quan_dat_test <- indiv_dat %>%
  mutate(lipid_type = ifelse(str_detect(lipid_ID, "HC"), "HC", "AsSugPL")) %>%
  group_by(lipid_type, sampleID) %>%
  summarise(total = sum(pMolAs_enviro)) %>%
  pivot_wider(id_cols = sampleID, names_from = lipid_type, values_from = total) %>%
  rename(Sample = sampleID)
quan_dat5 <- quan_dat4 %>%
  left_join(quan_dat_test, by = "Sample")





 table.latex <- xtable(quan_dat3, sanitize.text.function = identity)
caption(table.latex) <- "\\label{QuanSummaryTable}Summary of quantitative results of particulate arsenic."
align(table.latex) <- c("p{0.1cm}","p{1.2cm}", "c{2cm}", "c{2cm}")

print(table.latex, type="latex", 
      file="Tables/ManuscriptReady/QuantitativeSummary.tex",
      size="\\fontsize{7pt}{8pt}\\selectfont",
      sanitize.text.function = identity, include.rownames=FALSE)
