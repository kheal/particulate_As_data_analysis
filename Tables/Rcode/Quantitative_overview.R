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
  select(Sample, PC_uMperL_mean, sampleID, Latitude, Longitude)
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
  left_join(meta_dat %>% select(Sample, PC_uMperL_mean), by = "Sample") %>%
  mutate(nMAsperuMC = enviro_nM/PC_uMperL_mean) %>%
  group_by(Sample) %>%
  summarise(`ave` = round(mean(enviro_nM)*1000, digits = 1),
            stdev = round(sd(enviro_nM)*1000,  digits = 1),
            `aveperC` = round(mean(nMAsperuMC)*1000, digits = 1),
            stdevperC = round(sd(nMAsperuMC)*1000,  digits = 1)) %>%
  mutate(`\\makecell{particulate As \\\\ (pM)}` = paste0(ave, "  $\\pm$ ",stdev)) %>%
  mutate(`\\makecell{particulate As/C \\\\ ($x$10\\textsuperscript{6})}`=
           paste0(aveperC, "  $\\pm$ ", stdevperC)) %>%
  select(-(ave:stdevperC))

# Adding in under hump data
quan_dat4 <- quan_dat3 %>%
  left_join(under_hump_dat %>% 
              rename(Sample = `Sample ID`) %>%
              select(Sample, pMolAs_enviro), by = "Sample") %>%
  mutate(pMolAs_enviro = round(pMolAs_enviro, digits =1)) %>%
  rename(`\\makecell{total As lipids \\\\ (pM)}` = pMolAs_enviro)

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
  mutate(AsSugPL = round(AsSugPL, digits =2),
         HC = round(HC, digits =2))%>%
  rename(`\\makecell{identified  \\\\ AsSugPL (pM)}`= AsSugPL,
         `\\makecell{identified  \\\\ AsHC (pM)}` = HC)
  
# Adding Lat and Long
quan_dat6 <- quan_dat5 %>%
  left_join(meta_dat %>% select(Sample, Latitude, Longitude), by = "Sample") %>%
  select(Sample, Latitude, Longitude, everything())


table.latex <- xtable(quan_dat6, sanitize.text.function = identity)
caption(table.latex) <- "\\label{QuanSummaryTable}Summary of quantitative results of particulate arsenic and arsenolipids."
align(table.latex) <- c("p{0.1cm}","p{1.2cm}", "c{1.1cm}", "c{1.1cm}", "c{2cm}", 
                        "c{2cm}", "c{1.6cm}", "c{1.6cm}", "c{1.6cm}")

print(table.latex, type="latex", 
      file="Tables/ManuscriptReady/QuantitativeSummary.tex",
      size="\\fontsize{7pt}{8pt}\\selectfont",
      sanitize.text.function = identity, include.rownames=FALSE)
