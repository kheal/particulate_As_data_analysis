library(here)
library(tidyverse)
library(xtable)
library(lubridate)

# Name files
meta_dat_file <- "MetaData/SampleMetaDat.csv"
quan_dat_file <- "Intermediates/bulk_quantification.csv"

# Load files
meta_dat <- read_csv(meta_dat_file)
quan_dat <- read_csv(quan_dat_file)

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
  mutate(Sample = Sample %>% str_replace("_", "-"))
  mutate(replicate = str_extract(descript_ID, ".bulk") %>%
           str_remove("bulk")) 

quan_dat3 <- quan_dat2 %>%
  group_by(Sample) %>%
  summarise(`ave` = round(mean(enviro_nM)*1000, digits = 1),
            stdev = round(sd(enviro_nM)*1000,  digits = 1)) %>%
  mutate(`particulate As (pM)` = paste0(ave, "  $\\pm$", stdev)) %>%
  select(-ave, -stdev) 

table.latex <- xtable(quan_dat3, sanitize.text.function = identity)
caption(table.latex) <- "\\label{QuanSummaryTable}Summary of quantitative results of particulate arsenic."
align(table.latex) <- c("p{0.1cm}","p{1.2cm}", "p{2cm}")

print(table.latex, type="latex", 
      file="Tables/ManuscriptReady/QuantitativeSummary.tex",
      size="\\fontsize{7pt}{8pt}\\selectfont",
      sanitize.text.function = identity, include.rownames=FALSE)
