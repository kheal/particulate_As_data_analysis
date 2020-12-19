####### This code finds peaks on ICPMS data and calculates each peak area, based on Jiwoon's code

# Load appropriate libraries ---
library(baseline)
library(signal)
library(tidyverse)
detach("package:dplyr")
library(dplyr)
library(DescTools)


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

write_csv(total_areas, "Intermediates/total_75As_areas.csv")
## TO DO NEXT - get area of Cobalt in good samples




