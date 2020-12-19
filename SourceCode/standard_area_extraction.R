####### This code finds peaks on ICPMS data and calculates each peak area, based on Jiwoon's code

# Load appropriate libraries ---
library(tidyverse)

# Files to read
files <- c("2020_10_26_Heal_AsLipidswB12_1.csv", 
             "2020_10_26_Heal_AsLipidswB12_2.csv", 
             "2020_10_26_Heal_AsLipidswB12_3.csv")
folder <- "RawDat/20201026_icapdata_secondLipidRun/Pivoted/"


# Read file
i=3
dat <- read_csv(paste0(folder, files[i])) %>%
  filter(element == "59Co") 
%>%
  filter(str_detect(sampleID, "Std_"))









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




