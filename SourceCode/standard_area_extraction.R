####### This code finds peaks on ICPMS data and calculates each peak area, based on Jiwoon's code

# Load appropriate libraries ---
library(tidyverse)

# Files to read
file<- c("2020_10_26_Heal_AsLipidswB12_3.csv")
folder <- "RawDat/20201026_icapdata_secondLipidRun/Pivoted/"


# Read file
dat <- read_csv(paste0(folder, file)) %>%
  filter(element == "59Co") %>%
  filter(str_detect(sampleID, "Std_"))

# Split into the two injections
dat1 <- dat %>%
  filter(str_detect(sampleID, "a")) %>% 
  filter(time < 130 & time > 80)

dat2 <- dat %>%
  filter(str_detect(sampleID, "b"))%>% 
  filter(time < 130 & time > 80)

dat1_integrated <- AUC(dat1$time, dat1$intenstiy)
dat2_integrated <- AUC(dat2$time, dat2$intenstiy)

RF <- tibble(area = mean(dat1_integrated, dat2_integrated), 
             element = "Co", concentration = 400,
             percentB = 30)

write_csv(RF, "Intermediates/CoRF.csv")




