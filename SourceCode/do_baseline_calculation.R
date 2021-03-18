# Load appropriate libraries ---
library(baseline)
library(signal)
library(tidyverse)
detach("package:dplyr")
library(dplyr)

# Load functions ---
g <- function (x) {
  fx <- f(x)
  ifelse(fx > 0, fx, 0)
}

# Set files to analyze names ----
files_to_analyze <- c("2020_10_26_Heal_AsLipidswB12_1.csv", 
                      "2020_10_26_Heal_AsLipidswB12_2.csv", 
                    "2020_10_26_Heal_AsLipidswB12_3.csv")
element_to_analyze <- "75As" #This code only works for 75As
smp_tags <- c("_crude", "elute", "KM") #Only gets peak areas of crude samples
blank_tags <- c("blank", "prefilter")
co_rf_dat <- read_csv("Intermediates/CoRF.csv")
co_to_as <- 0.75 #this is the observed ratio in Co / As signal at > 30% B
meta_dat <- read_csv("MetaData/SampleMetaDat.csv") %>%
  select(`Sample ID`, sampleID, lipid_fraction_reconst_vol_mL, L_extracted)


# call up the pivoted ICPMS data and make into one file
ICPdata_all <- read_csv(paste0("RawDat/20201026_icapdata_secondLipidRun/Pivoted/", files_to_analyze[1]))  %>%
  rbind(read_csv(paste0("RawDat/20201026_icapdata_secondLipidRun/Pivoted/", files_to_analyze[2]))) %>%
  rbind(read_csv(paste0("RawDat/20201026_icapdata_secondLipidRun/Pivoted/", files_to_analyze[3])))

samples <- ICPdata_all %>% select(sampleID) %>% unique() %>%
  filter(str_detect(sampleID, paste(smp_tags, collapse = '|')))%>%
  filter(!str_detect(sampleID, paste(blank_tags, collapse = '|')))

peaks_df_all <- list()

# call out one sample at a time
# loop around each sample-----
for (j in 1:length(samples$sampleID)){
ICPdat_onesamp <- ICPdata_all %>%
  filter(sampleID == samples$sampleID[j])

ICPdat_onesamp_onelement <- ICPdat_onesamp %>%
  filter(element == element_to_analyze) %>%
  filter(time < 2500)
ICPdat_onesamp_onelement.t <- t(ICPdat_onesamp_onelement$intenstiy) %>% as.matrix()

# check out your baseline and save out the results-----
bc.als <- baseline(ICPdat_onesamp_onelement.t, 
                   lambda = 6, p = .01, maxit = 20, method='als') #l6, p0.01, maxit20
#plot(bc.als, xlim=c(0, 6500))
#text(x = 6000, y = 500, label = paste0(samples$sampleID[j], "\n", element_to_analyze))

# Add your baseline spectra to the same df, save out this information somewhere, but not now-----
ICPdat_onesamp_onelement_wbl <- ICPdat_onesamp_onelement %>%
  mutate(intenstiy_baselined = c(getBaseline(bc.als)),
         intensity_corrected =  c(getCorrected(bc.als))) %>%
  mutate(intensity_smoothed = sgolayfilt(intenstiy, p = 15),
         intensity_corrected_smoothed = sgolayfilt(intensity_corrected, p = 15))
write_csv(ICPdat_onesamp_onelement_wbl, paste0("Intermediates/Baseline_CSVs/", samples$sampleID[j], "_ICPwithbaseline.csv"))
