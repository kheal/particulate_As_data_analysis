#Load your libraries
library(tidyverse)
library(modelr)


# direct to folder with pivoted data
filenames <- list.files("RawDat/ICAP_run_20200901/Pivoted/", full.names = TRUE)
meta_dat_filename <- "MetaData/SampInfoNew.csv"

#load data, combine into one big DF
dat <- list()
for (i in 1:length(filenames)){
  dat[[i]] <- read_csv(filenames[i])
}
dat <- do.call(rbind, dat)

meta_dat <- read_csv(meta_dat_filename)%>%
  mutate(descript_ID = paste0(`Sample ID`, Replicate, `Fraction ID`)) %>%
  select(LC_ID, descript_ID, samplerun) %>%
  rename(sampleID = LC_ID)

meta_dat_more <- read_csv(meta_dat_filename)%>%
  mutate(descript_ID = paste0(`Sample ID`, Replicate, `Fraction ID`)) %>%
  rename(sampleID = LC_ID) %>%
  select(descript_ID, samplerun, sampleID, `Sample size (L)`, Reconst_volume_mL)

descript_ID <- unique(meta_dat$descript_ID)

dat_2 <- dat %>%
  left_join(meta_dat, by = c("sampleID", "samplerun"))


#TO DO- will need to calculate RFs at some point
# extract 75As signal
dat_as <- dat_2 %>%
  filter(element == "75As") %>%
  filter(time > 40 & time < 95) %>%
  filter(str_detect(sampleID, "Smp") | str_detect(sampleID, "Std"))

dat_as_pooled <- dat_as %>%
  group_by(sampleID, samplerun, descript_ID) %>%
  summarise(As_75_intensity = sum(intenstiy))

# extract 77'Se' signal
dat_77 <- dat_2 %>%
  filter(element == "77Se") %>%
  filter(time > 40 & time < 95) %>%
  filter(str_detect(sampleID, "Smp") | str_detect(sampleID, "Std"))

dat_77_pooled <- dat_77 %>%
  group_by(sampleID, samplerun, descript_ID) %>%
  summarise(Se_77_intensity = sum(intenstiy))

# correct 75As signal with 77Se signal
dat_as_pooled_77_corrected <- dat_as_pooled %>%
  left_join(dat_77_pooled, by = c(
    "sampleID", "samplerun", "descript_ID"))%>%
  mutate(Signal = As_75_intensity-Se_77_intensity*(0.7576/.2424)) 

# extract just the standards and calculate RFs
dat_corrected_stds <- dat_as_pooled_77_corrected %>%
  filter(str_detect(sampleID, "Std")) %>%
  mutate(sampleID = sampleID %>% str_replace_all("1uM", "1000nM")) %>%
  mutate(concentration_nM = sampleID %>% 
           str_extract("[^Std_nM]+") %>% as.numeric()) %>%
  filter(concentration_nM == 0 | concentration_nM >= 1)

# Calculating RF and intercept
df <- dat_corrected_stds %>%
  group_by(samplerun) %>%
  summarise(slope = coef(lm(Signal ~ concentration_nM, data = .))[["concentration_nM"]][1],
     intercept = lm(Signal ~ concentration_nM, data = .)[[1]][[1]][[1]]) %>%
  ungroup() %>%
  summarise(slope = mean(slope),
            intercept = mean(intercept))


# calculate concentrations
dat_quan_1 <- dat_as_pooled_77_corrected %>%
  mutate(nM_in_vial = (Signal - df$intercept[1])/df$slope[1]) %>%
  filter(str_detect(descript_ID, "ulk" )) %>%
  left_join(meta_dat_more, by = c("sampleID", "samplerun", "descript_ID")) %>%
  mutate(enviro_nM = nM_in_vial*Reconst_volume_mL/1000/`Sample size (L)`)
