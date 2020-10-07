library(tidyverse)

# Set file names ----
files_to_pivot <- c("2020_09_01_KRH_AsTotal_1c.csv", 
                    "2020_09_01_KRH_AsTotal_2.csv", 
                    "2020_09_01_KRH_AsTotal_3.csv", 
                    "2020_09_01_KRH_AsTotal_4.csv", 
                    "2020_09_01_KRH_AsTotal_5.csv")
elements_to_analyze = c("75As", "77Se", "78Se","103Rh")


#Loop through each file to pivot and clean up a bit
for (i in 1:length(files_to_pivot)){
  
  dat_filename <-paste0("RawDat/ICAP_run_20200901/FixedNames/", files_to_pivot[i])
  
  # Read in first data file
  dat <- read_csv(dat_filename, skip = 1)
  
  # Mundge into a usable shape, add some much needed column names
  dat2 <- dat %>%
    select(-X1) %>%
    rename(scan = X2,
           element = X3,
           value = X4)
  
  #Initialize an empty list that we'll fill with dfs of each element we're analyzing
  dat_pivoted = list()
  
  #Loop through each of the elements to analyze
  for (j in 1:length(elements_to_analyze)){
    element_test = elements_to_analyze[j]
    
    #Get the intensity values for the element, pivot
    dat2_sub <- dat2 %>%
      filter(element == element_test) %>%
      filter(value == "Y") %>%
      select(-value) %>%
      pivot_longer(-(scan:element), names_to = "sampleID", values_to = "intenstiy")
    
    #Get the time values for the element, pivot
    dat2_sub_time <-  dat2 %>%
      filter(element == element_test) %>%
      filter(value == "Time") %>%
      select(-value)%>%
      pivot_longer(-(scan:element), names_to = "sampleID", values_to = "time")
    
    dat2_sub_long <- dat2_sub %>%
      left_join(dat2_sub_time, by = c("scan", "element", "sampleID"))
    dat_pivoted[[j]] <- dat2_sub_long
  }
  dat_pivoted_df <- do.call(rbind, dat_pivoted) %>%
    mutate(samplerun = files_to_pivot[i]) %>%
    mutate(sampleID = str_replace_all(sampleID, "-", "_"))
  
  output_filename = paste0("RawDat/ICAP_run_20200901/Pivoted/", files_to_pivot[i])
  write_csv(dat_pivoted_df, output_filename)
}

