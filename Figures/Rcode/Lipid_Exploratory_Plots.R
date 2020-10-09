library(tidyverse)
library(here)
library(cowplot)
theme_set(theme_cowplot())
library(stats)

# Files to examine
files <- c("RawDat/20200116_icapdata_firstLipidRun/2020_01_16_Heal_AsLipidsRP_3.csv", 
           "RawDat/20200116_icapdata_firstLipidRun/2020_01_16_Heal_AsLipidsRP_4.csv",
           "RawDat/20200116_icapdata_firstLipidRun/2020_01_16_Heal_AsLipidsRP_5.csv")

samps_from_files1 <- c("ALOHA_crude", "ALOHA_elute", "ALOHAPrefilter_crude", "ALOHAPrefilter_elute",
                       "PS107_crude", "PS107_elute_1", "ALOHABlank_crude_1")

samps_from_files2 <- c("PS2_crude", "PS2_elute", "PS3_crude", "PS3_elute",
                       "BATS07_crude", "BATS07_elute", "BATS03_crude")

# Plot out files from first fileset -----
icapdat1 <- read_delim(files[1], 
                       ";", escape_double = FALSE, trim_ws = TRUE, 
                       skip = 1)
plots <- list()

for (i in 1:length(samps_from_files1)){
dat_all <- icapdat1 %>% select(X2, X3, X4, samps_from_files1[i]) %>%
  rename(scan = X2, analyte  = X3, type = X4)
  
dat_as <- dat_all %>%
  filter(analyte == "75As" & type == "Y") %>%
  rename(As.intensity = samps_from_files1[i]) %>% select(-analyte, -type)

dat_time <- dat_all %>%
  filter(analyte == "75As" & type == "Time") %>%
  rename(As.time = samps_from_files1[i]) %>% select(-analyte, -type)

dat_time_as <- left_join(dat_as, dat_time, by = "scan")  %>%
  filter(As.time > 200 & As.time < 2000) %>% select(-scan)

dat_to_plot <- tibble("Time" = runmed(dat_time_as$As.time, k = 7),
                       "Intensity" = runmed(dat_time_as$As.intensity, k = 7))

plots[[i]] <- ggplot(data = dat_to_plot, aes(x = Time,  y = Intensity)) +
  geom_line() + 
  ggtitle(samps_from_files1[i])
}



# Plot out files from second fileset -----
icapdat1 <- read_delim(files[2], 
                       ";", escape_double = FALSE, trim_ws = TRUE, 
                       skip = 1)
plots_2 <- list()

for (i in 1:length(samps_from_files2)){
  dat_all <- icapdat1 %>% select(X2, X3, X4, samps_from_files2[i]) %>%
    rename(scan = X2, analyte  = X3, type = X4)
  
  dat_as <- dat_all %>%
    filter(analyte == "75As" & type == "Y") %>%
    rename(As.intensity = samps_from_files2[i]) %>% select(-analyte, -type)
  
  dat_time <- dat_all %>%
    filter(analyte == "75As" & type == "Time") %>%
    rename(As.time = samps_from_files2[i]) %>% select(-analyte, -type)
  
  dat_time_as <- left_join(dat_as, dat_time, by = "scan")  %>%
    filter(As.time > 200 & As.time < 2000) %>% select(-scan)
  
  dat_to_plot <- tibble("Time" = runmed(dat_time_as$As.time, k = 7),
                        "Intensity" = runmed(dat_time_as$As.intensity, k = 7))
  
  plots_2[[i]] <- ggplot(data = dat_to_plot, aes(x = Time,  y = Intensity)) +
    geom_line() + 
    ggtitle(samps_from_files2[i])
}

# grid out the plots
plot_save_1 <- plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], ncol = 1)
plot_save_2 <- plot_grid(plots_2[[1]], plots_2[[2]], plots_2[[3]], plots_2[[4]], plots_2[[5]], plots_2[[6]], ncol = 1)

# crudes only
plot_crudes <- plot_grid(plots_2[[7]], plots[[1]], plots[[5]], plots_2[[1]], plots_2[[3]], plots_2[[5]], ncol = 1)
plot_crudes

save_plot("Figures/Preliminary/As_traces.pdf", plot_crudes, base_width = 8, base_height = 12, units = "in")


# Print out crudes only