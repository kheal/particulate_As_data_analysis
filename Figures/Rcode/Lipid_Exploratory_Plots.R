library(tidyverse)
library(here)
library(cowplot)
theme_set(theme_cowplot())
library(stats)

# Files to examine-----
files <- c("RawDat/20200116_icapdata_firstLipidRun/2020_01_16_Heal_AsLipidsRP_3.csv", 
           "RawDat/20200116_icapdata_firstLipidRun/2020_01_16_Heal_AsLipidsRP_4.csv",
           "RawDat/20200116_icapdata_firstLipidRun/2020_01_16_Heal_AsLipidsRP_5.csv")

samps_from_files1_crude <- c("ALOHA_crude","ALOHAPrefilter_crude", 
                             "PS107_crude")

samps_from_files1_elute <- c( "ALOHA_elute",  "ALOHAPrefilter_elute",
                       "PS107_elute_1", "ALOHABlank_crude_1")

samps_from_files2_crude <- c("PS2_crude","PS3_crude",
                       "BATS07_crude")

samps_from_files2_elute <- c( "PS2_elute", "PS3_elute",
                        "BATS07_elute")

# Plot out files from first fileset -----
icapdat1 <- read_delim(files[1], 
                       ";", escape_double = FALSE, trim_ws = TRUE, 
                       skip = 1)
plots <- list()

for (i in 1:length(samps_from_files1_crude)){
dat_all <- icapdat1 %>% select(X2, X3, X4, samps_from_files1_crude[i]) %>%
  rename(scan = X2, analyte  = X3, type = X4)
dat_as <- dat_all %>%
  filter(analyte == "75As" & type == "Y") %>%
  rename(As.intensity = samps_from_files1_crude[i]) %>% select(-analyte, -type)
dat_time <- dat_all %>%
  filter(analyte == "75As" & type == "Time") %>%
  rename(As.time = samps_from_files1_crude[i]) %>% select(-analyte, -type)
dat_time_as <- left_join(dat_as, dat_time, by = "scan")  %>%
  filter(As.time > 200 & As.time < 2000) %>% select(-scan)
dat_to_plot_crude <- tibble("Time" = runmed(dat_time_as$As.time, k = 7),
                       "Intensity" = runmed(dat_time_as$As.intensity, k = 7))

dat_all <- icapdat1 %>% select(X2, X3, X4, samps_from_files1_elute[i]) %>%
  rename(scan = X2, analyte  = X3, type = X4)
dat_as <- dat_all %>%
  filter(analyte == "75As" & type == "Y") %>%
  rename(As.intensity = samps_from_files1_elute[i]) %>% select(-analyte, -type)
dat_time <- dat_all %>%
  filter(analyte == "75As" & type == "Time") %>%
  rename(As.time = samps_from_files1_elute[i]) %>% select(-analyte, -type)
dat_time_as <- left_join(dat_as, dat_time, by = "scan")  %>%
  filter(As.time > 200 & As.time < 2000) %>% select(-scan)
dat_to_plot_elute <- tibble("Time" = runmed(dat_time_as$As.time, k = 7),
                            "Intensity" = runmed(dat_time_as$As.intensity, k = 7))

dat_to_plot <- dat_to_plot_elute %>% 
  mutate(ID = samps_from_files1_elute[i]) %>%
  bind_rows(dat_to_plot_crude
            %>% mutate(ID = samps_from_files1_crude[i]))

plots[[i]] <- ggplot(data = dat_to_plot, aes(x = Time,  y = Intensity, color = ID)) +
  geom_line() + 
    ggtitle(samps_from_files1_crude[i])+
  theme(legend.position = "none")
}



# Plot out files from second fileset -----
icapdat1 <- read_delim(files[2], 
                       ";", escape_double = FALSE, trim_ws = TRUE, 
                       skip = 1)
plots_2 <- list()

for (i in 1:length(samps_from_files2_crude)){
  dat_all <- icapdat1 %>% select(X2, X3, X4, samps_from_files2_crude[i]) %>%
    rename(scan = X2, analyte  = X3, type = X4)
  dat_as <- dat_all %>%
    filter(analyte == "75As" & type == "Y") %>%
    rename(As.intensity = samps_from_files2_crude[i]) %>% select(-analyte, -type)
  dat_time <- dat_all %>%
    filter(analyte == "75As" & type == "Time") %>%
    rename(As.time = samps_from_files2_crude[i]) %>% select(-analyte, -type)
  dat_time_as <- left_join(dat_as, dat_time, by = "scan")  %>%
    filter(As.time > 200 & As.time < 2000) %>% select(-scan)
  dat_to_plot_crude <- tibble("Time" = runmed(dat_time_as$As.time, k = 7),
                              "Intensity" = runmed(dat_time_as$As.intensity, k = 7))
  
  dat_all <- icapdat1 %>% select(X2, X3, X4, samps_from_files2_elute[i]) %>%
    rename(scan = X2, analyte  = X3, type = X4)
  dat_as <- dat_all %>%
    filter(analyte == "75As" & type == "Y") %>%
    rename(As.intensity = samps_from_files2_elute[i]) %>% select(-analyte, -type)
  dat_time <- dat_all %>%
    filter(analyte == "75As" & type == "Time") %>%
    rename(As.time = samps_from_files2_elute[i]) %>% select(-analyte, -type)
  dat_time_as <- left_join(dat_as, dat_time, by = "scan")  %>%
    filter(As.time > 200 & As.time < 2000) %>% select(-scan)
  dat_to_plot_elute <- tibble("Time" = runmed(dat_time_as$As.time, k = 7),
                              "Intensity" = runmed(dat_time_as$As.intensity, k = 7))
  
  dat_to_plot <- dat_to_plot_elute %>% 
    mutate(ID = samps_from_files2_elute[i]) %>%
    bind_rows(dat_to_plot_crude
              %>% mutate(ID = samps_from_files2_crude[i]))
  
  plots_2[[i]] <- ggplot(data = dat_to_plot, aes(x = Time,  y = Intensity, color = ID)) +
    geom_line() + 
    ggtitle(samps_from_files2_crude[i])+
    theme(legend.position = "none")
}



# grid out the plots
plot_save_1 <- plot_grid(plots[[1]], plots[[2]], plots[[3]], ncol = 1)
plot_save_2 <- plot_grid(plots_2[[1]], plots_2[[2]], plots_2[[3]], ncol = 1)

save_plot("Figures/Preliminary/As_traces_1.pdf", plot_save_1, base_width = 4, base_height = 8, units = "in")
save_plot("Figures/Preliminary/As_traces_2.pdf", plot_save_2, base_width = 4, base_height = 8, units = "in")

