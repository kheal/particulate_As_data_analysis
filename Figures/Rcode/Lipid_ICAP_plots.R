library(tidyverse)
library(here)
library(cowplot)
theme_set(theme_cowplot())
library(stats)

# Name files ----
found_icap_peaks_filename <- "Intermediates/AsLipids_ICPpeakareas.csv"
sample_matcher_filename <- "MetaData/QE_ICAP_samplematcher.csv"

# Plot out files from first fileset -----
#Load file of the samples to analyze, load in detected peaks
sample_matcher <- read_csv(sample_matcher_filename) %>%
  filter(str_detect(ICP_sampID, "crude")) %>%
  filter(!str_detect(ICP_sampID, "lank|filter")) %>%
  select(ICP_sampID) %>% unique()
peaks_all <- read_csv(found_icap_peaks_filename)

plot_list <- list()

#Plot ALOHA with its detected peaks
for (i in 1:5){
ICP_sampID <- sample_matcher$ICP_sampID[i]
ICPdata <- read_csv(paste0("Intermediates/Baseline_CSVs/", 
                           sample_matcher$ICP_sampID[i], "_ICPwithbaseline.csv"))
peaks <- peaks_all %>%
  filter(sampID == sample_matcher$ICP_sampID[i]) %>%
  filter(peak_times > 120)

g <- ggplot(data = ICPdata) +
  geom_vline(data = peaks, aes(xintercept = peak_times/60), 
             alpha = 0.2, linetype = "dashed")+
  geom_line(aes(x = time/60,  y = intensity_smoothed), ) + 
 # geom_line(aes(x = time,  y = intenstiy_baselined)) + 
  ggtitle(str_remove(str_remove(sample_matcher$ICP_sampID[i], "_crude"), "Smp_"))+
  scale_x_continuous("Retention time (minutes)", limits = c(7, 35),
                     expand = c(0,0))+
  scale_y_continuous("Intensity",
                     limits = c(0, max(ICPdata$intensity_smoothed)*1.1),
                     expand = c(0,0))+
  theme(legend.position = "none")
g
plot_list[[i]] <-g
}

g2 <- plot_grid(plot_list[[1]], plot_list[[2]], plot_list[[3]],
                plot_list[[4]], plot_list[[5]], nrow = 5)
g2

save_plot("Figures/Preliminary/As_traces_ALOHA_BATS_ENTP.pdf", g2,
          base_width = 10, base_height = 12, units = "in")





