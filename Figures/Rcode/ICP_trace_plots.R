library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(nationalparkcolors)

# Name files here
dat_file <-  "Intermediates/ICP_combined.csv"
meta_dat_file <- "MetaData/SampleMetaDat.csv"
peaks_to_integrate_filename <- "Intermediates/Quantified_Individual_Lipids_betterint.csv"


# Load data
dat <- read_csv(dat_file)
meta_dat <- read_csv(meta_dat_file)
peaks <- read_csv(peaks_to_integrate_filename) 


dat2 <- dat %>%
  left_join(meta_dat %>% select(`Sample ID`, sampleID), by = "sampleID") %>%
  filter(time > 600 & time < 2000)

max_intensity <- dat2 %>%
  group_by(`Sample ID`) %>%
  summarise(max_int = max(intensity_smoothed)*1.1)

peaks2 <- peaks %>% left_join(max_intensity, by = "Sample ID") %>%
  mutate(lipid_type = ifelse(str_detect(lipid_ID, "HC"), "AsHC",
                             ifelse(str_detect(lipid_ID, "AsSugPL"), "AsSugPL",
                                    ifelse(str_detect(lipid_ID, "AsSugPeL"), "AsSugPeL",
                                    ifelse(str_detect(lipid_ID, "mz"), "unknown, \nbut m/z known",
                                    "unknown")))))

peaks2$lipid_type = factor(peaks2$lipid_type, levels = c("AsHC", "AsSugPL", "AsSugPeL", "unknown, \nbut m/z known", "unknown"))

pal <- c(park_palette("Badlands", 4), 'grey80')
# Plot up the ICP traces
g <- ggplot()+
  geom_rect(data = peaks2, aes(xmin = rt_sec_lower/60, xmax = rt_sec_higher/60, 
                               ymin = 0, ymax = max_int,
                               fill = lipid_type),
            alpha = 0.4) +
  geom_area(data = dat2, aes(x = time/60, y = intenstiy_baselined), fill = "white")+
  geom_line(data = dat2, aes(x = time/60, y = intensity_smoothed)) +
  geom_hline(yintercept = 0)+
  facet_wrap(facets = vars(`Sample ID`), ncol = 1, scales = "free_y") +
  scale_y_continuous("Intensity", limits = c(0,NA), expand = c(0, NA)) +
  scale_x_continuous("Time (min)", limits = c(12.8, 33), expand = c(0, NA)) +
  scale_fill_manual(values = pal)+
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 7),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 7),
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.position = c(0.75, 0.95),
        legend.background = element_rect(fill="white", 
                                         size=0.5, linetype="solid", 
                                         colour ="white"))
  
g

save_plot("Figures/ManuscriptReady/ICP_traces.pdf", g, base_width = 4, base_height = 4.5)
