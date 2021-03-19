library(tidyverse)
library(xcms)
library(cowplot)
theme_set(theme_cowplot())
library(Rdisop)
library(fuzzyjoin)
library(patchwork)
library(ggrepel)
library(CluMSID)
round_any = function(x, accuracy, f=round){f(x/ accuracy) * accuracy}

# Name your files ----
icap_dat_file <- "Intermediates/ICP_combined.csv"
as_db_high_filename <- "MetaData/AsLipidDatabase/highmass_inclusionlist.csv"
as_db_low_filename <- "MetaData/AsLipidDatabase/lowmass_inclusionlist.csv"
MS2_library_filename <- "MetaData/AsLipidDatabase/MS2_Fragments.csv"

# Load dat used by all----
icap_dat <- read_csv(icap_dat_file)
as_bd <- read_csv(as_db_high_filename) %>%
  bind_rows(read_csv(as_db_low_filename))
MS2_library <- read_csv(MS2_library_filename)%>%
  mutate(mz = sapply(Formula, 
                     function(x)getMolecule(x)$exactmass)) %>%
  filter(str_detect(Formula, "As"))

# Set theme and layout for all----
my_theme <- theme(strip.background = element_blank(),
                  plot.title = element_text(size = 9, hjust = 0.5),
                  axis.title = element_text(size = 8),
                  axis.text.x = element_text(size = 7),
                  axis.text.y = element_blank(),
                  legend.position = "none",
                  plot.margin = unit(c(3,3,3,3), "pt"))

layout <- "
AAABBBB
AAABBBB
AAABBBB
"
# Plot 1 = unknown 971 late in PS3-----
sample_name <- "Smp_PS3_crude"
sample_name_to_plot <- "ETNP-PS3"
lipid_name_to_plot <- "AsSugPeL970"
rt_lipid = 26.4
ESI_dat_file <- "RawDat/20201026_QE_secondLipidRun/201028_Smp_PS3_elute_highmass.mzXML"
MS2_dat_file <- "RawDat/20201026_QE_secondLipidRun/201028_Smp_PS3_elute_highmass_mergedMS2.rds"
mass = 971.5540
icap_dat_sub <- icap_dat %>% filter(sampleID == sample_name) %>%
  filter(time > 8*60)

ESI_dat <- xcmsRaw(ESI_dat_file,profstep=0.01, profmethod="bin",
                   profparam=list(),
                   includeMSn=FALSE, mslevel=NULL, scanrange=NULL)
mzrange<-c(-0.01,0.01)+mass
EIC<-rawEIC(ESI_dat,mzrange)
rt_min<-ESI_dat@scantime/60
EIC_df<-data.frame(cbind(rt_min,unlist(EIC[[2]])))
colnames(EIC_df) <- c("rt_min", "intensity")

MS2_dat <- readRDS(MS2_dat_file)
MS2s_rt <- getSpectrum(MS2_dat, "rt", rt_lipid*60, rt.tol = 30)
MS2s_mz <- getSpectrum(MS2s_rt,  "precursor", 
                       mass, 
                       mz.tol = .2)
MS2_time <- MS2s_mz@rt/60
MS2_precursormass <- MS2s_mz@precursor
MS2s <- tibble(mz = as.numeric(MS2s_mz@spectrum[,1]),
               intensity = as.numeric(MS2s_mz@spectrum[,2])) %>%
  mutate(mz_round = round_any(mz, accuracy = 0.05)) %>%
  group_by(mz_round) %>%
  summarise(mz = weighted.mean(mz, intensity),
            intensity = sum(intensity))%>%
  distance_left_join(MS2_library %>%
                       select(mz, Formula) %>%
                       unique(), max_dist = 0.02, by = "mz") %>%
  mutate(mz = mz.x) %>%
  mutate(has_As = !is.na(Formula)) %>%
  mutate(rt = MS2_time)
MS2s_top <- MS2s %>%
  mutate(rank = rank(-intensity)) %>%
  filter(((has_As == TRUE) & (rank < 20) )| (rank < 5) | (mz> 450 & rank < 18))

p1a <- ggplot() +
  geom_line(data = EIC_df, aes(x = rt_min, y = intensity))+
  geom_line(data = icap_dat_sub, 
            aes(x = time/60, y = intensity_smoothed*1000), color = "grey80",
            alpha = 0.6, size = 0.3)+
  scale_x_continuous("Time (min)", limits = c(25, 32)) +
  scale_y_continuous("Intensity", limits = c(0, NA), expand = c(0, NA)) +
  labs(title = paste0(lipid_name_to_plot, " in ", sample_name_to_plot))+
  my_theme+
  theme(axis.title.x = element_blank())

p1b <- ggplot(data = MS2s, aes(x = mz,xend = mz, 
                               y = 0, yend = intensity)) +
  geom_segment(data = MS2s %>% 
                 filter(has_As == TRUE), 
               color = 'coral1', 
               size = 1)+
  geom_segment()+
  geom_hline(yintercept = 0)+
  geom_label_repel(data = MS2s_top, 
                   aes(x = mz, y = intensity, 
                       label = round(mz, digits = 3)), 
                   nudge_y = 0.2*max(MS2s_top$intensity),
                   size = 1.8, color = 'black',  min.segment.length = 0.35,
                   segment.color = 'grey80',
                   direction = "both")+
  scale_x_continuous("m/z", expand = c(0,NA), limits = c(50, 1002))+
  scale_y_continuous(element_blank(),  
                     limits=c(0, max(MS2s$intensity*1.15)), expand = c(0,NA))+  
  my_theme+
  theme(axis.title.x = element_blank())


p1 <- p1a + p1b + 
  plot_layout(design = layout)
p1


# Plot 2 = unknown 973 late in PS2-----
sample_name <- "Smp_PS3_crude"
sample_name_to_plot <- "ETNP-PS3"
lipid_name_to_plot <- "AsSugPeL972"
rt_lipid = 27
ESI_dat_file <- "RawDat/20201026_QE_secondLipidRun/201028_Smp_PS3_elute_highmass.mzXML"
MS2_dat_file <- "RawDat/20201026_QE_secondLipidRun/201028_Smp_PS3_elute_highmass_mergedMS2.rds"
mass = 973.5679
icap_dat_sub <- icap_dat %>% filter(sampleID == sample_name) %>%
  filter(time > 8*60)


ESI_dat <- xcmsRaw(ESI_dat_file,profstep=0.01, profmethod="bin",
                   profparam=list(),
                   includeMSn=FALSE, mslevel=NULL, scanrange=NULL)
mzrange<-c(-0.01,0.01)+mass
EIC<-rawEIC(ESI_dat,mzrange)
rt_min<-ESI_dat@scantime/60
EIC_df<-data.frame(cbind(rt_min,unlist(EIC[[2]])))
colnames(EIC_df) <- c("rt_min", "intensity")

MS2_dat <- readRDS(MS2_dat_file)
MS2s_rt <- getSpectrum(MS2_dat, "rt", rt_lipid*60, rt.tol = 20)
MS2s_mz <- getSpectrum(MS2s_rt,  "precursor", 
                       mass, 
                       mz.tol = .1)
MS2_time <- MS2s_mz@rt/60
MS2_precursormass <- MS2s_mz@precursor
MS2s <- tibble(mz = as.numeric(MS2s_mz@spectrum[,1]),
               intensity = as.numeric(MS2s_mz@spectrum[,2])) %>%
  mutate(mz_round = round_any(mz, accuracy = 0.05)) %>%
  group_by(mz_round) %>%
  summarise(mz = weighted.mean(mz, intensity),
            intensity = sum(intensity))%>%
  distance_left_join(MS2_library %>%
                       select(mz, Formula) %>%
                       unique(), max_dist = 0.02, by = "mz") %>%
  mutate(mz = mz.x) %>%
  mutate(has_As = !is.na(Formula)) %>%
  mutate(rt = MS2_time)
MS2s_top <- MS2s %>%
  mutate(rank = rank(-intensity)) %>%
  filter(((has_As == TRUE) & (rank < 20) )| (rank < 5) | (mz> 450 & rank < 20))

p2a <- ggplot() +
  geom_line(data = EIC_df, aes(x = rt_min, y = intensity))+
  geom_line(data = icap_dat_sub, 
            aes(x = time/60, y = intensity_smoothed*1000), color = "grey80",
            alpha = 0.6, size = 0.3)+
  scale_x_continuous("Time (min)", limits = c(25, 32)) +
  scale_y_continuous("Intensity", limits = c(0, NA), expand = c(0, NA)) +
  labs(title = paste0(lipid_name_to_plot, " in ", sample_name_to_plot))+
  my_theme+
  theme(axis.title.x = element_blank())

p2b <- ggplot(data = MS2s, aes(x = mz,xend = mz, 
                               y = 0, yend = intensity)) +
  geom_segment(data = MS2s %>% 
                 filter(has_As == TRUE), 
               color = 'coral1', 
               size = 1)+
  geom_segment()+
  geom_hline(yintercept = 0)+
  geom_label_repel(data = MS2s_top, 
                   aes(x = mz, y = intensity, 
                       label = round(mz, digits = 3)), 
                   nudge_y = 0.2*max(MS2s_top$intensity),
                   size = 1.8, color = 'black',  min.segment.length = 0.35,
                   segment.color = 'grey80',
                   direction = "both")+
  scale_x_continuous("m/z", expand = c(0,NA), limits = c(50, 1002))+
  scale_y_continuous(element_blank(),  
                     limits=c(0, max(MS2s$intensity*1.15)), expand = c(0,NA))+  
  my_theme+
  theme(axis.title.x = element_blank())


p2 <- p2a + p2b + 
  plot_layout(design = layout)
p2


# Plot 3 = unknown 999 late in PS2-----
sample_name <- "Smp_PS3_crude"
sample_name_to_plot <- "ETNP-PS3"
lipid_name_to_plot <- "AsSugPeL998"
rt_lipid = 27.6
ESI_dat_file <- "RawDat/20201026_QE_secondLipidRun/201028_Smp_PS3_elute_highmass.mzXML"
MS2_dat_file <- "RawDat/20201026_QE_secondLipidRun/201028_Smp_PS3_elute_highmass_mergedMS2.rds"
mass = 999.5884
icap_dat_sub <- icap_dat %>% filter(sampleID == sample_name) %>%
  filter(time > 8*60)


ESI_dat <- xcmsRaw(ESI_dat_file,profstep=0.01, profmethod="bin",
                   profparam=list(),
                   includeMSn=FALSE, mslevel=NULL, scanrange=NULL)
mzrange<-c(-0.01,0.01)+mass
EIC<-rawEIC(ESI_dat,mzrange)
rt_min<-ESI_dat@scantime/60
EIC_df<-data.frame(cbind(rt_min,unlist(EIC[[2]])))
colnames(EIC_df) <- c("rt_min", "intensity")

MS2_dat <- readRDS(MS2_dat_file)
MS2s_rt <- getSpectrum(MS2_dat, "rt", rt_lipid*60, rt.tol = 15)
MS2s_mz <- getSpectrum(MS2s_rt,  "precursor", 
                       mass, 
                       mz.tol = .1)
MS2_time <- MS2s_mz@rt/60
MS2_precursormass <- MS2s_mz@precursor
MS2s <- tibble(mz = as.numeric(MS2s_mz@spectrum[,1]),
               intensity = as.numeric(MS2s_mz@spectrum[,2])) %>%
  mutate(mz_round = round_any(mz, accuracy = 0.05)) %>%
  group_by(mz_round) %>%
  summarise(mz = weighted.mean(mz, intensity),
            intensity = sum(intensity))%>%
  distance_left_join(MS2_library %>%
                       select(mz, Formula) %>%
                       unique(), max_dist = 0.02, by = "mz") %>%
  mutate(mz = mz.x) %>%
  mutate(has_As = !is.na(Formula)) %>%
  mutate(rt = MS2_time)
MS2s_top <- MS2s %>%
  mutate(rank = rank(-intensity)) %>%
  filter(((has_As == TRUE) & (rank < 20) )| (rank < 5) | (mz> 450 & rank < 20))

p3a <- ggplot() +
  geom_line(data = EIC_df, aes(x = rt_min, y = intensity))+
  geom_line(data = icap_dat_sub, 
            aes(x = time/60, y = intensity_smoothed*1000), color = "grey80",
            alpha = 0.6, size = 0.3)+
  scale_x_continuous("Time (min)", limits = c(25, 32)) +
  scale_y_continuous("Intensity", limits = c(0, NA), expand = c(0, NA)) +
  labs(title = paste0(lipid_name_to_plot, " in ", sample_name_to_plot))+
  my_theme+
  theme(axis.title.x = element_blank())

p3b <- ggplot(data = MS2s, aes(x = mz,xend = mz, 
                               y = 0, yend = intensity)) +
  geom_segment(data = MS2s %>% 
                 filter(has_As == TRUE), 
               color = 'coral1', 
               size = 1)+
  geom_segment()+
  geom_hline(yintercept = 0)+
  geom_label_repel(data = MS2s_top, 
                   aes(x = mz, y = intensity, 
                       label = round(mz, digits = 3)), 
                   nudge_y = 0.2*max(MS2s_top$intensity),
                   size = 1.8, color = 'black',  min.segment.length = 0.35,
                   segment.color = 'grey80',
                   direction = "both")+
  scale_x_continuous("m/z", expand = c(0,NA), limits = c(50, 1002))+
  scale_y_continuous(element_blank(),  
                     limits=c(0, max(MS2s$intensity*1.15)), expand = c(0,NA))+  
  my_theme+
  theme(axis.title.x = element_blank())


p3 <- p3a + p3b + 
  plot_layout(design = layout)
p3


# Plot 4 = unknown 1001 late in PS2-----
sample_name <- "Smp_PS3_crude"
sample_name_to_plot <- "ETNP-PS3"
lipid_name_to_plot <- "AsSugPeL1000"
rt_lipid = 28.2
ESI_dat_file <- "RawDat/20201026_QE_secondLipidRun/201028_Smp_PS3_elute_highmass.mzXML"
MS2_dat_file <- "RawDat/20201026_QE_secondLipidRun/201028_Smp_PS3_elute_highmass_mergedMS2.rds"
mass = 1001.6018
icap_dat_sub <- icap_dat %>% filter(sampleID == sample_name) %>%
  filter(time > 8*60)


ESI_dat <- xcmsRaw(ESI_dat_file,profstep=0.01, profmethod="bin",
                   profparam=list(),
                   includeMSn=FALSE, mslevel=NULL, scanrange=NULL)
mzrange<-c(-0.01,0.01)+mass
EIC<-rawEIC(ESI_dat,mzrange)
rt_min<-ESI_dat@scantime/60
EIC_df<-data.frame(cbind(rt_min,unlist(EIC[[2]])))
colnames(EIC_df) <- c("rt_min", "intensity")

MS2_dat <- readRDS(MS2_dat_file)
MS2s_rt <- getSpectrum(MS2_dat, "rt", rt_lipid*60, rt.tol = 15)
MS2s_mz <- getSpectrum(MS2s_rt,  "precursor", 
                       mass, 
                       mz.tol = .1)
MS2_time <- MS2s_mz@rt/60
MS2_precursormass <- MS2s_mz@precursor
MS2s <- tibble(mz = as.numeric(MS2s_mz@spectrum[,1]),
               intensity = as.numeric(MS2s_mz@spectrum[,2])) %>%
  mutate(mz_round = round_any(mz, accuracy = 0.05)) %>%
  group_by(mz_round) %>%
  summarise(mz = weighted.mean(mz, intensity),
            intensity = sum(intensity))%>%
  distance_left_join(MS2_library %>%
                       select(mz, Formula) %>%
                       unique(), max_dist = 0.02, by = "mz") %>%
  mutate(mz = mz.x) %>%
  mutate(has_As = !is.na(Formula)) %>%
  mutate(rt = MS2_time)
MS2s_top <- MS2s %>%
  mutate(rank = rank(-intensity)) %>%
  filter(((has_As == TRUE) & (rank < 20) )| (rank < 5) | (mz> 450 & rank < 20))

p4a <- ggplot() +
  geom_line(data = EIC_df, aes(x = rt_min, y = intensity))+
  geom_line(data = icap_dat_sub, 
            aes(x = time/60, y = intensity_smoothed*1000), color = "grey80",
            alpha = 0.6, size = 0.3)+
  scale_x_continuous("Time (min)", limits = c(25, 32)) +
  scale_y_continuous("Intensity", limits = c(0, NA), expand = c(0, NA)) +
  labs(title = paste0(lipid_name_to_plot, " in ", sample_name_to_plot))+
  my_theme

p4b <- ggplot(data = MS2s, aes(x = mz,xend = mz, 
                               y = 0, yend = intensity)) +
  geom_segment(data = MS2s %>% 
                 filter(has_As == TRUE), 
               color = 'coral1', 
               size = 1)+
  geom_segment()+
  geom_hline(yintercept = 0)+
  geom_label_repel(data = MS2s_top, 
                   aes(x = mz, y = intensity, 
                       label = round(mz, digits = 3)), 
                   nudge_y = 0.2*max(MS2s_top$intensity),
                   size = 1.8, color = 'black',  min.segment.length = 0.35,
                   segment.color = 'grey80',
                   direction = "both")+
  scale_x_continuous("m/z", expand = c(0,NA), limits = c(50, 1002))+
  scale_y_continuous(element_blank(),  
                     limits=c(0, max(MS2s$intensity*1.15)), expand = c(0,NA))+  
  my_theme+
  theme(axis.title.x = element_blank())

p4 <- p4a + p4b + 
  plot_layout(design = layout)


# Put them all together ----
p6 <- p1/p2/p3/p4
p6

save_plot("Figures/ManuscriptReady/PS3_unknowns.pdf", p6, 
          base_width = 6, base_height = 7)

