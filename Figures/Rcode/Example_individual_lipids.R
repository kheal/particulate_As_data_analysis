library(tidyverse)
library(xcms)
library(cowplot)
theme_set(theme_cowplot())
library(Rdisop)
library(fuzzyjoin)
library(patchwork)


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
                  legend.position = "none")

layout <- "
AAAABB
AAAABB
AAAABB
"

# Plot 1 = AsHC332 in PS3-----
sample_name <- "Smp_PS3_crude"
sample_name_to_plot <- "ETNP-PS3"
lipid_name_to_plot <- "AsHC332"
rt_lipid = 14
ESI_dat_file <- "RawDat/20201026_QE_secondLipidRun/201028_Smp_PS3_elute_lowmass.mzXML"
MS2_dat_file <- "RawDat/20201026_QE_secondLipidRun/201028_Smp_PS3_elute_lowmass_mergedMS2.rds"
mass = as_bd %>% filter(Lipid_Name == lipid_name_to_plot) %>% select(Lipid_Name, mz)
icap_dat_sub <- icap_dat %>% filter(sampleID == sample_name)


ESI_dat <- xcmsRaw(ESI_dat_file,profstep=0.01, profmethod="bin",
                                profparam=list(),
                                includeMSn=FALSE, mslevel=NULL, scanrange=NULL)
mzrange<-c(-0.01,0.01)+mass$mz[1]
EIC<-rawEIC(ESI_dat,mzrange)
rt_min<-ESI_dat@scantime/60
EIC_df<-data.frame(cbind(rt_min,unlist(EIC[[2]])))
colnames(EIC_df) <- c("rt_min", "intensity")

MS2_dat <- readRDS(MS2_dat_file)
MS2s_rt <- getSpectrum(MS2_dat, "rt", rt_lipid*60, rt.tol = 20)
MS2s_mz <- getSpectrum(MS2s_rt,  "precursor", 
                       mass$mz[1], 
                       mz.tol = .1)
MS2_time <- MS2s_mz@rt/60
MS2_precursormass <- MS2s_mz@precursor
MS2s <- tibble(mz = as.numeric(MS2s_mz@spectrum[,1]),
                    intensity = as.numeric(MS2s_mz@spectrum[,2])) %>%
  distance_left_join(MS2_library %>%
                       select(mz, Formula) %>%
                       unique(), max_dist = 0.02, by = "mz") %>%
  mutate(mz = mz.x) %>%
  mutate(has_As = !is.na(Formula)) %>%
  mutate(rt = MS2_time) %>%
  filter(mz < 500)
MS2s_top <- MS2s %>%
  mutate(rank = rank(-intensity)) %>%
  filter(((has_As == TRUE) & (rank < 6) )| (rank < 3))

p1a <- ggplot() +
   geom_line(data = icap_dat_sub, 
              aes(x = time/60, y = intensity_smoothed*700), color = "coral1",
             alpha = 0.5)+
  geom_line(data = EIC_df, aes(x = rt_min+.2, y = intensity))+
  scale_x_continuous("Time (min)", limits = c(10, 32)) +
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
  geom_text(data = MS2s_top, 
            aes(x = mz, y = intensity*1.1, 
                label = round(mz, digits = 4)), 
            size = 1.8, color = 'black',
            check_overlap = TRUE)+
  scale_x_continuous("m/z", expand = c(0,NA), limits = c(0, 400))+
  scale_y_continuous(element_blank(),  
                     limits=c(0, max(MS2s$intensity*1.15)), expand = c(0,NA))+  
  my_theme+
  theme(axis.title.x = element_blank())

p1 <- p1a + p1b + 
    plot_layout(design = layout)
p1


# Plot 2 = AsSugPL954 in PS1----
sample_name <- "Smp_PS1_crude"
sample_name_to_plot <- "ETNP-PS1"
lipid_name_to_plot <- "AsSugPL954"
rt_lipid = 23.1
ESI_dat_file <- "RawDat/20201026_QE_secondLipidRun/201028_Smp_PS1_elute_highmass.mzXML"
MS2_dat_file <- "RawDat/20201026_QE_secondLipidRun/201028_Smp_PS1_elute_highmass_mergedMS2.rds"
mass = as_bd %>% filter(Lipid_Name == lipid_name_to_plot) %>% select(Lipid_Name, mz)
icap_dat_sub <- icap_dat %>% filter(sampleID == sample_name)


ESI_dat <- xcmsRaw(ESI_dat_file,profstep=0.01, profmethod="bin",
                   profparam=list(),
                   includeMSn=FALSE, mslevel=NULL, scanrange=NULL)
mzrange<-c(-0.01,0.01)+mass$mz[1]
EIC<-rawEIC(ESI_dat,mzrange)
rt_min<-ESI_dat@scantime/60
EIC_df<-data.frame(cbind(rt_min,unlist(EIC[[2]])))
colnames(EIC_df) <- c("rt_min", "intensity")

MS2_dat <- readRDS(MS2_dat_file)
MS2s_rt <- getSpectrum(MS2_dat, "rt", rt_lipid*60, rt.tol = 20)
MS2s_mz <- getSpectrum(MS2s_rt,  "precursor", 
                       mass$mz[1], 
                       mz.tol = .1)
MS2_time <- MS2s_mz@rt/60
MS2_precursormass <- MS2s_mz@precursor
MS2s <- tibble(mz = as.numeric(MS2s_mz@spectrum[,1]),
               intensity = as.numeric(MS2s_mz@spectrum[,2])) %>%
  distance_left_join(MS2_library %>%
                       select(mz, Formula) %>%
                       unique(), max_dist = 0.02, by = "mz") %>%
  mutate(mz = mz.x) %>%
  mutate(has_As = !is.na(Formula)) %>%
  mutate(rt = MS2_time)
MS2s_top <- MS2s %>%
  mutate(rank = rank(-intensity)) %>%
  filter(((has_As == TRUE) & (rank < 20) )| (rank < 3))
p2a <- ggplot() +
  geom_line(data = icap_dat_sub, 
            aes(x = time/60, y = intensity_smoothed*1500), color = "coral1",
            alpha = 0.5)+
  geom_line(data = EIC_df, aes(x = rt_min+.1, y = intensity))+
  scale_x_continuous("Time (min)", limits = c(10, 32)) +
  scale_y_continuous("Intensity", limits = c(0, NA), expand = c(0, NA)) +
  labs(title = paste0(lipid_name_to_plot, " in ", sample_name_to_plot))+
  my_theme +
  theme(axis.title.x = element_blank())

p2b <- ggplot(data = MS2s, aes(x = mz,xend = mz, 
                               y = 0, yend = intensity)) +
  geom_segment(data = MS2s %>% 
                 filter(has_As == TRUE), 
               color = 'coral1', 
               size = 1)+
  geom_segment()+
  geom_hline(yintercept = 0)+
  geom_text(data = MS2s_top, 
            aes(x = mz, y = intensity+5000, 
                label = round(mz, digits = 4)), 
            size = 2, color = 'black',
            check_overlap = TRUE)+
  scale_x_continuous("m/z", limits = c(0, 450),
                     expand = c(0,0))+
  scale_y_continuous(element_blank(),  
                     limits=c(0, max(MS2s$intensity*1.05)), expand = c(0,NA))+  
  my_theme+
  theme(axis.title.x = element_blank())

p2 <- p2a + p2b + 
  plot_layout(design = layout)


# Plot 3 = AsSugPL1010 in BATS----
sample_name <- "Smp_BATS_crude"
sample_name_to_plot <- "BATS"
lipid_name_to_plot <- "AsSugPL1010"
rt_lipid = 25.7
ESI_dat_file <- "RawDat/20201026_QE_secondLipidRun/201028_Smp_BATS_elute_highmass.mzXML"
MS2_dat_file <- "RawDat/20201026_QE_secondLipidRun/201028_Smp_BATS_elute_highmass_mergedMS2.rds"
mass = as_bd %>% filter(Lipid_Name == lipid_name_to_plot) %>% select(Lipid_Name, mz)
icap_dat_sub <- icap_dat %>% filter(sampleID == sample_name) %>%
  filter(time > 8*60)


ESI_dat <- xcmsRaw(ESI_dat_file,profstep=0.01, profmethod="bin",
                   profparam=list(),
                   includeMSn=FALSE, mslevel=NULL, scanrange=NULL)
mzrange<-c(-0.01,0.01)+mass$mz[1]
EIC<-rawEIC(ESI_dat,mzrange)
rt_min<-ESI_dat@scantime/60
EIC_df<-data.frame(cbind(rt_min,unlist(EIC[[2]])))
colnames(EIC_df) <- c("rt_min", "intensity")

MS2_dat <- readRDS(MS2_dat_file)
MS2s_rt <- getSpectrum(MS2_dat, "rt", rt_lipid*60, rt.tol = 20)
MS2s_mz <- getSpectrum(MS2s_rt,  "precursor", 
                       mass$mz[1], 
                       mz.tol = .1)
MS2_time <- MS2s_mz@rt/60
MS2_precursormass <- MS2s_mz@precursor
MS2s <- tibble(mz = as.numeric(MS2s_mz@spectrum[,1]),
               intensity = as.numeric(MS2s_mz@spectrum[,2])) %>%
  distance_left_join(MS2_library %>%
                       select(mz, Formula) %>%
                       unique(), max_dist = 0.02, by = "mz") %>%
  mutate(mz = mz.x) %>%
  mutate(has_As = !is.na(Formula)) %>%
  mutate(rt = MS2_time)
MS2s_top <- MS2s %>%
  mutate(rank = rank(-intensity)) %>%
  filter(((has_As == TRUE) & (rank < 20) )| (rank < 3))

p3a <- ggplot() +
  geom_line(data = icap_dat_sub, 
            aes(x = time/60, y = intensity_smoothed*1500), color = "coral1",
            alpha = 0.5)+
  geom_line(data = EIC_df, aes(x = rt_min+.1, y = intensity))+
  scale_x_continuous("Time (min)", limits = c(10, 32)) +
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
  geom_text(data = MS2s_top, 
            aes(x = mz, y = intensity+5000, 
                label = round(mz, digits = 4)), 
            size = 2, color = 'black',
            check_overlap = TRUE)+
  scale_x_continuous("m/z", limits = c(0, 450),
                     expand = c(0,0))+
  scale_y_continuous(element_blank(),  
                     limits=c(0, max(MS2s$intensity*1.05)), expand = c(0,NA))+  
  my_theme+
  theme(axis.title.x = element_blank())

p3 <- p3a + p3b + 
  plot_layout(design = layout)


# Plot 4 = unknown early in ALOHA -----
sample_name <- "Smp_ALOHA_crude"
sample_name_to_plot <- "ALOHA"
lipid_name_to_plot <- "unknown_mz1003.5009"
rt_lipid = 21.8
ESI_dat_file <- "RawDat/20201026_QE_secondLipidRun/201028_Smp_ALOHA_crude_highmass.mzXML"
MS2_dat_file <- "RawDat/20201026_QE_secondLipidRun/201028_Smp_ALOHA_crude_highmass_mergedMS2.rds"
mass = 1003.5009
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
  distance_left_join(MS2_library %>%
                       select(mz, Formula) %>%
                       unique(), max_dist = 0.02, by = "mz") %>%
  mutate(mz = mz.x) %>%
  mutate(has_As = !is.na(Formula)) %>%
  mutate(rt = MS2_time)
MS2s_top <- MS2s %>%
  mutate(rank = rank(-intensity)) %>%
  filter(((has_As == TRUE) & (rank < 20) )| (rank < 3))

p4a <- ggplot() +
  geom_line(data = icap_dat_sub, 
            aes(x = time/60, y = intensity_smoothed*1000), color = "coral1",
            alpha = 0.5)+
  geom_line(data = EIC_df, aes(x = rt_min+.1, y = intensity))+
  scale_x_continuous("Time (min)", limits = c(10, 32)) +
  scale_y_continuous("Intensity", limits = c(0, NA), expand = c(0, NA)) +
  labs(title = paste0(lipid_name_to_plot, " in ", sample_name_to_plot))+
  my_theme+
  theme(axis.title.x = element_blank())

p4b <- ggplot(data = MS2s, aes(x = mz,xend = mz, 
                               y = 0, yend = intensity)) +
  geom_segment(data = MS2s %>% 
                 filter(has_As == TRUE), 
               color = 'coral1', 
               size = 1)+
  geom_segment()+
  geom_hline(yintercept = 0)+
  geom_text(data = MS2s_top, 
            aes(x = mz, y = intensity+5000, 
                label = round(mz, digits = 4)), 
            size = 2, color = 'black',
            check_overlap = TRUE)+
  scale_x_continuous("m/z", limits = c(0, 450),
                     expand = c(0,0))+
  scale_y_continuous(element_blank(),  
                     limits=c(0, max(MS2s$intensity*1.05)), expand = c(0,NA))+  
  my_theme+
  theme(axis.title.x = element_blank())

p4 <- p4a + p4b + 
  plot_layout(design = layout)



# Plot 5 = unknown late in PS2-----
sample_name <- "Smp_PS2_crude"
sample_name_to_plot <- "ETNP-PS2"
lipid_name_to_plot <- "unknown_mz999.5884"
rt_lipid = 27.6
ESI_dat_file <- "RawDat/20201026_QE_secondLipidRun/201028_Smp_PS2_elute_highmass.mzXML"
MS2_dat_file <- "RawDat/20201026_QE_secondLipidRun/201028_Smp_PS2_elute_highmass_mergedMS2.rds"
mass = 999.5872
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
  distance_left_join(MS2_library %>%
                       select(mz, Formula) %>%
                       unique(), max_dist = 0.02, by = "mz") %>%
  mutate(mz = mz.x) %>%
  mutate(has_As = !is.na(Formula)) %>%
  mutate(rt = MS2_time)
MS2s_top <- MS2s %>%
  mutate(rank = rank(-intensity)) %>%
  filter(((has_As == TRUE) & (rank < 20) )| (rank < 3))

p5a <- ggplot() +
  geom_line(data = icap_dat_sub, 
            aes(x = time/60, y = intensity_smoothed*1000), color = "coral1",
            alpha = 0.5)+
  geom_line(data = EIC_df, aes(x = rt_min+.1, y = intensity))+
  scale_x_continuous("Time (min)", limits = c(10, 32)) +
  scale_y_continuous("Intensity", limits = c(0, NA), expand = c(0, NA)) +
  labs(title = paste0(lipid_name_to_plot, " in ", sample_name_to_plot))+
  my_theme

p5b <- ggplot(data = MS2s, aes(x = mz,xend = mz, 
                               y = 0, yend = intensity)) +
  geom_segment(data = MS2s %>% 
                 filter(has_As == TRUE), 
               color = 'coral1', 
               size = 1)+
  geom_segment()+
  geom_hline(yintercept = 0)+
  geom_text(data = MS2s_top, 
            aes(x = mz, y = intensity+5000, 
                label = round(mz, digits = 4)), 
            size = 2, color = 'black',
            check_overlap = TRUE)+
  scale_x_continuous("m/z", expand = c(0,NA), limits = c(0,450))+
  scale_y_continuous(element_blank(),  
                     limits=c(0, max(MS2s$intensity*1.05)), expand = c(0,NA))+  
  my_theme

p5 <- p5a + p5b + 
  plot_layout(design = layout)


# Put them all together ----
p6 <- p1/p2/p3/p4/p5
p6

save_plot("Figures/ManuscriptReady/Example_lipids_chromats_and_spectra.pdf", p6, 
          base_width = 6, base_height = 8)
