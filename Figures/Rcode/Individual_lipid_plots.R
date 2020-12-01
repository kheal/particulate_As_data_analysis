library(tidyverse)
library(here)
library(cowplot)
theme_set(theme_cowplot())
library(stats)

# Name files ----
found_icap_peaks_filename <- "Intermediates/AsLipids_ICPpeakareas.csv"
sample_matcher_filename <- "MetaData/QE_ICAP_samplematcher.csv"
lipids_to_plot_filename <- "Intermediates/Targeted_firstPass.csv"
as_db_high_filename <- "MetaData/AsLipidDatabase/highmass_inclusionlist.csv"
as_db_low_filename <- "MetaData/AsLipidDatabase/lowmass_inclusionlist.csv"
MS2_library_filename <- "MetaData/AsLipidDatabase/MS2_Fragments.csv"

# Load files
sample_matcher <- read_csv(sample_matcher_filename) %>%
#  filter(str_detect(ICP_sampID, "crude")) %>%
  filter(!str_detect(ICP_sampID, "KM|lank|filter"))
peaks_all <- read_csv(found_icap_peaks_filename)
as_bd <- read_csv(as_db_high_filename) %>%
  bind_rows(read_csv(as_db_low_filename))
location_of_QEfiles <- "RawDat/20201026_QE_secondLipidRun"
lipids_to_plot <- read_csv(lipids_to_plot_filename, comment = "#")
MS2_library <- read_csv(MS2_library_filename)%>%
  mutate(mz = sapply(Formula, 
                     function(x)getMolecule(x)$exactmass))



# Get ICP data for crude injection of each sample
ICP_data <- list()
ICP_crude_samps <- sample_matcher %>% filter(str_detect(ICP_sampID, "crude")) 
for (i in 1:5){
  ICP_sampID <- unique(ICP_crude_samps$ICP_sampID)[i]
  ICP_data[[i]] <- read_csv(paste0("Intermediates/Baseline_CSVs/", 
                                   ICP_sampID, "_ICPwithbaseline.csv"))}

# Get MS1 data for each sample and MS2 
ESI_dat_highmass <- list()
MS2_dat_highmass <- list()
ESI_elute_samps <- sample_matcher %>% 
  filter(str_detect(QE_filename, "elute") & str_detect(QE_filename, "highmas")) 
for (i in 1:5){
  QE_file <- paste0(location_of_QEfiles, "/", ESI_elute_samps$QE_filename[i])
  ESI_dat_highmass[[i]] <- xcmsRaw(QE_file,profstep=0.01, profmethod="bin",
                        profparam=list(),
                        includeMSn=FALSE, mslevel=NULL, scanrange=NULL)
  MS2_dat_highmass[[i]] <- readRDS(str_replace(QE_file, ".mzXML", '_mergedMS2.rds'))
}

# Enter in name of lipid and pull out EIC from each sample and plot faceted
lipids_to_plot_withmass <- lipids_to_plot %>%
  mutate(Lipid_Name = Lipid_to_replot) %>% 
  select(Lipid_Name, Time) %>% 
  left_join(as_bd, by = "Lipid_Name") %>%
  filter(mz > 900)


# Get shorter names of each sample to plot
names <- ESI_elute_samps %>% select(QE_filename) %>%
  separate(QE_filename, into = c("Date", "Type", "SampID_short"), sep = "_") %>%
  mutate(sample = SampID_short)

for (j in 1:length(lipids_to_plot_withmass$mz)){

EICs <- list()
for (i in 1:5){
  mzrange<-c(-0.01,0.01)+lipids_to_plot_withmass$mz[j]
  EIC<-rawEIC(ESI_dat_highmass[[i]],mzrange)
  times<-ESI_dat_highmass[[i]]@scantime/60
  EIC_df<-data.frame(cbind(times,unlist(EIC[[2]])))
  colnames(EIC_df) <- c("times", "int")
  EIC_df$sample <- names$SampID_short[i]
  EICs[[i]] <- EIC_df
}
EICs_unpack <- do.call(rbind, EICs)

MS2s <- list()
for (i in 1:5){
  MS2s_mz <- c()
  MS2s_rt <- getSpectrum(MS2_dat_highmass[[i]], "rt", 
                         lipids_to_plot_withmass$Time[j]*60, 
                         rt.tol = 20)
   MS2s_mz <- getSpectrum(MS2s_rt,  "precursor", 
                          lipids_to_plot_withmass$mz[j], 
                          mz.tol = .1)
   if (length(MS2s_mz)==0){MS2s[[i]] = tibble(mz = NA, intensity = 0, has_As = FALSE) %>%
     mutate(sample = names$SampID_short[i])} else{
     if (length(MS2s_mz)>1){MS2s_mz <- MS2s_mz[[1]]} 
     MS2_time <- MS2s_mz@rt/60
     MS2_precursormass <- MS2s_mz@precursor
     MS2s[[i]] <- tibble(mz = as.numeric(MS2s_mz@spectrum[,1]),
                         intensity = as.numeric(MS2s_mz@spectrum[,2])) %>%
     distance_left_join(MS2_library %>% select(mz, Formula) %>% unique(), max_dist = 0.02) %>%
     mutate(mz = mz.x) %>%
     mutate(has_As = !is.na(Formula)) %>%
     mutate(sample = names$SampID_short[i])}
   
   
   }
MS2s_unpack <- do.call(bind_rows, MS2s)


# EICs of the samples
my_theme <- theme(strip.background = element_blank(),
                  plot.title = element_text(size = 9),
                  axis.text.y = element_blank(),
                  axis.title = element_text(size = 8),
                  axis.text.x = element_text(size = 7),
                  legend.position = "none")

g <- ggplot(data = EICs_unpack, aes(x = times, y = int)) +
  geom_line()+
  facet_wrap(~  sample, ncol = 1, scales = "free") +
  labs(title = paste0(lipids_to_plot_withmass$Lipid_Name[j], ", m/z = ",
                      round(lipids_to_plot_withmass$mz[j], digits = 4)))+
  scale_x_continuous("Retention time (minutes)", limits = c(7, 35),
                     expand = c(0,0))+
  scale_y_continuous("Intensity",
                     expand = c(0,0))+
  my_theme +
  theme(strip.text = element_text(hjust = 0, size = 8))
g

g2 <- ggplot(data = MS2s_unpack, aes(x = mz,xend = mz, y = 0, 
                                     yend = intensity, color = has_As)) +
  geom_segment()+
  facet_wrap(~  sample, ncol = 1, scales = "free")+
  scale_x_continuous("m/z", limits = c(0, 800),
                     expand = c(0,0))+
  scale_y_continuous("Intensity",
                     expand = c(0,0))+  
  scale_color_manual(values = c("black", "coral1"))+
  my_theme +
  theme(strip.text = element_blank())
g2

g3 <- plot_grid(g, g2, align = 'h', ncol =2)
g3

savefilename <- paste0("Figures/Preliminary/Each_individual_lipid/",
                       lipids_to_plot_withmass$Lipid_Name[j],
                       "_", round(lipids_to_plot_withmass$mz[j], digits = 4), ".pdf")
save_plot(savefilename, g3, base_height = 8, base_width = 6, units = "in")
}



