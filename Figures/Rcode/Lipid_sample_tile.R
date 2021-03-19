library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(nationalparkcolors)
library(ggtext)


# Name your files ----
dat_filename <-"Intermediates/Quantified_Individual_Lipids_betterint.csv"
dat_noicp_filename <- "MetaData/Targeted_visually_IDd_2.csv"
as_db_high_filename <- "MetaData/AsLipidDatabase/highmass_inclusionlist.csv"
as_db_low_filename <- "MetaData/AsLipidDatabase/lowmass_inclusionlist.csv"
as_bd <- read_csv(as_db_high_filename) %>%
  bind_rows(read_csv(as_db_low_filename))

# Set your theme
my_theme <- theme(strip.background = element_blank(),
                  plot.title = element_text(size = 9, hjust = 0.5),
                  axis.title = element_text(size = 8),
                  axis.text.x = element_text(size = 7),
                  axis.text.y = element_text(size = 7),
                  legend.position = "none")


# Load dat used by all----
dat <- read_csv(dat_filename) %>%
  select(`Sample ID`, lipid_ID, rt_sec, pMolAs_enviro) %>%
  mutate(quan = TRUE)

dat2 <- read_csv(dat_noicp_filename) %>%
  filter(!is.na(notes)) %>%
  mutate(rt_sec = ifelse(is.na(rt_sec), rt_min*60, rt_sec)) %>%
  rename(`Sample ID`  = sampleID) %>%
  select(`Sample ID`, lipid_ID, rt_sec) %>%
  mutate(quan = FALSE)

dat3 <- bind_rows(dat, dat2) %>%
  left_join(as_bd %>% select(Lipid_Name, mz) %>% unique() %>% 
              rename(lipid_ID = Lipid_Name), by = "lipid_ID") %>%
  mutate(mz = ifelse(str_detect(lipid_ID, "mz"),
                     str_extract_all(lipid_ID, "\\d+.\\d+"), mz)) 

dat4 <- dat3 %>%
  group_by(lipid_ID) %>%
  summarise(rt_min_ave = mean(rt_sec, na.rm = TRUE)/60)

dat5 <- dat3 %>%
  left_join(dat4, by = "lipid_ID") %>%
  mutate(lipid_ID_forprint = ifelse(str_detect(lipid_ID, "unknown"),
                                    str_extract(lipid_ID, "unknown_."),
                                    lipid_ID)) %>%
  mutate(long_name = ifelse(!is.na(mz), 
                            paste0("**", lipid_ID_forprint, "**, *", 
                                   round(as.numeric(mz), digits = 4), "*, ", 
                            round(rt_min_ave, digits = 1), " min"), 
                            paste0("**",lipid_ID_forprint, "**, ", 
                                   round(rt_min_ave, digits = 1), " min")))

dat6 <- dat5 %>% arrange(desc(rt_min_ave)) %>% 
  select(long_name) %>% unique()

dat7 <- dat5 %>% 
  mutate(long_name = factor(long_name, 
                                    levels = dat6$long_name)) %>%
  mutate(lipid_type = ifelse(str_detect(lipid_ID, "HC"), "AsHC",
                             ifelse(str_detect(lipid_ID, "AsSugPL"), "AsSugPL",
                                    ifelse(str_detect(lipid_ID, "AsSugPeL"), "AsSugPeL",
                                           ifelse(str_detect(lipid_ID, "mz"), "unknown, but m/z known",
                                                  "unknown")))))


dat7$lipid_type = factor(dat7$lipid_type, levels = c("AsHC", "AsSugPL", "AsSugPeL", "unknown, but m/z known", "unknown"))

# Plot up the data, these are in fMol/L
# To do: color by type of lipid
# Order samples
# Order lipids by RT
pal <- c(park_palette("Badlands", 4), 'grey80')
dat7$`Sample ID` <- str_wrap(dat7$`Sample ID`, width = 5) 
dat7 <- dat7 %>%
  filter(!(`Sample ID` == "ETNP-\nPS2" &
             lipid_ID == "AsSugPL982" &
             is.na(pMolAs_enviro)))

g <- ggplot(data = dat7, aes(x =`Sample ID`,  y = long_name, 
                             fill = lipid_type, 
                             label = round(pMolAs_enviro*1000, digits = 0))) +
  geom_tile(alpha = 0.5, color = "grey30") +
  geom_text(size = 2) +
  scale_fill_manual(values=pal)+
  my_theme+
  theme(axis.title = element_blank(), 
        axis.text.y = ggtext::element_markdown())
g


save_plot("Figures/ManuscriptReady/Results_tiles.pdf", g, base_height = 4.5, base_width = 4)
