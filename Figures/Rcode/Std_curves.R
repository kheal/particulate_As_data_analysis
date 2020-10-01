library(here)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(broom)
library(wesanderson)

#TO DO: Check the slopes etc on the different sig_types from the %HCl
#TO DO: Check units - is this really ppb?  Or is it nMol

dat.file <- "Intermediates/Std_curves_variableHClandMeOH.csv"

#Read in data file of cps of arsenic standard curves----
dat <- read_csv(dat.file) 

#Plots exploring the relationship between 75As concentration, 75As signal, and %HCl-----
dat.HCl <- dat %>%
  filter(`%MeOH` == 2.5)

dat.noHCl <- dat.HCl %>%
  filter(`%HCl` == 0) %>%
  mutate(`75As_noHCl` = `75As`) %>%
  select(`[As]_ppb`, `75As_noHCl`, Mode)%>%
  group_by(`[As]_ppb`, Mode) %>%
  summarise(`75As_noHCl` = mean(`75As_noHCl`)) %>% ungroup()

#Get 77Se corrected 75As signals
dat.77correct <- dat.HCl %>%
  mutate(sig_type = "77Correction") %>%
  mutate(Signal = `75As`-`77Se`*(0.7576/.2424)) %>% 
  select(`[As]_ppb`, Mode, `%HCl`, sig_type, Signal)

#Get raw 75As signals
dat.raw <- dat.HCl %>%
  mutate(sig_type = "raw") %>%
  mutate(Signal = `75As`) %>% 
  select(`[As]_ppb`, Mode, `%HCl`, sig_type, Signal)

#Get blank corrected 75As signals
dat.HCl.blanks <- dat %>%
  filter(`%MeOH` == 2.5) %>% filter(`[As]_ppb` == 0) %>% 
  group_by(Mode, `%HCl`) %>%
  summarise(AsBlank = mean(`75As`))

dat.blankcorrect <- dat.HCl %>%
  left_join(dat.HCl.blanks, by = c("Mode", "%HCl")) %>%
  mutate(sig_type = "blankCorrected") %>%
  mutate(Signal = `75As`-AsBlank ) %>%
  select(`[As]_ppb`, Mode, `%HCl`, sig_type, Signal) 

#Connect
dat.curves.all <- rbind(dat.raw, dat.77correct) %>%
  rbind(dat.blankcorrect) 

dat.curves.all <- dat.curves.all %>%
  mutate(sig_type = factor(dat.curves.all$sig_type, levels = c("raw", "77Correction", "blankCorrected"))) %>%
  mutate(Mode = factor(dat.curves.all$Mode, levels = c("Std", "KED")))

#Get nice names for facets
sigtype.labs <- c("uncorrected \n 75As", "corrected by \n77Se signal", "corrected by \nblank subtraction")
names(sigtype.labs) <- c("raw", "77Correction",
                         "blankCorrected")


#Plot it up!
g.HCl.curves.log <- ggplot(dat.curves.all%>%
                             filter(Signal > 0) %>%
                             filter(`[As]_ppb` > 0.005), aes(x =  `[As]_ppb`, y = Signal, color = factor(`%HCl`))) +
  stat_smooth(geom='line', alpha=0.5, se=FALSE, method = "lm", size = 1.5)+
  geom_point(size = 2)+
  scale_color_manual(name="% HCl",
                     values = pal)+
  facet_grid(rows = vars(Mode), cols = vars(sig_type), scales = "free", 
             labeller = labeller(sig_type = sigtype.labs))+
  scale_x_log10()+scale_y_log10()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6))

g.HCl.curves.log


#Plots exploring the relationship between 75As concentration, 75As signal, and %MeOH-----
dat.MeOH <- dat %>%
  filter(`%HCl` == 0)

#Get 77Se corrected 75As signals
dat.77correct.MeOH <- dat.MeOH %>%
  mutate(sig_type = "77Correction") %>%
  mutate(Signal = `75As`-`77Se`*(0.7576/.2424)) %>% 
  select(`[As]_ppb`, Mode, `%MeOH`, sig_type, Signal)

#Get raw 75As signals
dat.raw.MeOH <- dat.MeOH %>%
  mutate(sig_type = "raw") %>%
  mutate(Signal = `75As`) %>% 
  select(`[As]_ppb`, Mode, `%MeOH`, sig_type, Signal)

#Get blank corrected 75As signals
dat.MeOH.blanks <- dat %>%
  filter(`%HCl` == 0) %>% filter(`[As]_ppb` == 0) %>% 
  group_by(Mode, `%MeOH`) %>%
  summarise(AsBlank = mean(`75As`))

dat.blankcorrect.MeOH <- dat.MeOH %>%
  left_join(dat.MeOH.blanks, by = c("Mode", "%MeOH")) %>%
  mutate(sig_type = "blankCorrected") %>%
  mutate(Signal = `75As`-AsBlank ) %>%
  select(`[As]_ppb`, Mode, `%MeOH`, sig_type, Signal)

#Connect
dat.curves.all.MeOH<- rbind(dat.raw.MeOH, dat.77correct.MeOH) %>%
  rbind(dat.blankcorrect.MeOH) 

dat.curves.all.MeOH <- dat.curves.all.MeOH %>%
  mutate(sig_type = factor(dat.curves.all.MeOH$sig_type, levels = c("raw", "77Correction", "blankCorrected"))) %>%
  mutate(Mode = factor(dat.curves.all.MeOH$Mode, levels = c("Std", "KED")))

#Plot it up!-------
g.MeOH.curves.log <- ggplot(dat.curves.all.MeOH%>%
                             filter(Signal > 0) %>%
                             filter(`[As]_ppb` > 0.005), aes(x =  `[As]_ppb`, y = Signal, color = factor(`%MeOH`))) +
  stat_smooth(geom='line', alpha=0.5, se=FALSE, method = "lm", size = 1.5)+
  geom_point(size = 2)+
  scale_color_manual(name="% MeOH",
                     values = c("grey20", pal))+
  facet_grid(rows = vars(Mode), cols = vars(sig_type), scales = "free", 
             labeller = labeller(sig_type = sigtype.labs))+
  scale_x_log10()+scale_y_log10()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6))

g.MeOH.curves.log


#Get summary statistics - Mode, correction type, %HCl, R2, slope, error on slope
dat.HCl.stats <- dat.curves.all %>%
  ungroup() %>%
  rename(percentHCl = `%HCl`) %>%
  group_by(Mode, sig_type, percentHCl) %>%
  do(fitSig = lm(`[As]_ppb`~Signal, .))

dfHClCoef <- tidy(dat.HCl.stats, fitSig) %>%
  filter(term == "Signal") %>% 
  select(Mode, sig_type, percentHCl, estimate, std.error) %>%
  rename(slope = estimate, slope.std.error = std.error)

dfHClSumm <- glance(dat.HCl.stats, fitSig) %>%
  select(Mode, sig_type, percentHCl, r.squared) 
  
dat.HCl.summary.stats <- dfHClCoef %>% 
  left_join(dfHClSumm, 
            by = c("Mode", "percentHCl", "sig_type"))
  

#Save out plots
save_plot("Figures/ManuscriptReady/HCl_curves_loglog.pdf", g.HCl.curves.log, device = "pdf", base_width =  5, base_height = 3, units = "in")


save_plot("Figures/ManuscriptReady/MeOH_curves_loglog.pdf", g.MeOH.curves.log, device = "pdf", base_width =  5, base_height = 3, units = "in")

  